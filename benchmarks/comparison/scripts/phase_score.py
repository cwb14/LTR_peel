#!/usr/bin/env python3
"""Score an ltrharvest_phasing10.py run against the authors' subgenome labels.

Truth: data/<ds>/truth.tsv (seq = normalized ID, label = A|B): P. annua from the authors' chromosome names (Pa<n>A/B);
A. suecica from whole-chromosome alignment to TAIR10 (A. thaliana-derived contigs match ~0.8 of aligned bases, A.
arenosa-derived ~0.2; this also labels the 3 arenosa chromosomes absent from the v2 B FASTA).
Metrics (one TSV row):
  n_sg            subgenomes inferred (truth 2)
  accuracy        share of truth sequences assigned to the subgenome that the majority of their label went to
                  (unassigned truth sequences count as wrong); correct = n_sg == 2 and accuracy == 1
  replicability   root split's half-library adjusted Rand (structure.tsv)
  strength        root split's within-group excess; margin = strength / strongest other candidate split
  oe_ratio_med    median over truth sequences of oe_own / oe_other (chromosomes.tsv); support_min
  n_ltr, n_clades from summary.json
"""
import argparse
import csv
import json
import math
import os
import statistics
import sys
from collections import Counter, defaultdict


def fai_names(path):
    with open(path) as f:
        return [l.split("\t")[0] for l in f if l.strip()]


def num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Score a phasing run against A/B subgenome labels.")
    ap.add_argument("rundir", help="phasing --outdir")
    ap.add_argument("--truth", required=True, help="truth TSV: seq (normalized ID), label (A|B)")
    ap.add_argument("--label", default="", help="label column value")
    ap.add_argument("--out", required=True, help="output TSV (one row; appended if it exists)")
    ap.add_argument("--allow-failed", action="store_true", help="if the run wrote no summary.json, record a failed row "
                    "(n_sg 0, correct 0) instead of stopping")
    a = ap.parse_args(argv)
    truth = {r["seq"]: r["label"] for r in csv.DictReader(open(a.truth), delimiter="\t")}
    if not truth:
        sys.exit(f"error: empty truth table {a.truth}")
    if a.allow_failed and not os.path.exists(os.path.join(a.rundir, "summary.json")):
        row = dict(label=a.label, n_sg=0, verdict="failed", n_truth=len(truth), n_assigned=0, accuracy=0.0, correct=0,
                   replicability=float("nan"), strength=float("nan"), margin=float("nan"), oe_ratio_med=float("nan"),
                   support_min=float("nan"), n_ltr=float("nan"), n_clades=float("nan"))
        new = not os.path.exists(a.out)
        with open(a.out, "a") as f:
            w = csv.DictWriter(f, fieldnames=list(row), delimiter="\t")
            if new:
                w.writeheader()
            w.writerow(row)
        print(f"phase_score: {a.label} failed run recorded", file=sys.stderr)
        return
    chrom = {r["sequence"]: r for r in csv.DictReader(open(os.path.join(a.rundir, "chromosomes.tsv")), delimiter="\t")}
    summ = json.load(open(os.path.join(a.rundir, "summary.json")))
    struct = list(csv.DictReader(open(os.path.join(a.rundir, "structure.tsv")), delimiter="\t"))
    sg = {s: chrom[s]["subgenome"] for s in truth if s in chrom and chrom[s].get("status") == "assigned"}
    # map each inferred subgenome to the truth label most of its truth sequences carry
    votes = defaultdict(Counter)
    for s, g in sg.items():
        votes[g][truth[s]] += 1
    g2lab = {g: c.most_common(1)[0][0] for g, c in votes.items()}
    correct_seq = sum(1 for s in truth if s in sg and g2lab.get(sg[s]) == truth[s])
    acc = correct_seq / len(truth)
    labs_used = Counter(g2lab.values())
    n_sg = int(summ.get("n_subgenomes", 0) or 0)
    root = next((r for r in struct if r["level"] == "root"), None)
    others = [num(r["strength"]) for r in struct if r["level"] != "root"]
    others = [x for x in others if math.isfinite(x)]
    strength = num(root["strength"]) if root else float("nan")
    rep = num(root["replicability"]) if root else float("nan")
    margin = strength / max(others) if others and max(others) > 0 else float("nan")
    ratios = [num(chrom[s]["oe_own"]) / num(chrom[s]["oe_other"]) for s in truth
              if s in chrom and num(chrom[s].get("oe_other")) > 0 and math.isfinite(num(chrom[s].get("oe_own")))]
    sup = [num(chrom[s]["support"]) for s in truth if s in chrom and math.isfinite(num(chrom[s].get("support")))]
    row = dict(label=a.label, n_sg=n_sg, verdict=summ.get("verdict", ""), n_truth=len(truth), n_assigned=len(sg),
               accuracy=round(acc, 4), correct=int(n_sg == 2 and acc == 1.0 and len(labs_used) == 2),
               replicability=round(rep, 4), strength=round(strength, 4), margin=round(margin, 3),
               oe_ratio_med=round(statistics.median(ratios), 3) if ratios else float("nan"),
               support_min=min(sup) if sup else float("nan"), n_ltr=summ.get("n_ltr"),
               n_clades=summ.get("n_clades"))
    new = not os.path.exists(a.out)
    with open(a.out, "a") as f:
        w = csv.DictWriter(f, fieldnames=list(row), delimiter="\t")
        if new:
            w.writeheader()
        w.writerow(row)
    print("\t".join(f"{k}={v}" for k, v in row.items()), file=sys.stderr)


if __name__ == "__main__":
    main()
