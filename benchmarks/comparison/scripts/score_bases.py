#!/usr/bin/env python3
"""Base-level LTR-RT annotation accuracy of a RepeatMasker/MaskerAide .out against PrinTE's exact truth.

Truth LTR bases = union of PrinTE features in categories LTR_intact, LTR_piece, LTR_frag, LTR_solo (printe_truth.py
OUT.features.bed, original names). Test LTR bases = union of .out hits whose class/family matches lib-test.py's LTR
pattern (RLG|RLC|RLB|RLR|RLE|LTR|RLX|Gypsy|Copia). N bases are excluded from the genome, as lib-test.py does.
Writes one TSV row: TP FP FN TN sens spec prec F1 (same definitions as lib-test.py).
"""
import argparse
import csv
import re
import sys
from collections import defaultdict

LTR = re.compile("RLG|RLC|RLB|RLR|RLE|LTR|RLX|Gypsy|Copia", re.I)
TRUTH_CATS = {"LTR_intact", "LTR_piece", "LTR_frag", "LTR_solo"}


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def inter(a, b):
    i = j = tot = 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if e > s:
            tot += e - s
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def total(iv):
    return sum(e - s for s, e in iv)


def read_out(path, n2o):
    hits = defaultdict(list)
    with open(path) as f:
        for line in f:
            c = line.split()
            if len(c) < 11 or not c[0].isdigit():
                continue
            if LTR.search(c[10]):
                chrom = n2o.get(c[4], c[4])
                hits[chrom].append((int(c[5]) - 1, int(c[6])))
    return {k: merge(v) for k, v in hits.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description="Base-level LTR accuracy of a masking .out vs PrinTE truth.")
    ap.add_argument("--out-file", required=True, help="RepeatMasker-format .out (normalized IDs)")
    ap.add_argument("--features", required=True, help="printe_truth.py OUT.features.bed (original names)")
    ap.add_argument("--names", required=True, help="data/<ds>/names.tsv (lengths, N counts, name map)")
    ap.add_argument("--label", default="", help="label column")
    ap.add_argument("--out", required=True, help="output TSV row")
    a = ap.parse_args(argv)
    n2o, size, nN = {}, {}, 0
    with open(a.names) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            n2o[r["new_id"]] = r["original_id"]
            size[r["original_id"]] = int(r["length"])
            nN += int(r["n_count"])
    truth = defaultdict(list)
    with open(a.features) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[3] in TRUTH_CATS:
                truth[c[0]].append((int(c[1]), int(c[2])))
    truth = {k: merge(v) for k, v in truth.items()}
    test = read_out(a.out_file, n2o)
    tp = sum(inter(truth.get(c, []), test.get(c, [])) for c in size)
    t_all = sum(total(v) for v in truth.values())
    p_all = sum(total(v) for v in test.values())
    fn, fp = t_all - tp, p_all - tp
    tn = sum(size.values()) - nN - tp - fn - fp
    sens = tp / (tp + fn) if tp + fn else float("nan")
    spec = tn / (tn + fp) if tn + fp else float("nan")
    prec = tp / (tp + fp) if tp + fp else float("nan")
    f1 = 2 * tp / (2 * tp + fp + fn) if tp + fp + fn else float("nan")
    with open(a.out, "w") as f:
        f.write("label\tTP\tFP\tFN\tTN\tsens\tspec\tprec\tF1\tmasked_mb\ttruth_mb\n")
        f.write(f"{a.label}\t{tp}\t{fp}\t{fn}\t{tn}\t{sens:.4f}\t{spec:.4f}\t{prec:.4f}\t{f1:.4f}\t"
                f"{p_all / 1e6:.2f}\t{t_all / 1e6:.2f}\n")
    print(f"score_bases: {a.label} sens={sens:.4f} prec={prec:.4f} F1={f1:.4f}", file=sys.stderr)


if __name__ == "__main__":
    main()
