#!/usr/bin/env python3
"""Recall of all elements vs LTR divergence (aged series SD1_00..SD3_10), three genomes pooled, family-bootstrap CI.

Elements: the primary set of approach 1b (std truth elements whose PrinTE library entry has an LTR pair), every
protein class together, one row per element x tool in experiments/nonauto/div_elements.tsv (nonauto_sim.do_div;
runs with status ok or no_candidates; no_candidates = 0 calls, recall 0).
Per tool x timepoint: recall = matched / elements over the three genomes. 95% CI: 1,000 resamples of PrinTE library
families with replacement (numpy default_rng(1)), every copy of a drawn family in all three genomes and timepoints
carried along, percentiles 2.5 / 97.5 (nonauto_stats.ci): the procedure of nonauto_stats.div_d50 with the class split
removed (one family pool; the same draws for every tool and timepoint).
x = median LTR-LTR K2P of the pooled elements at the timepoint (each element once); the median of the three per-seed
medians in report/tables/sim_divergence.tsv is written beside it.
Check (fatal): pooled n and matched equal the sums of n_truth and tp over the three seeds in sim_divergence.tsv, and
each tool has all three seeds at every timepoint.
Output (--out): div_all_ci.tsv.
"""
import argparse
import csv
import os
import statistics
import sys
from collections import defaultdict

import numpy as np

sys.dont_write_bytecode = True
B = os.environ["BENCH_ROOT"]
NA = f"{B}/experiments/nonauto"
sys.path.insert(0, NA)
from nonauto_stats import NB, ci  # noqa: E402

TOOLS = ["ltrquest", "ltr_retriever", "ltrharvest", "ltr_finder", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]
DIAG = ["ltrquest_hc"]
NT = 11


def die(msg):
    sys.exit(f"error: {msg}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)), help="output dir (default: script dir)")
    a = ap.parse_args(argv)
    agg = defaultdict(lambda: np.zeros((NT, 2)))          # (tool, family) -> timepoint x (n, hit)
    per_ds = defaultdict(lambda: [0, 0])                  # (ds, tool) -> n, hit
    seeds = defaultdict(set)                              # (tool, timepoint) -> seeds
    k2p = defaultdict(dict)                               # timepoint -> {(ds, element): k2p}
    for r in csv.DictReader(open(f"{NA}/div_elements.tsv"), delimiter="\t"):
        if r["tool"] not in TOOLS + DIAG:
            continue
        t = int(r["timepoint"])
        agg[(r["tool"], r["family"])][t] += (1, int(r["hit"]))
        per_ds[(r["ds"], r["tool"])][0] += 1
        per_ds[(r["ds"], r["tool"])][1] += int(r["hit"])
        seeds[(r["tool"], t)].add(r["seed"])
        k2p[t][(r["ds"], r["element"])] = float(r["k2p"])
    sd = list(csv.DictReader(open(f"{B}/report/tables/sim_divergence.tsv"), delimiter="\t"))
    bad = []
    tab = defaultdict(lambda: [0, 0])
    seed_med = defaultdict(list)
    for r in sd:
        if r["tool"] not in TOOLS + DIAG:
            continue
        t = int(r["timepoint"])
        tab[(r["tool"], t)][0] += int(r["n_truth"])
        tab[(r["tool"], t)][1] += int(r["tp"])
        if (r["ds"], r["tool"]) in per_ds and tuple(per_ds[(r["ds"], r["tool"])]) != (int(r["n_truth"]), int(r["tp"])):
            bad.append(f"{r['ds']} {r['tool']}: div_elements {per_ds[(r['ds'], r['tool'])]} vs sim_divergence "
                       f"{r['n_truth']}/{r['tp']}")
        if r["tool"] == "ltrquest":
            seed_med[t].append(float(r["median_k2p"]))
    for tool in TOOLS + DIAG:
        for t in range(NT):
            if seeds[(tool, t)] != {"1", "2", "3"}:
                bad.append(f"{tool} t{t}: seeds {sorted(seeds[(tool, t)])}")
    fams = sorted({f for (_, f) in agg})
    draws = np.random.default_rng(1).integers(0, len(fams), size=(NB, len(fams)))
    rows = []
    for tool in TOOLS + DIAG:
        M = np.stack([agg[(tool, f)] for f in fams])     # fam x timepoint x (n, hit)
        n, h = M[..., 0].sum(0), M[..., 1].sum(0)
        S = M[draws]                                      # draw x fam x timepoint x 2
        bs = S[..., 1].sum(1) / np.maximum(S[..., 0].sum(1), 1)
        for t in range(NT):
            if (int(n[t]), int(h[t])) != tuple(tab[(tool, t)]):
                bad.append(f"{tool} t{t}: pooled {int(n[t])}/{int(h[t])} vs sim_divergence sum {tab[(tool, t)]}")
            lo, hi = ci(bs[:, t])
            rows.append(dict(tool=tool, timepoint=t, median_k2p=statistics.median(k2p[t].values()),
                             median_k2p_seeds=statistics.median(seed_med[t]), n_truth=int(n[t]), tp=int(h[t]),
                             recall=h[t] / n[t], ci_lo=lo, ci_hi=hi, n_fam=len(fams), n_seeds=3,
                             role="diagnostic" if tool in DIAG else "primary"))
    if bad:
        die("checks failed:\n  " + "\n  ".join(bad[:40]))
    print(f"checks passed: per genome x tool and pooled n/tp equal sim_divergence.tsv ({len(rows)} tool x timepoint "
          f"rows); {len(fams)} families", file=sys.stderr)
    cols = ["tool", "timepoint", "median_k2p", "median_k2p_seeds", "n_truth", "tp", "recall", "ci_lo", "ci_hi", "n_fam",
            "n_seeds", "role"]
    with open(f"{a.out}/div_all_ci.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
