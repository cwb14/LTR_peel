#!/usr/bin/env python3
"""Family-level bootstrap statistics for the non-autonomous benchmark (simulations).

Unit of independence = PrinTE library family (copies of one family share sequence). 1,000 resamples (numpy
default_rng(1)) of families with replacement, within each class; every copy of a drawn family is carried along.
  sim_gap.tsv  per set (S1, S2; 3 seeds pooled) and tool: recall of coding and protein-less elements, the gap
               (protein-less - coding), and its split into a length part and a residual. Length part = coding recall
               reweighted to the protein-less length distribution (7 length bins) minus coding recall; residual part =
               protein-less recall minus that reweighted coding recall. Bins with no coding element are left out of
               the reweighting (their protein-less weight is dropped and the rest renormalized).
  div_d50.tsv  divergence series (3 seeds pooled): per tool and class, d50 = LTR K2P at which pooled recall falls to
               half its young-genome value (div_score.d50 interpolation on the class's median K2P per timepoint),
               with family-bootstrap 95% CI, and the difference protein-less - coding.
"""
import csv
import os
import sys
from collections import defaultdict

import numpy as np

D = os.path.dirname(os.path.abspath(__file__))
TOOLS = ["ltr_peel", "ltrharvest", "ltr_finder", "ltr_retriever", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]
LB = [0, 2000, 3000, 4000, 5000, 6000, 8000, 10 ** 9]
NB = 1000


def ci(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return (np.percentile(x, 2.5), np.percentile(x, 97.5)) if len(x) else (np.nan, np.nan)


def f(x):
    return "NA" if x is None or not np.isfinite(x) else f"{x:.4f}"


def sim_gap():
    # per (set, tool, cls, family, lenbin): [n, hit]
    agg = defaultdict(lambda: np.zeros((len(LB) - 1, 2)))
    for r in csv.DictReader(open(f"{D}/sim_elements.tsv"), delimiter="\t"):
        if r["tool"] not in TOOLS or r["cls"] == "partial":
            continue
        b = np.searchsorted(LB, int(r["len"]), side="right") - 1
        a = agg[(r["ds"][:2], r["tool"], r["cls"], r["family"])]
        a[b, 0] += 1
        a[b, 1] += int(r["hit"])
    rng = np.random.default_rng(1)
    out = []
    for st in ("S1", "S2"):
        fams = {c: sorted({k[3] for k in agg if k[0] == st and k[2] == c}) for c in ("coding", "proteinless")}
        draws = {c: rng.integers(0, len(v), size=(NB, len(v))) for c, v in fams.items()}
        for tool in TOOLS:
            M = {c: np.stack([agg[(st, tool, c, fm)] for fm in fams[c]]) for c in fams}   # fam x bin x (n, hit)

            def est(Mc, Mp):
                nc, hc = Mc[..., 0].sum(0), Mc[..., 1].sum(0)
                npb, hp = Mp[..., 0].sum(0), Mp[..., 1].sum(0)
                rc, rp = hc.sum() / nc.sum(), hp.sum() / npb.sum()
                ok = nc > 0
                w = npb * ok
                rcb = np.divide(hc, nc, out=np.zeros_like(hc), where=ok)
                rcm = (w * rcb).sum() / w.sum()
                return rc, rp, rp - rc, rcm - rc, rp - rcm, rcm
            pt = est(M["coding"], M["proteinless"])
            bs = np.array([est(M["coding"][draws["coding"][i]], M["proteinless"][draws["proteinless"][i]])
                           for i in range(NB)])
            row = dict(set=st, tool=tool, n_fam_coding=len(fams["coding"]), n_fam_proteinless=len(fams["proteinless"]),
                       n_coding=int(M["coding"][..., 0].sum()), n_proteinless=int(M["proteinless"][..., 0].sum()))
            for j, nm in enumerate(["recall_coding", "recall_proteinless", "gap", "length_part", "residual_part",
                                    "recall_coding_lenmatched"]):
                lo, hi = ci(bs[:, j])
                row[nm], row[nm + "_lo"], row[nm + "_hi"] = f(pt[j]), f(lo), f(hi)
            out.append(row)
    return out


def d50(ks, rs):
    if not rs[0] > 0:
        return np.nan
    half = rs[0] / 2
    for i in range(1, len(rs)):
        if rs[i] <= half:
            return ks[i - 1] + (ks[i] - ks[i - 1]) * (rs[i - 1] - half) / (rs[i - 1] - rs[i]) if rs[i - 1] != rs[i] \
                else ks[i]
    return np.inf          # never fell to half within the series


def div_d50():
    # per (tool, cls, family): array timepoint x (n, hit); class median K2P per timepoint from div_class.tsv
    agg = defaultdict(lambda: np.zeros((11, 2)))
    for r in csv.DictReader(open(f"{D}/div_elements.tsv"), delimiter="\t"):
        if r["tool"] not in TOOLS or r["cls"] == "partial":
            continue
        a = agg[(r["tool"], r["cls"], r["family"])]
        a[int(r["timepoint"]), 0] += 1
        a[int(r["timepoint"]), 1] += int(r["hit"])
    med = defaultdict(list)
    for r in csv.DictReader(open(f"{D}/div_class.tsv"), delimiter="\t"):
        if r["tool"] == "ltr_peel":
            med[(r["cls"], int(r["timepoint"]))].append(float(r["cls_median_k2p"]))
    K = {c: np.array([np.median(med[(c, t)]) for t in range(11)]) for c in ("coding", "proteinless")}
    rng = np.random.default_rng(1)
    fams = {c: sorted({k[2] for k in agg if k[1] == c}) for c in ("coding", "proteinless")}
    draws = {c: rng.integers(0, len(v), size=(NB, len(v))) for c, v in fams.items()}
    out = []
    for tool in TOOLS:
        res = {}
        for c in ("coding", "proteinless"):
            M = np.stack([agg[(tool, c, fm)] for fm in fams[c]])
            if M[..., 0].sum(0).min() == 0:
                res[c] = None
                continue
            r0 = M[..., 1].sum(0) / M[..., 0].sum(0)
            bs = np.array([d50(K[c], (lambda S: S[..., 1].sum(0) / np.maximum(S[..., 0].sum(0), 1))(M[draws[c][i]]))
                           for i in range(NB)])
            res[c] = (d50(K[c], r0), bs, r0)
        for c in ("coding", "proteinless"):
            if res[c] is None:
                out.append(dict(tool=tool, cls=c, n_fam=len(fams[c]), recall_t0="NA", d50="NA", d50_lo="NA",
                                d50_hi="NA", d50_diff_vs_coding="NA", diff_lo="NA", diff_hi="NA"))
                continue
            pt, bs, r0 = res[c]
            lo, hi = ci(bs[np.isfinite(bs)]) if np.isfinite(bs).sum() > 0.9 * NB else (np.nan, np.nan)
            row = dict(tool=tool, cls=c, n_fam=len(fams[c]), recall_t0=f(r0[0]),
                       d50=f(pt) if np.isfinite(pt) else (">0.37" if pt == np.inf else "NA"), d50_lo=f(lo), d50_hi=f(hi))
            if c == "proteinless" and res["coding"] is not None and r0[0] > 0:
                dif = bs - res["coding"][1]
                dlo, dhi = ci(dif[np.isfinite(dif)]) if np.isfinite(dif).sum() > 0.9 * NB else (np.nan, np.nan)
                row.update(d50_diff_vs_coding=f(pt - res["coding"][0]), diff_lo=f(dlo), diff_hi=f(dhi))
            else:
                row.update(d50_diff_vs_coding="NA", diff_lo="NA", diff_hi="NA")
            out.append(row)
    return out


def write(path, rows):
    with open(path, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    write(f"{D}/sim_gap.tsv", sim_gap())
    write(f"{D}/div_d50.tsv", div_d50())
    print("nonauto_stats: done", file=sys.stderr)
