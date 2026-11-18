#!/usr/bin/env python3
"""Robustness checks of approach 4's library age-bias result (post hoc diagnostics, not a primary metric).

Pre-registered test (approach 4): loci holding an LTR-RT protein coding unit, pooled over the 9 primary
callers (union-find at reciprocal overlap >= 0.9); a locus is library-supported when >= half of its member calls have
>= 80% of their span covered by curated LTR hits; supported share by K2P bin; logistic slope on log10(K2P + 0.001).
Checks (each answers one way the decline could be an artifact):
  V0 pre-registered pool (reproduces report/tables/library_age_bias.tsv)
  V1 precise boundaries only: member calls with an exact 4-6 bp TSD and TG..CA ends, so the span is one whole
     insertion and the 80% rule cannot fail through overextended ends
  V2 strict callers only (LTR_retriever, HiTE, LTR_Finder, DANTE_LTR: structure- or domain-filtered call sets)
  V3 without LTRquest (the tool under test)
  V4 per caller: its own coding calls (no pooling), supported share by age bin
  V5 unsupported coding loci by curated LTR coverage (max over members): 'none' < 10%, 'partial' 10-80%
     -> family absent from the library vs copy only partly annotated
usage: age_bias_checks.py [-o experiments/a4_diagnostics] [-v]
"""
import argparse
import math
import os
import statistics
import sys
from collections import defaultdict

import numpy as np

B = os.environ["BENCH_ROOT"]
sys.path.insert(0, f"{B}/scripts")
import adj_summary as su  # noqa: E402
import proxies  # noqa: E402

TOOLS = ["ltrquest", "ltr_retriever", "ltrharvest", "ltr_finder", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]
STRICT = ["ltr_retriever", "hite", "ltr_finder", "dante_ltr"]
BINS = [("<0.01", 0, 0.01), ("0.01-0.03", 0.01, 0.03), ("0.03-0.1", 0.03, 0.1), (">=0.1", 0.1, math.inf)]


def k2p_bin(k):
    if k is None:
        return "NA"
    for name, lo, hi in BINS:
        if lo <= k < hi:
            return name
    return "NA"


def pooled(rb, tools, keep=lambda r: True):
    """coding-unit loci from the given callers' calls passing keep(); -> list of (supported, k2p, members)"""
    mem = [(r["chrom"], r["start"], r["end"], (t, k)) for t in tools for k, r in enumerate(rb.get(t, []))
           if r["coding"] == "1" and keep(r)]
    lab, _ = su.loci(mem)
    by = defaultdict(list)
    for c, s, e, (t, k) in mem:
        by[lab[(t, k)]].append(rb[t][k])
    out = []
    for rs in by.values():
        sup = sum(r["supported"] == "1" for r in rs) >= len(rs) / 2
        ks = [float(r["k2p"]) for r in rs if r["k2p"] != "NA"]
        out.append((sup, statistics.median(ks) if ks else None, rs))
    return out


def summarize(items):
    """items: [(supported, k2p)] -> per-bin rows and logistic slope"""
    rows = []
    for name, lo, hi in BINS:
        sel = [s for s, k in items if k is not None and lo <= k < hi]
        n, kk = len(sel), sum(sel)
        lo_, hi_ = su.wilson(kk, n)
        rows.append((name, n, kk / n if n else math.nan, lo_, hi_))
    xy = [(math.log10(k + 0.001), float(s)) for s, k in items if k is not None]
    if len(xy) >= 20 and 0 < sum(y for _, y in xy) < len(xy):
        slope, ci = su.logistic(np.array([x for x, _ in xy]), np.array([y for _, y in xy]))
    else:
        slope, ci = math.nan, (math.nan, math.nan)
    return rows, slope, ci


def main(argv=None):
    ap = argparse.ArgumentParser(description="Robustness checks of the approach-4 library age-bias result.")
    ap.add_argument("-o", "--outdir", default=f"{B}/experiments/a4_diagnostics", help="output directory")
    ap.add_argument("-v", "--verbose", action="store_true", help="print progress")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    precise = lambda r: r["tsd"] == "1" and r["tgca"] == "1"  # noqa: E731
    out_bins = [("ds", "check", "age_bin", "n", "supported_share", "wilson_lo", "wilson_hi", "logit_slope",
                 "slope_lo", "slope_hi")]
    out_cov = [("ds", "age_bin", "n_unsupported_coding_loci", "cov_none_lt10", "cov_partial_10_80",
                "median_max_cov")]
    for ds in ["arab", "rice", "finch", "maize", "gorilla"]:
        rb = su.load_ds(f"{B}/eval/adj/{ds}", TOOLS)
        if len(rb) < 5:
            continue
        pool = [t for t in TOOLS if t in rb]
        checks = {
            "V0_preregistered": pooled(rb, pool),
            "V1_precise_boundaries": pooled(rb, pool, precise),
            "V2_strict_callers": pooled(rb, [t for t in STRICT if t in rb]),
            "V3_without_ltrquest": pooled(rb, [t for t in pool if t != "ltrquest"]),
        }
        for name, loci in checks.items():
            rows, slope, ci = summarize([(s, k) for s, k, _ in loci])
            for r in rows:
                out_bins.append((ds, name) + r + (slope, ci[0], ci[1]))
        for t in pool:                                         # V4: each caller's own coding calls
            items = [(r["supported"] == "1", float(r["k2p"]) if r["k2p"] != "NA" else None)
                     for r in rb[t] if r["coding"] == "1"]
            rows, slope, ci = summarize(items)
            for r in rows:
                out_bins.append((ds, f"V4_caller_{t}") + r + (slope, ci[0], ci[1]))
        # V5: curated LTR coverage of unsupported coding loci (max over member calls)
        merged = proxies.ref_ltr(f"{B}/mask/{ds}/curated/{ds}.fa.out")
        starts = {c: [s for s, _ in v] for c, v in merged.items()}
        cov_by_bin = defaultdict(list)
        for sup, k, rs in checks["V0_preregistered"]:
            if sup:
                continue
            m = max(proxies.covered(merged, starts, r["chrom"], r["start"], r["end"]) / max(1, r["end"] - r["start"])
                    for r in rs)
            cov_by_bin[k2p_bin(k)].append(m)
        for name, _, _ in BINS + [("NA", 0, 0)]:
            v = cov_by_bin.get(name, [])
            if v:
                out_cov.append((ds, name, len(v), sum(x < 0.1 for x in v) / len(v),
                                sum(0.1 <= x < 0.8 for x in v) / len(v), statistics.median(v)))
        if a.verbose:
            print(f"age_bias_checks: {ds} done ({len(pool)} callers)", file=sys.stderr)

    def fmt(x):
        return f"{x:.4f}" if isinstance(x, float) else str(x)

    for fn, rows in (("age_bias_checks.tsv", out_bins), ("unsupported_coverage.tsv", out_cov)):
        with open(os.path.join(a.outdir, fn), "w") as f:
            f.write("".join("\t".join(fmt(x) for x in r) + "\n" for r in rows))
    print(f"age_bias_checks: {a.outdir}/age_bias_checks.tsv, unsupported_coverage.tsv")


if __name__ == "__main__":
    main()
