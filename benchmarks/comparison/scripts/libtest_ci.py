#!/usr/bin/env python3
"""Block-bootstrap 95% CIs for EDTA lib-test's LTR row (masking F1 vs the curated annotation).

Recomputes lib-test's LTR-category base counts exactly (reference and test RepeatMasker .out hits whose class matches
lib-test's 'ltr' regex; TP = both, FP = test only, FN = reference only), in genome blocks of --block bp. Blocks are
resampled with replacement (the same draws for every library, so differences are paired) and F1 = 2TP/(2TP+FP+FN) is
recomputed per draw. Blocks rather than whole sequences: A. thaliana has 5 nuclear chromosomes, too few to resample;
LTR-RTs cluster at Mb scale (pericentromeres), so blocks are large (default 5 Mb, as in genome-wide block jackknifes).

Every library under mask/<ds>/ with an eval/lib/<ds>/<lib>.tsv is included (superseded *.uniform skipped); totals are
checked against lib-test's own TP/FP/FN (column 'check').
usage: libtest_ci.py DS [--block 5000000] [--reps 2000] [--ref-lib ltr_peel.intact]
"""
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np

B = os.environ["BENCH_ROOT"]
LTR_RE = re.compile("|".join(["RLG", "RLC", "RLB", "RLR", "RLE", "LTR", "RLX", "Gypsy", "Copia"]), re.I)  # lib-test


def ltr_hits(path):
    """(chrom, start0, end) of LTR-category hits in a RepeatMasker .out, as lib-test selects them."""
    out = {}
    with open(path) as f:
        for line in f:
            t = line.split()
            if len(t) < 11 or not t[0].isdigit():
                continue
            if t[10] != "Unspecified" and LTR_RE.search(t[10]):
                out.setdefault(t[4], []).append((int(t[5]) - 1, int(t[6])))
    return out


def cover(n, ivs):
    a = np.zeros(n, bool)
    for s, e in ivs or ():
        a[s:e] = True
    return a


def block_sums(x, block):
    """sum of a boolean array per block (last block may be shorter); no copy of x."""
    return np.array([np.count_nonzero(x[i:i + block]) for i in range(0, len(x), block)], np.int64)


def libtest_totals(path):
    """TP, FP, FN from the lib-test row written by the workflow (first data row)."""
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["TP"] != "TP":
                return int(r["TP"]), int(r["FP"]), int(r["FN"])
    return None


def f1(tp, fp, fn):
    d = 2 * tp + fp + fn
    return np.where(d > 0, 2 * tp / np.maximum(d, 1), np.nan)


def main(argv=None):  # noqa: C901
    ap = argparse.ArgumentParser(description="Block-bootstrap CIs for lib-test's LTR F1, paired across libraries.")
    ap.add_argument("ds", help="dataset with mask/<ds>/curated/<ds>.fa.out and eval/lib/<ds>/*.tsv")
    ap.add_argument("--block", type=int, default=5_000_000, help="block length in bp (default 5 Mb)")
    ap.add_argument("--reps", type=int, default=2000, help="bootstrap draws (default 2000)")
    ap.add_argument("--seed", type=int, default=1, help="random seed (default 1)")
    ap.add_argument("--ref-lib", default="ltr_peel.intact", help="library the paired differences are taken against")
    ap.add_argument("-o", "--out", help="output TSV (default eval/lib/<ds>.ci.tsv)")
    ap.add_argument("--base", default=None, help="benchmark root (default: this benchmark)")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-library progress")
    a = ap.parse_args(argv)
    base = a.base or B
    out = a.out or f"{base}/eval/lib/{a.ds}.ci.tsv"
    libs = sorted(os.path.basename(p)[:-4] for p in glob.glob(f"{base}/eval/lib/{a.ds}/*.tsv")
                  if not p.endswith(".uniform.tsv") and os.path.getsize(p) > 0
                  and os.path.exists(f"{base}/mask/{a.ds}/{os.path.basename(p)[:-4]}/{a.ds}.fa.out"))
    if not libs:
        sys.exit(f"error: no scored libraries for {a.ds}")
    lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(f"{base}/data/{a.ds}/genome.fa.fai")}
    ref = ltr_hits(f"{base}/mask/{a.ds}/curated/{a.ds}.fa.out")
    refcov = {c: cover(n, ref.get(c)) for c, n in lens.items()}
    counts = []                                          # per library: (n_blocks, 3) TP FP FN
    big = max(lens.values())
    tbuf, xbuf = np.zeros(big, bool), np.zeros(big, bool)  # reused per sequence (avoids re-zeroing fresh pages)
    for lib in libs:
        tst = ltr_hits(f"{base}/mask/{a.ds}/{lib}/{a.ds}.fa.out")
        cols = []
        for c, n in lens.items():
            r, t, x = refcov[c], tbuf[:n], xbuf[:n]
            t[:] = False
            for s, e in tst.get(c, ()):
                t[s:e] = True
            k = [block_sums(np.logical_and(r, t, out=x), a.block),   # TP
                 block_sums(np.greater(t, r, out=x), a.block),       # FP: test and not reference
                 block_sums(np.greater(r, t, out=x), a.block)]       # FN: reference and not test
            cols.append(np.stack(k, 1))
        counts.append(np.concatenate(cols))
        if a.verbose:
            print(f"{a.ds} {lib}: TP/FP/FN {counts[-1].sum(0).tolist()}", file=sys.stderr)
    C = np.stack(counts).astype(np.int64)                 # libs x blocks x 3
    nb = C.shape[1]
    rng = np.random.default_rng(a.seed)
    M = np.stack([np.bincount(rng.integers(0, nb, nb), minlength=nb) for _ in range(a.reps)])  # reps x blocks
    S = np.einsum("rb,lbk->lrk", M, C)                    # libs x reps x 3
    F = f1(S[..., 0], S[..., 1], S[..., 2])               # libs x reps
    tot = C.sum(1)
    point = f1(tot[:, 0], tot[:, 1], tot[:, 2])
    ri = libs.index(a.ref_lib) if a.ref_lib in libs else None
    bad = 0
    with open(out, "w") as fo:
        fo.write("ds\tlib\tTP\tFP\tFN\tprec\tsens\tF1\tF1_lo\tF1_hi\tref_lib\tdF1\tdF1_lo\tdF1_hi\tn_blocks\tblock_bp"
                 "\tcheck\n")
        for i, lib in enumerate(libs):
            tp, fp, fn = (int(x) for x in tot[i])
            lt = libtest_totals(f"{base}/eval/lib/{a.ds}/{lib}.tsv")
            chk = "na" if lt is None else ("ok" if lt == (tp, fp, fn) else f"MISMATCH libtest={lt}")
            bad += chk.startswith("MISMATCH")
            lo, hi = np.nanpercentile(F[i], [2.5, 97.5])
            if ri is not None:
                d = F[i] - F[ri]
                dl, dh = np.nanpercentile(d, [2.5, 97.5])
                dd = f"{a.ref_lib}\t{point[i] - point[ri]:.4f}\t{dl:.4f}\t{dh:.4f}"
            else:
                dd = ".\t.\t.\t."
            fo.write(f"{a.ds}\t{lib}\t{tp}\t{fp}\t{fn}\t{tp / max(tp + fp, 1):.4f}\t{tp / max(tp + fn, 1):.4f}\t"
                     f"{point[i]:.4f}\t{lo:.4f}\t{hi:.4f}\t{dd}\t{nb}\t{a.block}\t{chk}\n")
    print(f"libtest_ci: {a.ds}: {len(libs)} libraries, {nb} blocks of {a.block} bp -> {out}"
          + (f"; {bad} totals differ from lib-test" if bad else "; all totals match lib-test"), file=sys.stderr)


if __name__ == "__main__":
    main()
