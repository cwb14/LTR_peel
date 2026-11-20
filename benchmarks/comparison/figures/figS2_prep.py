#!/usr/bin/env python3
"""Data preparation for Fig. S2 (figS2_masking): the rice example window of panel A.
Cache written to this directory's tmp/.

Imported by figS2_masking.py (which calls init_bench(), then windows() when the cached table is missing); can also be
run alone (needs BENCH_ROOT):
  figS2_prep.py windows              rice 200-kb windows on assembled chromosomes, ranked for the example (panel A)
Coordinates: RepeatMasker .out 1-based closed -> 0-based half-open on read (start - 1, end); eval/adj is 0-based
half-open; MSU gff3 1-based closed -> 0-based half-open on read.
"""
import csv
import sys as _sys
_sys.dont_write_bytecode = True          # no __pycache__ next to the read-only imported modules
import gzip
import re
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TMP = HERE / "tmp"
B = libtest_ci = evidence = None        # set by init_bench()


def init_bench():
    """point B at $BENCH_ROOT and import its read-only scripts: libtest_ci (LTR_RE and ltr_hits reproduce lib-test's
    LTR class selection) and evidence (TE_NOTE_V2 for transposon-derived genes). Returns B."""
    global B, libtest_ci, evidence
    if B is None:
        if "BENCH_ROOT" not in os.environ:
            sys.exit("error: BENCH_ROOT is not set (the full dataset is needed unless --from-source-data)")
        B = Path(os.environ["BENCH_ROOT"])
        sys.path.insert(0, str(B / "scripts"))
        import libtest_ci as _lt
        import evidence as _ev
        libtest_ci, evidence = _lt, _ev
    TMP.mkdir(exist_ok=True)
    return B


DATASETS = ["arab", "rice", "finch", "maize", "gorilla"]
TOOLS = ["ltr_peel", "ltr_retriever", "ltrharvest", "ltr_finder", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]


def log(msg):
    print(f"[figS2_prep] {msg}", file=sys.stderr, flush=True)


def fai(ds):
    return [(l.split("\t")[0], int(l.split("\t")[1]), int(l.split("\t")[2]), int(l.split("\t")[3]),
             int(l.split("\t")[4])) for l in open(B / f"data/{ds}/genome.fa.fai")]


def merge(iv):
    """sorted, merged (start, end) array from an iterable of 0-based half-open intervals"""
    a = np.array(sorted(iv), dtype=np.int64).reshape(-1, 2)
    if len(a) == 0:
        return a
    out = [list(a[0])]
    for s, e in a[1:]:
        if s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return np.array(out, dtype=np.int64)


def intersect(a, b):
    """intersection of two merged interval arrays"""
    i = j = 0
    out = []
    while i < len(a) and j < len(b):
        s, e = max(a[i, 0], b[j, 0]), min(a[i, 1], b[j, 1])
        if s < e:
            out.append((s, e))
        if a[i, 1] < b[j, 1]:
            i += 1
        else:
            j += 1
    return np.array(out, dtype=np.int64).reshape(-1, 2)


def cum_cover(iv, x):
    """covered bp of merged intervals iv in [0, x) for each x (vectorized)"""
    x = np.asarray(x, dtype=np.int64)
    if len(iv) == 0:
        return np.zeros(len(x), np.int64)
    lens = iv[:, 1] - iv[:, 0]
    csum = np.concatenate([[0], np.cumsum(lens)])
    k = np.searchsorted(iv[:, 0], x, side="right")       # intervals starting before x
    full = csum[np.maximum(k - 1, 0)]                      # intervals fully before interval k-1
    part = np.where(k > 0, np.clip(x - iv[np.maximum(k - 1, 0), 0], 0, lens[np.maximum(k - 1, 0)]), 0)
    return np.where(k > 0, full + part, 0)


# ---------------------------------------------------------------- rice example window (panel A)
WIN = 200_000


def ren_map(ds):
    return {r["original_id"]: r["new_id"] for r in csv.DictReader(open(B / f"data/{ds}/names.tsv"), delimiter="\t")}


def genes(ds="rice"):
    """non-transposon genes (0-based half-open): chrom, start, end, strand, id; transposon-derived genes dropped by
    evidence.TE_NOTE_V2 on the gene's Note (the approach-4 rule)"""
    ren = ren_map(ds)
    out = []
    p = evidence.gene_gff(str(B / f"data/{ds}/genes"))
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9 or t[2] != "gene":
                continue
            if evidence.te_hit(evidence.TE_NOTE_V2, t[8]):
                continue
            c = ren.get(t[0])
            if c:
                m = re.search(r"ID=([^;]+)", t[8])
                out.append((c, int(t[3]) - 1, int(t[4]), t[6], m.group(1) if m else "."))
    return pd.DataFrame(out, columns=["chrom", "start", "end", "strand", "id"])


def windows(force=False):
    """rice 200-kb windows on c0001-c0012, with LTR_peel calls (by midpoint), their unsupported share, gene count
    and LTR_peel masking precision; distance to the genome-wide values (selection rule)"""
    path = TMP / "rice_windows.tsv"
    if path.exists() and not force:
        return pd.read_csv(path, sep="\t")
    ds = "rice"
    seqs = [s for s in fai(ds) if s[0] <= "c0012"]
    adj = pd.read_csv(B / f"eval/adj/{ds}/ltr_peel.calls.tsv", sep="\t", usecols=["chrom", "start", "end",
                                                                                     "supported"])
    g = genes(ds)
    ref = libtest_ci.ltr_hits(str(B / f"mask/{ds}/curated/{ds}.fa.out"))
    tst = libtest_ci.ltr_hits(str(B / f"mask/{ds}/ltr_peel.intact/{ds}.fa.out"))
    rows = []
    for c, n, *_ in seqs:
        nw = n // WIN
        edges = np.arange(nw + 1) * WIN
        a = adj[adj.chrom == c]
        mid = ((a.start + a.end) // 2).to_numpy()
        kc = np.minimum(mid // WIN, nw)
        ncall = np.bincount(kc, minlength=nw + 1)[:nw]
        nuns = np.bincount(kc, weights=(a.supported == 0).to_numpy(), minlength=nw + 1)[:nw]
        gm = ((g[g.chrom == c].start + g[g.chrom == c].end) // 2).to_numpy()
        ngene = np.bincount(np.minimum(gm // WIN, nw), minlength=nw + 1)[:nw]
        r, t = merge(ref.get(c, ())), merge(tst.get(c, ()))
        both = intersect(r, t)
        mt, mb = np.diff(cum_cover(t, edges)), np.diff(cum_cover(both, edges))
        for k in range(nw):
            rows.append((c, k * WIN, (k + 1) * WIN, int(ncall[k]), int(nuns[k]), int(ngene[k]), int(mt[k]),
                         int(mb[k])))
    w = pd.DataFrame(rows, columns=["chrom", "start", "end", "n_calls", "n_unsup", "n_genes", "lq_masked",
                                    "lq_TP"])
    w["unsup_share"] = w.n_unsup / w.n_calls.where(w.n_calls > 0)
    w["lq_prec"] = w.lq_TP / w.lq_masked.where(w.lq_masked > 0)
    # targets: medians over windows for counts; genome-wide (pooled) values for the two shares, because per-window
    # shares of ~3 calls are coarse (0, 1/3, ...); the panel shows masking, so LTR_peel's masking precision in the
    # window is a fourth criterion (added to the three above). Scale: interquartile range over windows.
    tgt = {"n_calls": w.n_calls.median(), "n_genes": w.n_genes.median(),
           "unsup_share": w.n_unsup.sum() / w.n_calls.sum(), "lq_prec": w.lq_TP.sum() / w.lq_masked.sum()}
    d = 0
    for k, v in tgt.items():
        sc = w[k].quantile(0.75) - w[k].quantile(0.25) or 1.0
        d = d + ((w[k] - v).abs() / sc).fillna(np.inf)
    w["dist"] = d
    w = w.sort_values(["dist", "chrom", "start"]).reset_index(drop=True)
    for k, v in tgt.items():
        w.attrs[k] = v
    with open(TMP / "rice_window_targets.tsv", "w") as f:
        f.write("metric\ttarget\n" + "".join(f"{k}\t{v}\n" for k, v in tgt.items()))
    w.to_csv(path, sep="\t", index=False)
    return w


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] != "windows":
        sys.exit(__doc__)
    init_bench()
    print(windows(force=True).head(15).to_string())
