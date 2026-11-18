#!/usr/bin/env python3
"""Approach 4 (post hoc): per-call library support and library-independent evidence, with a local null.

For each call of one tool on one real genome (rules frozen before any per-tool number):
  supported    >= 80% of the span covered by curated LTR hits (the proxies.py rule)
  k2p, age     Kmer2LTR K2P (status pass, finite) and its bin: <0.01, 0.01-0.03, 0.03-0.1, >=0.1, NA
  coding       >= 2 distinct REXdb Class_I|LTR domain names of one (strand, superfamily) fully inside the internal
               region (Kmer2LTR ltr5_end..ltr3_start; else the middle 80% of the span)
  gene         >= 300 bp non-TE protein-coding CDS in the span      (NA if the genome has no gene annotation)
  nonltr       >= 50% of the LTR bases (span outside the internal region) curated DNA/RC/LINE/SINE
  sat          >= 50% of the span curated satellite/simple/centromeric
  trf          >= 50% of the span in the genome-wide TRF track
  junk         any of gene, nonltr, sat, trf
  tgca, tsd    TG..CA ends; exact 4-6 bp TSD
Null: the span (and its internal region) moved by d = +-U(20 kb, 200 kb) on the same sequence (default_rng(seed), calls
in file order); redrawn (<= 100 tries) when off-sequence, overlapping the call or > 10% N; the same flags on it.
usage: adjudicate.py --ds DS --tool TOOL --out eval/adj/DS/TOOL.calls.tsv [--seed 1] [-v]
"""
import argparse
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence as ev  # noqa: E402
import proxies  # noqa: E402

B = os.environ["BENCH_ROOT"]
BINS = [(0.0, 0.01, "<0.01"), (0.01, 0.03, "0.01-0.03"), (0.03, 0.1, "0.03-0.1"), (0.1, math.inf, ">=0.1")]


def age_bin(k2p):
    if k2p is None:
        return "NA"
    for lo, hi, name in BINS:
        if lo <= k2p < hi:
            return name
    return "NA"


def read_k2p(path):
    """{(chrom, start0, end): (k2p or None, internal (a, b) offsets from start or None)} from Kmer2LTR annot.tsv"""
    out = {}
    if not path or not os.path.exists(path):
        return out
    for r in csv.DictReader(open(path), delimiter="\t"):
        c, rng = r["seq_id"].rsplit(":", 1)
        s1, e = rng.split("-")
        key = (c, int(s1) - 1, int(e))
        k2p = None
        if r["status"] == "pass":
            try:
                v = float(r["k2p"])
                k2p = v if math.isfinite(v) else None
            except ValueError:
                pass
        inner = None
        if r["status"] in ("pass", "weak_pair"):
            try:
                a, b = int(r["ltr5_end"]), int(r["ltr3_start"]) - 1
                if 0 < a < b:
                    inner = (a, b)
            except ValueError:
                pass
        out[key] = (k2p, inner)
    return out


def bed_mask(path, lens):
    arr = {c: np.zeros(n, bool) for c, n in lens.items()}
    if path and os.path.exists(path):
        for line in open(path):
            c, s, e = line.split("\t")[:3]
            if c in arr:
                arr[c][int(s):int(e)] = True
    return arr


class Genome:
    def __init__(self, ds, root=B):
        self.path = f"{root}/data/{ds}/genome.fa"
        self.idx = proxies.fai(self.path)
        self.lens = {c: v[0] for c, v in self.idx.items()}
        ref = f"{root}/mask/{ds}/curated/{ds}.fa.out"
        self.ltr_iv = proxies.ref_ltr(ref)
        self.ltr_starts = {k: [s for s, _ in v] for k, v in self.ltr_iv.items()}
        self.cls = ev.curated_classes(ref, self.lens)
        gff = ev.gene_gff(f"{root}/data/{ds}/genes")
        self.cds = ev.cds_mask(gff, f"{root}/data/{ds}/names.tsv", self.lens, ev.TE_NOTE_V2) if gff else None
        self.trf = bed_mask(f"{root}/eval/trf/{ds}/trf.bed", self.lens)
        self.doms = ev.dante_domains(f"{root}/runs/{ds}/dante_ltr/dante.gff3")
        self.fh = open(self.path, "rb")

    def seq(self, c, s, e):
        return proxies.fetch(self.fh, self.idx, c, s, e)

    def flags(self, c, s, e, a, b):
        """evidence for span [s, e) with internal region [a, b)"""
        n = e - s
        cl = self.cls[c]
        ltr_bases = np.concatenate([cl[s:a], cl[b:e]])
        f_nonltr = float(np.isin(ltr_bases, (ev.DNA, ev.RC, ev.LINE, ev.SINE)).mean()) if ltr_bases.size else 0.0
        f_sat = float((cl[s:e] == ev.SAT).mean())
        f_trf = float(self.trf[c][s:e].mean())
        gene_bp = int(self.cds[c][s:e].sum()) if self.cds is not None else None
        sup = proxies.covered(self.ltr_iv, self.ltr_starts, c, s, e) >= 0.8 * n
        r = dict(supported=int(sup), coding=int(ev.coding_unit(self.doms, c, a, b)),
                 gene_bp="NA" if gene_bp is None else gene_bp,
                 gene="NA" if gene_bp is None else int(gene_bp >= 300),
                 nonltr_frac=round(f_nonltr, 3), nonltr=int(f_nonltr >= 0.5), sat_frac=round(f_sat, 3),
                 sat=int(f_sat >= 0.5), trf_frac=round(f_trf, 3), trf=int(f_trf >= 0.5))
        r["junk"] = int(any(r[k] == 1 for k in ("gene", "nonltr", "sat", "trf")))
        return r

    def hallmarks(self, c, s, e):
        tg = int(self.seq(c, s, s + 2) == "TG" and self.seq(c, e - 2, e) == "CA")
        tsd = 0
        for k in (4, 5, 6):
            left, right = self.seq(c, s - k, s), self.seq(c, e, e + k)
            if len(left) == k and left == right and "N" not in left:
                tsd = 1
                break
        return tg, tsd


def draw_null(g, rng, c, s, e, tries=100):
    L = g.lens[c]
    for _ in range(tries):
        d = int(rng.uniform(20_000, 200_000)) * (1 if rng.random() < 0.5 else -1)
        ns, ne = s + d, e + d
        if ns < 0 or ne > L or (ns < e and s < ne):
            continue
        q = g.seq(c, ns, ne)
        if q.count("N") > 0.1 * (ne - ns):
            continue
        return d
    return None


FLAGS = ["supported", "coding", "gene_bp", "gene", "nonltr_frac", "nonltr", "sat_frac", "sat", "trf_frac", "trf", "junk"]


def main(argv=None):
    ap = argparse.ArgumentParser(description="Per-call library support, independent evidence and local null.")
    ap.add_argument("--ds", required=True, help="real dataset")
    ap.add_argument("--tool", required=True, help="tool (calls/<ds>/<tool>.ids.tsv)")
    ap.add_argument("--calls", help="calls TSV (default calls/<ds>/<tool>.ids.tsv)")
    ap.add_argument("--k2p", help="Kmer2LTR annot.tsv (default eval/k2p/<ds>/<tool>/annot.tsv)")
    ap.add_argument("--out", required=True, help="output TSV, one row per call")
    ap.add_argument("--seed", type=int, default=1, help="null seed (default 1)")
    ap.add_argument("--root", default=B, help="benchmark root (default: the bench directory)")
    ap.add_argument("-v", "--verbose", action="store_true", help="progress")
    a = ap.parse_args(argv)
    calls = a.calls or f"{a.root}/calls/{a.ds}/{a.tool}.ids.tsv"
    k2p = read_k2p(a.k2p or f"{a.root}/eval/k2p/{a.ds}/{a.tool}/annot.tsv")
    g = Genome(a.ds, a.root)
    if a.verbose:
        print(f"adjudicate: {a.ds}/{a.tool} tracks loaded (genes: {g.cds is not None})", file=sys.stderr)
    rng = np.random.default_rng(a.seed)
    rows = list(csv.DictReader(open(calls), delimiter="\t"))
    cols = ["id", "chrom", "start", "end", "len", "k2p", "age", "inner", "int_start", "int_end", "tgca", "tsd"] + \
        FLAGS + ["null_shift"] + [f"null_{k}" for k in FLAGS]
    n_drop = 0
    with open(a.out + ".part", "w") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for i, r in enumerate(rows):
            c, s, e = r["chrom"], int(r["start"]), int(r["end"])
            kv, inner = k2p.get((c, s, e), (None, None))
            if inner and inner[1] < e - s:
                ia, ib, src = s + inner[0], s + inner[1], "kmer2ltr"
            else:
                ia, ib, src = s + (e - s) // 10, e - (e - s) // 10, "middle80"
            tg, tsd = g.hallmarks(c, s, e)
            out = dict(id=r["id"], chrom=c, start=s, end=e, len=e - s, k2p="NA" if kv is None else f"{kv:.5f}",
                       age=age_bin(kv), inner=src, int_start=ia, int_end=ib, tgca=tg, tsd=tsd, **g.flags(c, s, e, ia, ib))
            d = draw_null(g, rng, c, s, e)
            if d is None:
                n_drop += 1
                out.update(null_shift="NA", **{f"null_{k}": "NA" for k in FLAGS})
            else:
                nf = g.flags(c, s + d, e + d, ia + d, ib + d)
                out.update(null_shift=d, **{f"null_{k}": v for k, v in nf.items()})
            w.writerow(out)
            if a.verbose and (i + 1) % 10000 == 0:
                print(f"adjudicate: {i + 1}/{len(rows)}", file=sys.stderr)
    os.replace(a.out + ".part", a.out)
    print(f"adjudicate: {a.ds}/{a.tool} {len(rows)} calls -> {a.out} (nulls dropped: {n_drop})", file=sys.stderr)


if __name__ == "__main__":
    main()
