#!/usr/bin/env python3
"""Element-level quality proxies for calls on a real genome (no truth), computed identically for every tool.

  n_calls, median_len
  frac_tgca      span starts with TG and ends with CA (strand-free: TG...CA is its own reverse complement)
  frac_tsd       exact target-site duplication: the k bp left of the span equal the k bp right of it, k in 4..6
                 (chance ~ 4^-4 + 4^-5 + 4^-6 ~ 0.5%)
  frac_curated   >= 80% of the span masked as LTR by the curated reference annotation (RepeatMasker .out, lib-test LTR
                 pattern); a proxy for element-level precision
"""
import argparse
import bisect
import csv
import re
import statistics
import sys
from collections import defaultdict

LTR = re.compile("RLG|RLC|RLB|RLR|RLE|LTR|RLX|Gypsy|Copia", re.I)


def fai(path):
    out = {}
    with open(path + ".fai") as f:
        for line in f:
            c = line.split("\t")
            out[c[0]] = (int(c[1]), int(c[2]), int(c[3]), int(c[4]))
    return out


def fetch(fh, idx, chrom, s, e):
    L, off, lb, lw = idx[chrom]
    s, e = max(0, s), min(L, e)
    if e <= s:
        return ""
    a = off + (s // lb) * lw + s % lb
    b = off + ((e - 1) // lb) * lw + (e - 1) % lb + 1
    fh.seek(a)
    return fh.read(b - a).decode().replace("\n", "").upper()


def ref_ltr(path):
    iv = defaultdict(list)
    with open(path) as f:
        for line in f:
            c = line.split()
            if len(c) >= 11 and c[0].isdigit() and LTR.search(c[10]):
                iv[c[4]].append((int(c[5]) - 1, int(c[6])))
    merged = {}
    for k, v in iv.items():
        out = []
        for s, e in sorted(v):
            if out and s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[k] = out
    return merged


def covered(merged, starts, chrom, s, e):
    v = merged.get(chrom)
    if not v:
        return 0
    i = max(0, bisect.bisect_right(starts[chrom], s) - 1)
    tot = 0
    while i < len(v) and v[i][0] < e:
        tot += max(0, min(e, v[i][1]) - max(s, v[i][0]))
        i += 1
    return tot


def main(argv=None):
    ap = argparse.ArgumentParser(description="Sequence-based quality proxies for LTR-RT calls on a real genome.")
    ap.add_argument("--calls", required=True, help="calls TSV with normalized IDs (--keep-ids)")
    ap.add_argument("--genome", required=True, help="normalized genome FASTA (+ .fai)")
    ap.add_argument("--ref", required=True, help="curated-library .out on the same genome")
    ap.add_argument("--label", default="", help="label column")
    ap.add_argument("--out", required=True, help="output TSV row")
    a = ap.parse_args(argv)
    idx = fai(a.genome)
    merged = ref_ltr(a.ref)
    starts = {k: [s for s, _ in v] for k, v in merged.items()}
    rows = list(csv.DictReader(open(a.calls), delimiter="\t"))
    n = tg = tsd = cur = 0
    lens = []
    with open(a.genome, "rb") as fh:
        for r in rows:
            c, s, e = r["chrom"], int(r["start"]), int(r["end"])
            n += 1
            lens.append(e - s)
            tg += fetch(fh, idx, c, s, s + 2) == "TG" and fetch(fh, idx, c, e - 2, e) == "CA"
            for k in (4, 5, 6):
                left, right = fetch(fh, idx, c, s - k, s), fetch(fh, idx, c, e, e + k)
                if len(left) == k and left == right and "N" not in left:
                    tsd += 1
                    break
            cur += covered(merged, starts, c, s, e) >= 0.8 * (e - s)
    with open(a.out, "w") as f:
        f.write("label\tn_calls\tmedian_len\tfrac_tgca\tfrac_tsd\tfrac_curated\n")
        if n:
            f.write(f"{a.label}\t{n}\t{statistics.median(lens):.0f}\t{tg / n:.4f}\t{tsd / n:.4f}\t{cur / n:.4f}\n")
        else:
            f.write(f"{a.label}\t0\tnan\tnan\tnan\tnan\n")
    print(f"proxies: {a.label} n={n}", file=sys.stderr)


if __name__ == "__main__":
    main()
