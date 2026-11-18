#!/usr/bin/env python3
"""Coordinate sanity check for a calls table: share of calls whose span starts with TG and ends with CA at the
reported position and at +-1/+-2 bp shifts (TG...CA reads the same on both strands, so this is strand-free).
A correct parser peaks at shift 0. Uses the normalized genome (calls with --keep-ids) or any FASTA + .fai."""
import argparse
import csv
import sys


def fai(path):
    out = {}
    with open(path + ".fai") as f:
        for line in f:
            c = line.split("\t")
            out[c[0]] = (int(c[1]), int(c[2]), int(c[3]), int(c[4]))
    return out


def fetch(fh, idx, chrom, s, e):
    """0-based half-open slice via .fai."""
    L, off, lb, lw = idx[chrom]
    s, e = max(0, s), min(L, e)
    if e <= s:
        return ""
    a = off + (s // lb) * lw + s % lb
    b = off + ((e - 1) // lb) * lw + (e - 1) % lb + 1
    fh.seek(a)
    return fh.read(b - a).decode().replace("\n", "").upper()


def main(argv=None):
    ap = argparse.ArgumentParser(description="TG...CA rate at reported call ends and at small shifts.")
    ap.add_argument("calls", help="calls TSV (calls_to_bed.py --keep-ids)")
    ap.add_argument("genome", help="genome FASTA with .fai (normalized IDs)")
    ap.add_argument("--max", type=int, default=2000, help="calls to sample (default 2000)")
    a = ap.parse_args(argv)
    idx = fai(a.genome)
    rows = list(csv.DictReader(open(a.calls), delimiter="\t"))[: a.max]
    with open(a.genome, "rb") as fh:
        res = {}
        for sh in (-2, -1, 0, 1, 2):
            n = ok_s = ok_e = 0
            for r in rows:
                s, e = int(r["start"]) + sh, int(r["end"]) + sh
                n += 1
                ok_s += fetch(fh, idx, r["chrom"], s, s + 2) == "TG"
                ok_e += fetch(fh, idx, r["chrom"], e - 2, e) == "CA"
            res[sh] = (ok_s / max(n, 1), ok_e / max(n, 1))
    print("shift\tTG_at_start\tCA_at_end\tn=" + str(len(rows)))
    for sh, (x, y) in res.items():
        print(f"{sh:+d}\t{x:.3f}\t{y:.3f}")


if __name__ == "__main__":
    main()
