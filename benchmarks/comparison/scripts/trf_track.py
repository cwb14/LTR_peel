#!/usr/bin/env python3
"""Genome-wide tandem-repeat track (BED, merged, 0-based) with TRF at EDTA cleanup_tandem.pl's defaults
(trf 2 7 7 80 10 1000 2000 -ngs -h -l 6). The genome is cut into 5 Mb windows overlapping by 20 kb, run in parallel;
intervals are shifted back to sequence coordinates and merged. One track per genome serves every tool and the nulls.
usage: trf_track.py --genome G.fa --out trf.bed [-t 16] [-v]
"""
import argparse
import multiprocessing as mp
import os
import subprocess
import sys
import tempfile

WIN, OVL = 5_000_000, 20_000
TRF_ARGS = ["2", "7", "7", "80", "10", "1000", "2000", "-ngs", "-h", "-l", "6"]


def fai(path):
    out = {}
    with open(path + ".fai") as f:
        for line in f:
            c = line.split("\t")
            out[c[0]] = (int(c[1]), int(c[2]), int(c[3]), int(c[4]))
    return out


def fetch(path, idx, chrom, s, e):
    L, off, lb, lw = idx[chrom]
    s, e = max(0, s), min(L, e)
    a = off + (s // lb) * lw + s % lb
    b = off + ((e - 1) // lb) * lw + (e - 1) % lb + 1
    with open(path, "rb") as fh:
        fh.seek(a)
        return fh.read(b - a).decode().replace("\n", "")


def parse_ngs(text, offset):
    """TRF -ngs output -> [(start0, end)] shifted by offset"""
    out = []
    for line in text.splitlines():
        if not line or line.startswith("@"):
            continue
        f = line.split()
        if len(f) > 2 and f[0].isdigit() and f[1].isdigit():
            out.append((int(f[0]) - 1 + offset, int(f[1]) + offset))
    return out


def run_window(job):
    trf, genome, idx, chrom, s, e = job
    seq = fetch(genome, idx, chrom, s, e)
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as d:
        fa = os.path.join(d, "w.fa")
        with open(fa, "w") as f:
            f.write(f">w\n{seq}\n")
        r = subprocess.run([trf, fa] + TRF_ARGS, cwd=d, capture_output=True, text=True)
    return chrom, parse_ngs(r.stdout, s)


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Genome-wide TRF tandem-repeat track (merged BED).")
    ap.add_argument("--genome", required=True, help="genome FASTA (+ .fai)")
    ap.add_argument("--out", required=True, help="output BED")
    ap.add_argument("-t", "--threads", type=int, default=16, help="parallel TRF processes")
    ap.add_argument("--trf", default="trf", help="trf executable")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-sequence progress")
    a = ap.parse_args(argv)
    idx = fai(a.genome)
    jobs = [(a.trf, a.genome, idx, c, s, min(L, s + WIN + OVL)) for c, (L, *_) in idx.items()
            for s in range(0, L, WIN)]
    print(f"trf_track: {len(jobs)} windows, {a.threads} processes", file=sys.stderr)
    iv = {c: [] for c in idx}
    with mp.Pool(a.threads) as pool:
        for k, (c, hits) in enumerate(pool.imap_unordered(run_window, jobs), 1):
            iv[c] += hits
            if a.verbose and k % 50 == 0:
                print(f"trf_track: {k}/{len(jobs)}", file=sys.stderr)
    n = bp = 0
    with open(a.out + ".part", "w") as f:
        for c in idx:
            for s, e in merge(iv[c]):
                f.write(f"{c}\t{s}\t{e}\n")
                n += 1
                bp += e - s
    os.replace(a.out + ".part", a.out)
    print(f"trf_track: {a.out} {n} intervals, {bp / 1e6:.2f} Mb", file=sys.stderr)


if __name__ == "__main__":
    main()
