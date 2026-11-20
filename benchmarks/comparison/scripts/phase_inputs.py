#!/usr/bin/env python3
"""Phasing inputs for ltrharvest_phasing10.py from any tool's calls (the uniform path, same for every tool).

  calls (calls_to_bed.py --keep-ids) -> elements.fa (genome orientation, headers c0001:start-end, 1-based)
  -> Kmer2LTR (LTR pair + K2P per copy, from the sequence alone) -> annot.tsv
  -> elements.k.fa = copies Kmer2LTR paired ('pass' with a finite K2P); share kept -> stats.tsv
Optional --subsample N --seed S writes a random subset of the paired copies (power curve).
Run with the kmer2ltr env on PATH (Kmer2LTR).
"""
import argparse
import csv
import math
import os
import random
import subprocess
import sys


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
    a = off + (s // lb) * lw + s % lb
    b = off + ((e - 1) // lb) * lw + (e - 1) % lb + 1
    fh.seek(a)
    return fh.read(b - a).decode().replace("\n", "")


def read_fa(path):
    name, buf = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Uniform phasing inputs: calls -> FASTA -> Kmer2LTR K2P table.")
    ap.add_argument("--calls", help="calls TSV with normalized IDs (omit with --subsample-from)")
    ap.add_argument("--genome", help="normalized genome FASTA (+ .fai)")
    ap.add_argument("--outdir", required=True, help="output directory")
    ap.add_argument("-t", "--threads", type=int, default=16, help="Kmer2LTR worker processes (default 16)")
    ap.add_argument("--subsample-from", help="existing outdir: subsample its elements.k.fa instead of rebuilding")
    ap.add_argument("--subsample", type=int, help="number of paired copies to keep (random)")
    ap.add_argument("--seed", type=int, default=1, help="random seed for --subsample (default 1)")
    ap.add_argument("--max-len", type=int, default=0,
                    help="spans longer than this (bp) skip Kmer2LTR and get status too_long, K2P NA (default 0: off)")
    ap.add_argument("-v", "--verbose", action="store_true", help="(reserved)")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    if a.subsample_from:
        recs = list(read_fa(os.path.join(a.subsample_from, "elements.k.fa")))
        if a.subsample is None or a.subsample >= len(recs):
            sys.exit(f"error: --subsample must be < {len(recs)} paired copies")
        rnd = random.Random(a.seed)
        keep = rnd.sample(recs, a.subsample)
        with open(os.path.join(a.outdir, "elements.k.fa"), "w") as f:
            for n, s in keep:
                f.write(f">{n}\n{s}\n")
        for ext in ("tsv", "gff3"):                      # uniform path: Kmer2LTR TSV; LTR_peel native: GFF3
            src = os.path.join(a.subsample_from, f"annot.{ext}")
            dst = os.path.join(a.outdir, f"annot.{ext}")
            if os.path.exists(src) and not os.path.lexists(dst):
                os.symlink(os.path.abspath(src), dst)
        print(f"phase_inputs: subsample {a.subsample} of {len(recs)} (seed {a.seed})", file=sys.stderr)
        return
    if not (a.calls and a.genome):
        sys.exit("error: --calls and --genome are required unless --subsample-from")
    idx = fai(a.genome)
    rows = list(csv.DictReader(open(a.calls), delimiter="\t"))
    el = os.path.join(a.outdir, "elements.fa")
    seen = set()
    too_long = []
    with open(a.genome, "rb") as fh, open(el, "w") as f:
        for r in rows:
            s, e = int(r["start"]), int(r["end"])
            key = f"{r['chrom']}:{s + 1}-{e}"
            if key in seen:
                continue
            seen.add(key)
            if a.max_len and e - s > a.max_len:
                too_long.append((key, e - s))
                continue
            f.write(f">{key}\n{fetch(fh, idx, r['chrom'], s, e)}\n")
    tsv = os.path.join(a.outdir, "annot.tsv")
    r = subprocess.run(["Kmer2LTR", el, "-o", tsv, "-t", str(a.threads)], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"error: Kmer2LTR failed:\n{r.stderr[-2000:]}")
    if too_long:                                   # same columns as Kmer2LTR's rows; status too_long, everything NA
        with open(tsv) as f:
            cols = f.readline().rstrip("\n").split("\t")
        with open(tsv, "a") as f:
            for key, n in too_long:
                v = {c: "NA" for c in cols}
                v.update(seq_id=key, seq_len=str(n), status="too_long")
                f.write("\t".join(v[c] for c in cols) + "\n")
        print(f"phase_inputs: {len(too_long)} spans > {a.max_len} bp not run through Kmer2LTR (status too_long)",
              file=sys.stderr)
    ok = set()
    with open(tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                k = float(row["k2p"])
            except (TypeError, ValueError):
                continue
            if row.get("status") == "pass" and math.isfinite(k):
                ok.add(row["seq_id"])
    n = 0
    with open(os.path.join(a.outdir, "elements.k.fa"), "w") as f:
        for name, seq in read_fa(el):
            if name in ok:
                f.write(f">{name}\n{seq}\n")
                n += 1
    with open(os.path.join(a.outdir, "stats.tsv"), "w") as f:
        f.write(f"n_calls\tn_unique\tn_paired\tfrac_paired\n{len(rows)}\t{len(seen)}\t{n}\t{n / max(len(seen), 1):.4f}\n")
    print(f"phase_inputs: {len(seen)} copies, {n} with an LTR pair and K2P ({n / max(len(seen), 1):.1%})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
