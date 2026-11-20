#!/usr/bin/env python3
"""Normalize a genome FASTA identically for every benchmarked tool.

Streams the input (plain or .gz) and writes, in OUTDIR:
  genome.fa        upper-case, non-ACGTN bases -> N, 60-column lines, short IDs (<prefix>0001, ...; no '_')
  genome.fa.fai    exact samtools-style index
  names.tsv        new_id  original_id  length  n_count  nonacgtn_converted
  split/<id>.fa    one FASTA per sequence (with --split; Look4LTRs takes a directory of .fa files)

Why: tools restrict IDs (LTRpred splits on '_', LTR_retriever limits length, LTR_FINDER_parallel parses '_sub'),
treat soft-masked bases differently, or misread IUPAC codes (LTR_peel uses IUPAC letters as mask characters).
"""
import argparse
import gzip
import os
import sys

WIDTH = 60
ACGTN = set(b"ACGTN")
TRANS = bytes(c if c in ACGTN else ord("N") for c in range(256))


def opener(path):
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")


def records(path):
    """Yield (header_first_token, iterator-free sequence chunks list) one record at a time."""
    name, chunks = None, []
    with opener(path) as f:
        for line in f:
            if line.startswith(b">"):
                if name is not None:
                    yield name, chunks
                tok = line[1:].split()
                name, chunks = (tok[0].decode() if tok else ""), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, chunks


def main(argv=None):
    ap = argparse.ArgumentParser(description="Normalize a genome FASTA (upper-case, ACGTN, short IDs, .fai, name map).")
    ap.add_argument("fasta", help="input FASTA (.fa/.fna, plain or .gz)")
    ap.add_argument("outdir", help="output directory")
    ap.add_argument("--prefix", default="c", help="ID prefix; IDs are <prefix>0001... (default c)")
    ap.add_argument("--split", action="store_true", help="also write split/<id>.fa, one file per sequence")
    ap.add_argument("--min-len", type=int, default=0, help="drop sequences shorter than this (default 0 = keep all)")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-sequence progress")
    a = ap.parse_args(argv)
    if "_" in a.prefix or not a.prefix.isalnum():
        sys.exit("error: --prefix must be alphanumeric without '_'")
    os.makedirs(a.outdir, exist_ok=True)
    if a.split:
        os.makedirs(os.path.join(a.outdir, "split"), exist_ok=True)
    fa_path = os.path.join(a.outdir, "genome.fa")
    print(f"prep_genome: {a.fasta} -> {fa_path}", file=sys.stderr)
    n_seq = n_bp = n_conv = n_drop = 0
    seen = set()
    with open(fa_path + ".tmp", "wb") as fa, open(os.path.join(a.outdir, "genome.fa.fai.tmp"), "w") as fai, \
            open(os.path.join(a.outdir, "names.tsv.tmp"), "w") as nm:
        nm.write("new_id\toriginal_id\tlength\tn_count\tnonacgtn_converted\n")
        offset = 0
        for orig, chunks in records(a.fasta):
            if orig in seen:
                sys.exit(f"error: duplicate sequence ID {orig!r}")
            seen.add(orig)
            raw = b"".join(chunks).upper()
            seq = raw.translate(TRANS)
            conv = sum(1 for x, y in zip(raw, seq) if x != y) if raw != seq else 0
            if len(seq) == 0:
                sys.exit(f"error: empty sequence {orig!r}")
            if len(seq) < a.min_len:
                n_drop += 1
                continue
            n_seq += 1
            new = f"{a.prefix}{n_seq:04d}"
            head = f">{new}\n".encode()
            body = b"".join(seq[i:i + WIDTH] + b"\n" for i in range(0, len(seq), WIDTH))
            fa.write(head)
            fa.write(body)
            fai.write(f"{new}\t{len(seq)}\t{offset + len(head)}\t{WIDTH}\t{WIDTH + 1}\n")
            offset += len(head) + len(body)
            nm.write(f"{new}\t{orig}\t{len(seq)}\t{seq.count(b'N')}\t{conv}\n")
            if a.split:
                with open(os.path.join(a.outdir, "split", f"{new}.fa"), "wb") as sf:
                    sf.write(head)
                    sf.write(body)
            n_bp += len(seq)
            n_conv += conv
            if a.verbose:
                print(f"  {orig} -> {new} {len(seq):,} bp", file=sys.stderr)
    for f in ("genome.fa", "genome.fa.fai", "names.tsv"):
        os.replace(os.path.join(a.outdir, f + ".tmp"), os.path.join(a.outdir, f))
    print(f"prep_genome: {n_seq} sequences, {n_bp:,} bp, {n_conv:,} non-ACGTN -> N, {n_drop} dropped (<{a.min_len} bp)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
