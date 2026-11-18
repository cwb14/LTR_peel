#!/usr/bin/env python3
"""How far does cd-hit-est shrink a large candidate set? Runs build_library.py's steps 1-2 (extract, cleanup_tandem)
on one tool's calls, then cd-hit-est at the given identities, and reports counts and seconds (no cleanup_nested).
Run with the ltrr_bench env active.  ->  OUTDIR/sizes.tsv
"""
import argparse
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.join(os.environ["BENCH_ROOT"], "scripts"))
import build_library as bl  # noqa: E402


def count(fa):
    return sum(1 for line in open(fa) if line.startswith(">"))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Candidate counts after tandem cleanup and cd-hit-est (no cleanup_nested).")
    ap.add_argument("--calls", required=True, help="calls TSV with normalized IDs")
    ap.add_argument("--genome", required=True, help="normalized genome FASTA (+ .fai)")
    ap.add_argument("--outdir", required=True, help="output directory")
    ap.add_argument("--prefix", required=True, help="entry name prefix")
    ap.add_argument("--ids", type=float, nargs="+", default=[0.95], help="cd-hit-est identities (default 0.95)")
    ap.add_argument("-t", "--threads", type=int, default=16, help="cd-hit-est threads (default 16)")
    a = ap.parse_args(argv)
    a.outdir = os.path.abspath(a.outdir)
    os.makedirs(a.outdir, exist_ok=True)
    rows = []
    t0 = time.time()
    raw = os.path.join(a.outdir, "raw.fa")
    import csv
    calls = list(csv.DictReader(open(a.calls), delimiter="\t"))
    idx = bl.fai(a.genome)
    with open(a.genome, "rb") as fh:
        bl.write_fa(raw, [(f"{a.prefix}_{i + 1}#LTR/{bl.superfamily(r.get('class'))}",
                           bl.fetch(fh, idx, r["chrom"], int(r["start"]), int(r["end"]))) for i, r in enumerate(calls)])
    rows.append(("raw", count(raw), round(time.time() - t0)))
    t0 = time.time()
    cln = raw + ".cln"
    with open(cln, "w") as out, open(os.path.join(a.outdir, "tandem.log"), "w") as log:
        subprocess.run(["perl", f"{bl.EDTA_BIN}/cleanup_tandem.pl", "-misschar", "N", "-nc", "50000", "-nr", "0.9",
                        "-minlen", "100", "-minscore", "3000", "-trf", "1", "-cleanN", "1", "-f", raw],
                       stdout=out, stderr=log, cwd=a.outdir, check=True)
    rows.append(("tandem", count(cln), round(time.time() - t0)))
    for c in a.ids:
        t0 = time.time()
        o = os.path.join(a.outdir, f"cd{int(round(c * 100))}.fa")
        n = 5 if c < 0.85 else (8 if c < 0.9 else 10)    # cd-hit-est word size for the identity range
        with open(o + ".log", "w") as log:
            subprocess.run(["cd-hit-est", "-i", cln, "-o", o, "-c", str(c), "-n", str(n), "-T", str(a.threads),
                            "-M", "0", "-d", "0"], stdout=log, stderr=log, check=True)
        rows.append((f"cd{int(round(c * 100))}", count(o), round(time.time() - t0)))
    with open(os.path.join(a.outdir, "sizes.tsv"), "w") as f:
        f.write("step\tn_seq\tseconds\n" + "".join(f"{s}\t{n}\t{t}\n" for s, n, t in rows))
    print("\n".join(f"{s}\t{n}\t{t}" for s, n, t in rows), file=sys.stderr)


if __name__ == "__main__":
    main()
