#!/usr/bin/env python3
"""Masking input from a tool's intact LTR-RT calls, or a tool's native library normalized.

intact (--calls --intact ID; the benchmark's primary masking input):
  each call's full element sequence (genome orientation) -> <name>#LTR/<superfamily|unknown>, then cd-hit-est -c ID
  removes near-identical copies only (masking is unchanged by them; experiments/libreduce). Writes OUT.fa, OUT.log.
EDTA 3-step builder (--calls without --intact; superseded, kept for provenance): Ou et al. 2019, Genome Biol 20:275,
Methods 'non-redundant test libraries'.

uniform (calls -> library):
  1. extract each call's full element sequence (genome orientation) -> <name>#LTR/<superfamily|unknown>
  2. cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -cleanN 1
     (TRF '2 7 7 80 10 3000 2000 -ngs -h -l 6' inside; removes tandem-repeat and N-rich candidates)
  (--prereduce ID: cd-hit-est -c ID -M 0 -d 0 between steps 2 and 3, for datasets whose candidate sets make step 3
   infeasible; applied to every tool of such a dataset)
  3. cleanup_nested.pl -in X -threads T -minlen 80 -cov 0.95 (EDTA.pl's invocation; iterates automatically until
     saturation, the current form of the paper's five iterations) -> non-redundant library
native (tool library -> library): keep LTR entries only, headers '<name>#LTR/<superfamily|unknown>'.
Both write OUT.fa and OUT.log; uniform also leaves its intermediates in OUT.work/.
"""
import argparse
import csv
import os
import re
import shutil
import subprocess
import sys

# EDTA scripts (cleanup_tandem.pl, cleanup_nested.pl): EDTA_BIN, default <parent of BENCH_ROOT>/EDTA/bin
EDTA_BIN = os.environ.get("EDTA_BIN") or os.path.join(os.path.dirname(os.environ.get("BENCH_ROOT", "").rstrip("/")), "EDTA", "bin")


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


def superfamily(cls):
    c = (cls or "").lower()
    if "copia" in c or "ty1" in c:
        return "Copia"
    if "gypsy" in c or "ty3" in c:
        return "Gypsy"
    return "unknown"


def write_fa(path, recs):
    with open(path, "w") as f:
        for name, seq in recs:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


def read_fa(path):
    name, buf = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def uniform(a):
    a.out = os.path.abspath(a.out)
    work = a.out + ".work"
    os.makedirs(work, exist_ok=True)
    idx = fai(a.genome)
    rows = list(csv.DictReader(open(a.calls), delimiter="\t"))
    if not rows:
        sys.exit(f"error: no calls in {a.calls}")
    recs = []
    with open(a.genome, "rb") as fh:
        for i, r in enumerate(rows):
            seq = fetch(fh, idx, r["chrom"], int(r["start"]), int(r["end"]))
            recs.append((f"{a.prefix}_{i + 1}#LTR/{superfamily(r.get('class'))}", seq))
    raw = os.path.join(work, "raw.fa")
    write_fa(raw, recs)
    log = open(a.out + ".log", "w")
    log.write(f"{len(recs)} call sequences extracted\n")
    cln = os.path.join(work, "raw.fa.cln")
    with open(cln, "w") as out:
        r = subprocess.run(["perl", f"{EDTA_BIN}/cleanup_tandem.pl", "-misschar", "N", "-nc", "50000", "-nr", "0.9",
                            "-minlen", "100", "-minscore", "3000", "-trf", "1", "-cleanN", "1", "-f", raw],
                           stdout=out, stderr=log, cwd=work)
    if r.returncode != 0:
        sys.exit(f"error: cleanup_tandem.pl failed (see {a.out}.log)")
    n_cln = sum(1 for _ in read_fa(cln))
    log.write(f"{n_cln} after cleanup_tandem\n")
    nest_in = cln
    if a.prereduce:                      # fallback for the >2 Gb genomes: drop near-identical copies first
        nest_in = cln + ".cdhit"
        r = subprocess.run(["cd-hit-est", "-i", cln, "-o", nest_in, "-c", str(a.prereduce), "-T", str(a.threads),
                            "-M", "0", "-d", "0"], stdout=log, stderr=log, cwd=work)
        if r.returncode != 0 or not os.path.exists(nest_in):
            sys.exit(f"error: cd-hit-est failed (see {a.out}.log)")
        n_cd = sum(1 for _ in read_fa(nest_in))
        log.write(f"{n_cd} after cd-hit-est -c {a.prereduce}\n")
    r = subprocess.run(["perl", f"{EDTA_BIN}/cleanup_nested.pl", "-in", nest_in, "-threads", str(a.threads),
                        "-minlen", "80", "-cov", "0.95"], stdout=log, stderr=log, cwd=work)
    final = nest_in + ".cln"
    if r.returncode != 0 or not os.path.exists(final):
        sys.exit(f"error: cleanup_nested.pl failed (see {a.out}.log)")
    lib = [(n.split()[0], s) for n, s in read_fa(final)]
    write_fa(a.out + ".fa", lib)
    log.write(f"{len(lib)} after cleanup_nested\n")
    log.close()
    print(f"build_library: {a.prefix}: {len(recs)} calls -> {n_cln} after tandem cleanup -> {len(lib)} library entries",
          file=sys.stderr)


def cdhit_word(c):
    """cd-hit-est word size for identity c (cd-hit user guide)."""
    for lo, n in ((0.95, 10), (0.90, 8), (0.88, 7), (0.85, 6), (0.80, 5), (0.75, 4)):
        if c >= lo:
            return n
    sys.exit(f"error: cd-hit-est identity {c} < 0.75 is not supported")


CDHIT_MAX_SEQ = 1_000_000  # this cd-hit-est build (4.8.1, ltrr_bench): 1,000,000 bp ok, 1,000,001 warns, 2 Mb segfaults


def intact(a):
    """Primary masking input: the tool's intact call set as sequences. The only change is removing near-identical
    copies with cd-hit-est (the masking is the same with or without them; experiments/libreduce)."""
    a.out = os.path.abspath(a.out)
    work = a.out + ".work"
    os.makedirs(work, exist_ok=True)
    idx = fai(a.genome)
    rows = list(csv.DictReader(open(a.calls), delimiter="\t"))
    if not rows:
        sys.exit(f"error: no calls in {a.calls}")
    raw = os.path.join(work, "calls.fa")
    longest = 0
    with open(a.genome, "rb") as fh:
        recs = []
        for i, r in enumerate(rows):
            seq = fetch(fh, idx, r["chrom"], int(r["start"]), int(r["end"]))
            longest = max(longest, len(seq))
            recs.append((f"{a.prefix}_{i + 1}#LTR/{superfamily(r.get('class'))}", seq))
    # cd-hit-est handles sequences up to its compiled MAX_SEQ and crashes beyond it (finch HiRepEat has
    # 2.15 Mb calls); longer calls skip the redundancy reduction and go into the library unchanged
    longs = [x for x in recs if len(x[1]) > CDHIT_MAX_SEQ]
    write_fa(raw, [x for x in recs if len(x[1]) <= CDHIT_MAX_SEQ])
    red = os.path.join(work, "calls.cdhit.fa")
    with open(a.out + ".log", "w") as log:
        log.write(f"{len(recs)} call sequences extracted (longest {longest} bp); {len(longs)} longer than "
                  f"{CDHIT_MAX_SEQ} bp kept without cd-hit-est\n")
        log.flush()
        lib = []
        if len(recs) > len(longs):
            r = subprocess.run(["cd-hit-est", "-i", raw, "-o", red, "-c", str(a.intact), "-n",
                                str(cdhit_word(a.intact)), "-T", str(a.threads), "-M", "0", "-d", "0"],
                               stdout=log, stderr=log)
            if r.returncode != 0 or not os.path.exists(red):
                sys.exit(f"error: cd-hit-est failed (see {a.out}.log)")
            lib = list(read_fa(red))
        lib += longs
        if not lib:
            sys.exit(f"error: cd-hit-est wrote no sequences (see {a.out}.log)")
        write_fa(a.out + ".fa", lib)
        log.write(f"{len(lib)} after cd-hit-est -c {a.intact}\n")
    print(f"build_library: {a.prefix}: {len(recs)} intact calls -> {len(lib)} after cd-hit-est -c {a.intact}",
          file=sys.stderr)


def native(a):
    keep, drop = [], 0
    for name, seq in read_fa(a.native):
        tok = name.split()[0]
        nm, _, cls = tok.partition("#")
        if cls and not re.search(r"LTR|Gypsy|Copia|Ty1|Ty3|ERV", cls, re.I):
            drop += 1
            continue
        if re.search(r"nonLTR|LINE|SINE", cls, re.I):
            drop += 1
            continue
        sf = superfamily(cls)
        keep.append((f"{a.prefix}_{nm}#LTR/{sf}", seq.upper()))
    if not keep:
        sys.exit(f"error: no LTR entries in {a.native}")
    write_fa(a.out + ".fa", keep)
    with open(a.out + ".log", "w") as f:
        f.write(f"native {a.native}: {len(keep)} LTR entries kept, {drop} non-LTR dropped\n")
    print(f"build_library: {a.prefix}: native {len(keep)} LTR entries ({drop} non-LTR dropped)", file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description="EDTA-protocol test library from calls, or normalized native library.")
    ap.add_argument("--out", required=True, help="output prefix (writes OUT.fa, OUT.log)")
    ap.add_argument("--prefix", required=True, help="name prefix for library entries (e.g. tool name)")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--calls", help="calls TSV with normalized IDs (calls_to_bed.py --keep-ids) -> uniform builder")
    g.add_argument("--native", help="a tool's own library FASTA -> keep LTR entries, normalize headers")
    ap.add_argument("--genome", help="normalized genome FASTA (+ .fai), required with --calls")
    ap.add_argument("-t", "--threads", type=int, default=16, help="threads for cleanup_nested.pl (default 16)")
    ap.add_argument("--prereduce", type=float, default=None, metavar="ID",
                    help="cd-hit-est -c ID on the tandem-cleaned candidates before cleanup_nested.pl (removes "
                         "near-identical copies; fallback for very large candidate sets). Default: off")
    ap.add_argument("--intact", type=float, default=None, metavar="ID",
                    help="with --calls: the intact call set itself, reduced only by cd-hit-est -c ID (near-identical "
                         "copies); no tandem or nested cleanup. Without it, --calls runs the EDTA 3-step builder")
    ap.add_argument("-v", "--verbose", action="store_true", help="(reserved)")
    a = ap.parse_args(argv)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    need = ("cd-hit-est",) if a.intact else ("perl", "trf", "blastn", "makeblastdb") + \
        (("cd-hit-est",) if a.prereduce else ())
    for tool in need:
        if a.calls and not shutil.which(tool):
            sys.exit(f"error: {tool} not on PATH (activate ltrr_bench)")
    if a.calls:
        if not a.genome:
            sys.exit("error: --genome is required with --calls")
        intact(a) if a.intact else uniform(a)
    else:
        native(a)


if __name__ == "__main__":
    main()
