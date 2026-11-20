#!/usr/bin/env python3
"""Convert one tool's intact LTR-RT calls to a common table, with original sequence names.

Output TSV columns: chrom start end strand id parent lltr_end rltr_start tsd class
  start/end: 0-based half-open element span (genome orientation); lltr_end / rltr_start: end of the left LTR and start
  of the right LTR in the same system ('.' if the tool does not report them); parent: id of the call this element is
  nested in ('.' if none or not reported); class: the tool's classification ('.' if none).
Only full-length (two-LTR) calls are kept: solo LTRs, partial/domain-only elements are dropped.

Formats (per run directory runs/<ds>/<tool>/):
  ltr_peel       <ds>_LTRs_depth*_clean_ltr.tsv   (seq_id chrom:start-end#class, 1-based; nest_status links)
  ltr_peel_hc    same files, confidence == high only (post hoc diagnostic subset; rundir = the LTR_peel run)
  ltrharvest     ltrharvest.tabout                (1-based inclusive; seq-nr = 0-based FASTA index)
  ltr_finder     <ds>.fa.finder.combine.scn       (-w 2 table; Location s-e 1-based inclusive)
  ltr_retriever  <ds>.fa.mod.pass.list            (chr:s..e, 1-based; may be written e..s for minus)
  look4ltrs      out/Rtr/*.rtr                    (0-based; LeftStart..RightEnd; NestOut = host id; solo: RightStart NA)
  hirepeat       *_ltr_output/all_rounds/*_all_rounds_original.bed (seq_index 0-based FASTA index; 0-based)
  hite           out/intact_LTR.list + out/chr_name.map (pass.list format on renamed chr_N)
  ltrpred        *_ltrpred/*_LTRpred_DataSheet.tsv (1-based); none written when LTRpred finds < 3 loci = 0 calls
  dante_ltr      dante_ltr.gff3 transposable_element with Rank in DL, DLT, DLP, DLTP (1-based inclusive)
"""
import argparse
import csv
import glob
import os
import re
import sys


def names(datadir):
    """new_id -> original_id, plus FASTA order (index -> new_id)."""
    m, order = {}, []
    with open(os.path.join(datadir, "names.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            m[r["new_id"]] = r["original_id"]
            order.append(r["new_id"])
    return m, order


def one(pattern):
    hits = sorted(glob.glob(pattern))
    if not hits:
        sys.exit(f"error: no file matches {pattern}")
    return hits


def rec(chrom, s0, e, strand=".", id_=None, parent=".", lend=".", rstart=".", tsd=".", cls="."):
    return dict(chrom=chrom, start=int(s0), end=int(e), strand=strand if strand in "+-" else ".", id=id_,
                parent=parent, lltr_end=lend, rltr_start=rstart, tsd=tsd or ".", cls=cls or ".")


def p_ltr_peel(rd, ds, order, conf=None):
    out = []
    files = [f for f in one(f"{rd}/{ds}_LTRs_depth*_clean_ltr.tsv") if re.search(r"_depth\d+_clean_ltr\.tsv$", f)]
    for fn in files:
        with open(fn) as f:
            rdr = csv.DictReader((l for l in f if l.strip()), delimiter="\t")
            if conf and "confidence" not in (rdr.fieldnames or []):
                sys.exit(f"error: no confidence column in {fn}")
            for r in rdr:
                if conf and r["confidence"] != conf:
                    continue
                sid = r.get("seq_id") or r.get("#seq_id") or r[rdr.fieldnames[0]]
                loc, _, cls = sid.partition("#")
                chrom, _, rng = loc.rpartition(":")
                s, e = (int(x) for x in rng.split("-"))
                ns = r.get("nest_status", "") or ""
                parent = "."
                m = re.search(r"nest-inner:([^;,\s|]+)", ns)
                if m:
                    parent = m.group(1)
                l5 = int(r["ltr5_len"]) if r.get("ltr5_len", "").isdigit() else None
                l3 = int(r["ltr3_len"]) if r.get("ltr3_len", "").isdigit() else None
                out.append(rec(chrom, s - 1, e, r.get("strand", "."), loc, parent,
                               s - 1 + l5 if l5 else ".", e - l3 if l3 else ".", r.get("tsd", "."), cls))
    if conf:                                      # a subset keeps only links whose host is also in it
        ids = {x["id"] for x in out}
        for x in out:
            if x["parent"] not in ids:
                x["parent"] = "."
    return out


def p_ltr_peel_hc(rd, ds, order):
    """Post hoc developer diagnostic, not a primary row: LTR_peel calls with confidence == high (TG..CA ends and an
    exact TSD). rd is the LTR_peel run directory."""
    return p_ltr_peel(rd, ds, order, conf="high")


def p_ltrharvest(rd, ds, order):
    out = []
    with open(f"{rd}/ltrharvest.tabout") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.split()
            s, e, ls_, le, rs = int(c[0]), int(c[1]), int(c[3]), int(c[4]), int(c[6])
            chrom = order[int(c[10])]
            out.append(rec(chrom, s - 1, e, ".", f"{chrom}:{s}-{e}", ".", le, rs - 1))
    return out


def p_ltr_finder(rd, ds, order):
    out = []
    with open(f"{rd}/{ds}.fa.finder.combine.scn") as f:
        for line in f:
            if line.startswith("#") or line.startswith("index") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            chrom, (s, e) = c[1], (int(x) for x in c[2].split("-"))
            l5, l3 = (int(x) for x in c[3].split(","))
            out.append(rec(chrom, s - 1, e, c[12], f"{chrom}:{s}-{e}", ".", s - 1 + l5, e - l3,
                           c[5] if c[5] != "N" else "."))
    return out


def p_passlist(path, rename=None, cls_col=9):
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.split()
            chrom, _, rng = c[0].rpartition(":")
            a, b = (int(x) for x in rng.split(".."))
            s, e = min(a, b), max(a, b)
            if rename:
                chrom = rename.get(chrom, chrom)
            m = re.match(r"IN:(\d+)\.\.(\d+)", c[6]) if len(c) > 6 else None
            lend, rstart = (int(m.group(1)) - 1, int(m.group(2))) if m else (".", ".")
            tsd = c[3].split(":", 1)[1] if c[3].startswith("TSD:") else "."
            strand = c[8] if len(c) > 8 else "."
            cls = c[cls_col] if len(c) > cls_col else "."
            out.append(rec(chrom, s - 1, e, strand, f"{chrom}:{s}-{e}", ".", lend, rstart, tsd, cls))
    return out


def p_ltr_retriever(rd, ds, order):
    hits = glob.glob(f"{rd}/{ds}.fa.mod.pass.list") or glob.glob(f"{rd}/{ds}.fa.pass.list")
    if not hits:
        sys.exit(f"error: no pass.list in {rd}")
    return p_passlist(hits[0])


def p_hite(rd, ds, order):
    ren = {}
    with open(f"{rd}/out/chr_name.map") as f:
        for line in f:
            c = line.split()
            if len(c) >= 2:
                ren[c[0]] = c[1]                  # chr_name.map: renamed (chr_N) <TAB> original (checked on smoke)
    return p_passlist(f"{rd}/out/intact_LTR.list", ren, cls_col=11)


def p_look4ltrs(rd, ds, order):
    out = []
    for fn in one(f"{rd}/out/Rtr/*.rtr"):
        with open(fn) as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for line in f:
                r = dict(zip(hdr, line.rstrip("\n").split("\t")))
                if r.get("RightStart", "NA") == "NA" or r.get("RightEnd", "NA") == "NA":
                    continue                     # solo LTR
                s, e = int(r["LeftStart"]), int(r["RightEnd"])
                hid = f"{r['chrom']}#{r['ID']}"
                nout = (r.get("NestOut") or "NA").strip("{},")
                parent = f"{r['chrom']}#{nout.split(',')[0]}" if nout not in ("NA", "") else "."
                out.append(rec(r["chrom"], s, e, "-" if r.get("RC") == "1" else "+", hid, parent,
                               int(r["LeftEnd"]), int(r["RightStart"])))
    return out


def p_hirepeat(rd, ds, order):
    out = []
    for fn in one(f"{rd}/*_ltr_output/all_rounds/*_all_rounds_original.bed"):
        with open(fn) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                c = line.split()
                chrom = order[int(c[0])]
                out.append(rec(chrom, int(c[1]), int(c[2]), ".", f"r{c[4]}_{c[3]}"))
    return out


def ltrpred_too_few_loci(rd):
    """True if run.log shows LTRpred's designed no-prediction outcome: with fewer than 3 LTRdigest loci its
    read.prediction() reports "weren't enough predicted loci" and LTRpred() returns "No prediction for <genome>"
    without writing a DataSheet. A missing LTRdigest GFF also ends in "No prediction", but that is a failure."""
    few = none = False
    try:
        with open(f"{rd}/run.log", errors="replace") as f:
            for line in f:
                few = few or "weren't enough predicted loci" in line
                none = none or "No prediction for " in line
    except OSError:
        return False
    return few and none


def p_ltrpred(rd, ds, order):
    out = []
    pat = f"{rd}/*_ltrpred/*_LTRpred_DataSheet.tsv"
    if not glob.glob(pat) and ltrpred_too_few_loci(rd):
        print(f"calls_to_bed: LTRpred made no prediction for {ds} (too few loci; run.log): 0 calls", file=sys.stderr)
        return out
    for fn in one(pat):
        with open(fn) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                s, e = int(float(r["start"])), int(float(r["end"]))
                le = r.get("lLTR_end", "NA")
                rs = r.get("rLTR_start", "NA")
                out.append(rec(r["chromosome"], s - 1, e, r.get("strand", "."), f"{r['chromosome']}:{s}-{e}", ".",
                               int(float(le)) if le not in ("NA", "") else ".",
                               int(float(rs)) - 1 if rs not in ("NA", "") else ".",
                               r.get("lTSD_motif", ".")))
    return out


def p_dante_ltr(rd, ds, order):
    out, ltrs = [], {}
    rows = []
    with open(f"{rd}/dante_ltr.gff3") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            at = dict(x.split("=", 1) for x in c[8].split(";") if "=" in x)
            if c[2] == "transposable_element" and at.get("Rank") in ("DL", "DLT", "DLP", "DLTP"):
                rows.append((c, at))
            elif c[2] == "long_terminal_repeat" and "Parent" in at:
                ltrs.setdefault(at["Parent"], []).append((int(c[3]), int(c[4])))
    for c, at in rows:
        s, e = int(c[3]), int(c[4])
        lt = sorted(ltrs.get(at.get("ID"), []))
        lend, rstart = (lt[0][1], lt[1][0] - 1) if len(lt) == 2 else (".", ".")
        out.append(rec(c[0], s - 1, e, c[6], at.get("ID"), ".", lend, rstart, at.get("TSD", "."),
                       at.get("Final_Classification", ".")))
    return out


PARSERS = {k[2:]: v for k, v in globals().items() if k.startswith("p_") and k != "p_passlist"}


def main(argv=None):
    ap = argparse.ArgumentParser(description="Normalize a tool's intact LTR-RT calls (0-based, original names).")
    ap.add_argument("tool", choices=sorted(PARSERS), help="tool whose output to read")
    ap.add_argument("rundir", help="runs/<ds>/<tool> directory")
    ap.add_argument("datadir", help="data/<ds> directory (names.tsv)")
    ap.add_argument("out", help="output TSV")
    ap.add_argument("--keep-ids", action="store_true", help="keep normalized sequence IDs (c0001...)")
    ap.add_argument("-v", "--verbose", action="store_true", help="report counts")
    a = ap.parse_args(argv)
    ds = os.path.basename(os.path.normpath(a.datadir))
    m, order = names(a.datadir)
    calls = PARSERS[a.tool](os.path.abspath(a.rundir), ds, order)
    bad = [x for x in calls if x["chrom"] not in m]
    if bad:
        sys.exit(f"error: {len(bad)} calls on unknown sequences, e.g. {bad[0]['chrom']}")
    n0 = len(calls)
    calls = [x for x in calls if x["end"] > x["start"]]
    seen, uniq = set(), []
    for x in sorted(calls, key=lambda x: (x["chrom"], x["start"], x["end"])):
        k = (x["chrom"], x["start"], x["end"])
        if k not in seen:
            seen.add(k)
            uniq.append(x)
    cols = ["chrom", "start", "end", "strand", "id", "parent", "lltr_end", "rltr_start", "tsd", "class"]
    with open(a.out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for x in uniq:
            y = dict(x, cls=x["cls"])
            if not a.keep_ids:
                y["chrom"] = m[x["chrom"]]
            f.write("\t".join(str(y[c] if c != "class" else y["cls"]) for c in cols) + "\n")
    print(f"calls_to_bed: {a.tool} {ds}: {len(uniq)} calls ({n0 - len(uniq)} empty/duplicate dropped)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
