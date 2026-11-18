"""Shared I/O for the ltrkin_geno scripts: calls TSV (0-based half-open), names map, FASTA lengths, generic TSV."""
import csv
import os
import re
import sys

CALL_COLS = ["chrom", "start", "end", "strand", "id", "parent", "lltr_end", "rltr_start", "tsd", "class"]
_TSD = re.compile(r"[ACGT]{3,8}")


def die(msg):
    sys.exit(f"error: {msg}")


def read_tsv(path):
    """Read a tab-separated file with header into a list of dicts (all values str)."""
    if not os.path.isfile(path):
        die(f"file not found: {path}")
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(rows, path, cols):
    """Write rows (dicts) with the given column order; None -> '.'."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow(["." if r.get(c) is None else r[c] for c in cols])


def _int_or_none(v, path, row_no, col):
    if v is None or v in (".", "", "NA"):
        return None
    try:
        return int(v)
    except ValueError:
        die(f"{path}: row {row_no}: column {col!r} is not an integer: {v!r}")


def read_calls(path, tool):
    """Calls TSV (calls_to_bed.py format) -> list of dicts; coordinates 0-based half-open.

    Keys: chrom start end strand id parent lltr_end rltr_start tsd cls tool idx. tsd is upper-cased and kept only
    if it is 3-8 ACGT letters, else None. idx is the 0-based row number."""
    if not os.path.isfile(path):
        die(f"calls file not found: {path}")
    with open(path, newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        need = [c for c in CALL_COLS if c != "class"] + ["class"]
        missing = [c for c in need if c not in (rd.fieldnames or [])]
        if missing:
            die(f"{path}: missing column(s) {missing}; expected {CALL_COLS}")
        out = []
        for i, r in enumerate(rd):
            s, e = _int_or_none(r["start"], path, i, "start"), _int_or_none(r["end"], path, i, "end")
            if s is None or e is None or e <= s:
                die(f"{path}: row {i}: bad start/end {r['start']!r}/{r['end']!r}")
            tsd = (r["tsd"] or "").upper()
            out.append(dict(chrom=r["chrom"], start=s, end=e, strand=r["strand"], id=r["id"], parent=r["parent"],
                            lltr_end=_int_or_none(r["lltr_end"], path, i, "lltr_end"),
                            rltr_start=_int_or_none(r["rltr_start"], path, i, "rltr_start"),
                            tsd=tsd if _TSD.fullmatch(tsd) else None, cls=r["class"], tool=tool, idx=i))
    return out


def write_calls(calls, path):
    """Write calls in calls_to_bed.py column order; None -> '.'."""
    rows = [dict(c, **{"class": c.get("cls")}) for c in calls]
    write_tsv(rows, path, CALL_COLS)


def read_names(datadir):
    """<datadir>/names.tsv (new_id, original_id) -> {normalized: original}."""
    p = os.path.join(datadir, "names.tsv")
    rows = read_tsv(p)
    if rows and not {"new_id", "original_id"} <= set(rows[0]):
        die(f"{p}: needs columns new_id and original_id")
    return {r["new_id"]: r["original_id"] for r in rows}


def fai_lengths(fasta):
    """{contig: length} from <fasta>.fai; if absent, from one streaming pass over the FASTA (not indexed)."""
    fai = fasta + ".fai"
    if os.path.isfile(fai):
        with open(fai) as f:
            return {x[0]: int(x[1]) for x in (l.rstrip("\n").split("\t") for l in f if l.strip())}
    if not os.path.isfile(fasta):
        die(f"FASTA not found: {fasta}")
    out, name = {}, None
    with open(fasta) as f:
        for l in f:
            if l.startswith(">"):
                name = l[1:].split()[0]
                out[name] = 0
            elif name is not None:
                out[name] += len(l.strip())
    return out


def final_calls_path(calls_dir, tool):
    """Final calls TSV of a tool: <calls_dir>/hirepeat.ltr.tsv for hirepeat (Kmer2LTR edges added), else
    <calls_dir>/<tool>.tsv (also for ltrquest_hc)."""
    return os.path.join(calls_dir, "hirepeat.ltr.tsv" if tool == "hirepeat" else f"{tool}.tsv")


BG_COLS = CALL_COLS + ["tool", "idx", "gff_id"]


def read_background(path):
    """background.tsv -> list of dicts with the read_calls keys plus gff_id; tool and idx are the original call's.
    Fails fast on a missing file or column; warns and skips a row with a non-integer coordinate."""
    if not os.path.isfile(path):
        die(f"background file not found: {path}")
    out = []
    with open(path, newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        missing = [c for c in BG_COLS if c not in (rd.fieldnames or [])]
        if missing:
            die(f"{path}: missing column(s) {missing}; expected {BG_COLS}")
        for i, r in enumerate(rd, start=2):
            try:
                s, e, idx = int(r["start"]), int(r["end"]), int(r["idx"])
                le = None if r["lltr_end"] in (".", "", "NA") else int(r["lltr_end"])
                rs = None if r["rltr_start"] in (".", "", "NA") else int(r["rltr_start"])
            except ValueError:
                print(f"warning: {path}: row {i}: non-integer value, skipped", file=sys.stderr)
                continue
            tsd = (r["tsd"] or "").upper()
            out.append(dict(chrom=r["chrom"], start=s, end=e, strand=r["strand"], id=r["id"], parent=r["parent"],
                            lltr_end=le, rltr_start=rs, tsd=tsd if _TSD.fullmatch(tsd) else None, cls=r["class"],
                            tool=r["tool"], idx=idx, gff_id=r["gff_id"]))
    return out
