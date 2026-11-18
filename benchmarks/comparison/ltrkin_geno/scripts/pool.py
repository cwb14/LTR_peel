"""Pool of real insertion polymorphisms: nominate candidate loci from per-call truth tables, find each locus's
true element ends, and pool the loci whose insertion is confirmed.

Subcommands:
  candidates   read truth/<ds>/<tool>.tsv (hap_truth.py output), nominate calls, merge into candidate loci.
  true-ends    realign each locus's reference and other-haplotype segments for tool-independent true ends
               (candidates_te.tsv) and write the true elements as GFF3 for hap_truth.py.
  final        pool = loci whose true element is confirmed E; a locus without true ends (te_flag != ok) enters
               only through a member call that is itself confirmed E.

Coordinates inside this module are 0-based half-open. read_truth converts hap_truth.py output (start 1-based,
end inclusive, o_left/o_right 1-based) at the reader.
"""
import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pysam  # noqa: E402

from calls_to_gff import write_gff  # noqa: E402
from lg_config import (MAX_LOCUS, MIN_DUP, MIN_INTERNAL, NOM_ABS, NOM_REL, O_PAD, SEG_FLANK, SPREAD,  # noqa: E402
                       TOOLS, TSD_ASSUMED)
from lg_io import die, write_tsv  # noqa: E402
from true_ends import true_ends  # noqa: E402

TRUTH_REQUIRED = ["element_id", "chrom", "start", "end", "length", "ltr5_len", "ltr3_len", "tsd_len", "o_chrom",
                  "o_strand", "o_left", "o_right", "delta", "note", "state", "truth_source"]
CAND_COLS = ["locus_id", "chrom", "start", "end", "n_members", "members", "median_tool", "median_element",
             "median_length", "median_ltr5", "median_ltr3", "o_chrom", "o_strand", "o_left", "o_right", "n_pass2",
             "flags"]
TE_COLS = ["true_start", "true_end", "dup_len", "dup_seq", "window", "te_flag"]   # candidates_te.tsv = CAND + TE
TE_FLAGS = ["ok", "no_clean_gap", "too_long", "short_dup"]
POOL_COLS = ["locus_id", "chrom", "match_start", "match_end", "true_start", "true_end", "dup_len", "dup_seq",
             "window", "te_flag", "basis", "n_members", "n_pass2", "members", "flags"]
POOL_INT = ["match_start", "match_end", "true_start", "true_end", "dup_len", "window", "n_members", "n_pass2"]
BASES = ["true_confirmed", "member_confirmed"]


def warn(msg):
    print(f"warning: {msg}", file=sys.stderr)


def _int(v):
    return int(v)


def _int_or_none(v):
    return None if v in (None, "", ".") else int(v)


def read_truth(path, tool):
    """hap_truth.py table -> list of dicts (0-based half-open; o_left/o_right 0-based or None; delta int or None).

    Fails fast on a missing file/column; warns and skips a row with a non-integer start/end/length/tsd_len."""
    if not os.path.isfile(path):
        die(f"truth file not found: {path}")
    with open(path, newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        missing = [c for c in TRUTH_REQUIRED if c not in (rd.fieldnames or [])]
        if missing:
            die(f"{path}: missing required column(s): {', '.join(missing)}")
        rows = []
        for i, d in enumerate(rd, start=2):
            try:
                row = dict(tool=tool, element_id=d["element_id"], chrom=d["chrom"], start=_int(d["start"]) - 1,
                           end=_int(d["end"]), length=_int(d["length"]), ltr5_len=_int_or_none(d["ltr5_len"]),
                           ltr3_len=_int_or_none(d["ltr3_len"]), tsd_len=_int(d["tsd_len"]),
                           o_chrom=d["o_chrom"], o_strand=d["o_strand"],
                           o_left=None, o_right=None, delta=None, note=d["note"], state=d["state"],
                           truth_source=d["truth_source"])
                for k in ("o_left", "o_right"):
                    v = _int_or_none(d[k])
                    row[k] = None if v is None or v < 0 else v - 1   # -1 = absent; else 1-based -> 0-based
                row["delta"] = _int_or_none(d["delta"])
            except (ValueError, TypeError):
                warn(f"{path}: row {i}: non-integer value, skipped")
                continue
            rows.append(row)
    return rows


def _length(row):
    return row["length"] if row.get("length") is not None else row["end"] - row["start"]


def nominates(row):
    """True if the call's own length change fits an empty pre-insertion site on the other haplotype."""
    if row.get("o_chrom") in (None, ".", "") or row.get("delta") is None:
        return False
    L = _length(row)
    return abs(row["delta"] + L + row["tsd_len"]) <= NOM_ABS + NOM_REL * L


def _build_locus(chrom, group, ds, n):
    members = sorted(group, key=lambda r: (r["start"], r["end"], TOOLS.index(r["tool"]), r["element_id"]))
    by_len = sorted(group, key=lambda r: (_length(r), TOOLS.index(r["tool"]), r["start"]))
    med = by_len[(len(by_len) - 1) // 2]
    lens = [_length(r) for r in group]
    flags = ["spread"] if (min(lens) > 0 and (max(lens) - min(lens)) / min(lens) > SPREAD) else []
    return dict(locus_id=f"{ds}_L{n:05d}", chrom=chrom, start=min(r["start"] for r in group),
                end=max(r["end"] for r in group),
                members=[(r["tool"], r["element_id"], r["start"], r["end"]) for r in members],
                median_tool=med["tool"], median_element=med["element_id"], median_length=_length(med),
                median_ltr5=med.get("ltr5_len"), median_ltr3=med.get("ltr3_len"), o_chrom=med["o_chrom"],
                o_strand=med["o_strand"], o_left=med["o_left"], o_right=med["o_right"],
                n_pass2=sum(r.get("note") == "remapped" for r in group), flags=flags)


def candidate_loci(rows_by_tool, ds):
    """Merge nominating calls of the primary TOOLS into loci (union of overlapping spans, per chrom)."""
    by_chrom = defaultdict(list)
    for tool, rows in rows_by_tool.items():
        if tool not in TOOLS:
            continue
        for r in rows:
            if nominates(r):
                by_chrom[r["chrom"]].append(dict(r, tool=tool))
    loci, n = [], 0
    groups = []
    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda r: (r["start"], r["end"]))
        cur, cur_end = [], None
        for r in rows:
            if cur and r["start"] < cur_end:
                cur.append(r)
                cur_end = max(cur_end, r["end"])
            else:
                if cur:
                    groups.append((chrom, cur))
                cur, cur_end = [r], r["end"]
        if cur:
            groups.append((chrom, cur))
    groups.sort(key=lambda g: (g[0], min(r["start"] for r in g[1])))
    for chrom, g in groups:
        n += 1
        loci.append(_build_locus(chrom, g, ds, n))
    return loci


def _fmt(v):
    return "." if v is None else str(v)


def _fmt_members(members):
    return ",".join(f"{t}:{e}:{s}:{en}" for t, e, s, en in members)


def write_candidates(loci, path, te=False):
    """Candidates TSV (CAND_COLS); te=True adds TE_COLS (candidates_te.tsv). None -> '.'."""
    cols = CAND_COLS + (TE_COLS if te else [])
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for l in loci:
            row = dict(l, n_members=len(l["members"]), members=_fmt_members(l["members"]),
                       flags=",".join(l["flags"]) or ".")
            w.writerow([_fmt(row[c]) for c in cols])


def _read_te(d, path, i):
    """TE_COLS of one candidates_te row; fails fast on an unknown te_flag or ends that contradict it."""
    flag = d["te_flag"]
    if flag not in TE_FLAGS:
        die(f"{path}: row {i}: unknown te_flag {flag!r} (expected one of {', '.join(TE_FLAGS)})")
    te = dict(true_start=_int_or_none(d["true_start"]), true_end=_int_or_none(d["true_end"]),
              dup_len=_int_or_none(d["dup_len"]), dup_seq=None if d["dup_seq"] in ("", ".") else d["dup_seq"],
              window=_int_or_none(d["window"]), te_flag=flag)
    if (flag == "ok") != (te["true_start"] is not None and te["true_end"] is not None):
        die(f"{path}: row {i}: te_flag {flag} contradicts true ends {d['true_start']}/{d['true_end']}")
    return te


def read_candidates(path):
    """Candidates TSV -> loci dicts with the keys of candidate_loci (members as tuples); a candidates_te TSV also
    gives TE_COLS (ints or None; dup_seq str or None)."""
    if not os.path.isfile(path):
        die(f"candidates file not found: {path}")
    out = []
    with open(path, newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        fields = rd.fieldnames or []
        missing = [c for c in CAND_COLS if c not in fields]
        if missing:
            die(f"{path}: missing column(s): {', '.join(missing)}")
        has_te = [c in fields for c in TE_COLS]
        if any(has_te) and not all(has_te):
            die(f"{path}: incomplete true-end columns; missing: "
                f"{', '.join(c for c, h in zip(TE_COLS, has_te) if not h)}")
        for i, d in enumerate(rd, start=2):
            members = []
            for m in d["members"].split(","):
                tool, rest = m.split(":", 1)
                elem, s, e = rest.rsplit(":", 2)       # element ids could contain ':'
                members.append((tool, elem, int(s), int(e)))
            out.append(dict(
                locus_id=d["locus_id"], chrom=d["chrom"], start=int(d["start"]), end=int(d["end"]),
                members=members, median_tool=d["median_tool"], median_element=d["median_element"],
                median_length=int(d["median_length"]), median_ltr5=_int_or_none(d["median_ltr5"]),
                median_ltr3=_int_or_none(d["median_ltr3"]), o_chrom=d["o_chrom"], o_strand=d["o_strand"],
                o_left=_int_or_none(d["o_left"]), o_right=_int_or_none(d["o_right"]), n_pass2=int(d["n_pass2"]),
                flags=[] if d["flags"] == "." else d["flags"].split(",")))
            if all(has_te):
                out[-1].update(_read_te(d, path, i))
    return out


# ---- true ends per locus, true elements, final pool ----

_COMP = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn")


def revcomp(s):
    """Reverse complement; keeps case (soft-masking), complements IUPAC codes, N stays N."""
    return s.translate(_COMP)[::-1]


def segments(locus, ref, other):
    """(r_seq, r_offset, o_seq) for true_ends, or None when end - start > MAX_LOCUS (nothing is fetched).

    r_seq = ref[start - SEG_FLANK, end + SEG_FLANK), r_offset its start; o_seq = other[o_left - SEG_FLANK - O_PAD,
    o_right + 1 + SEG_FLANK + O_PAD) around the median member's lifted anchors, reverse-complemented when o_strand
    is '-' so that it reads like r_seq. Both are clipped to their contig (0-based half-open, pysam). o_seq is ''
    when the locus has no other-haplotype anchors (o_left or o_right None, e.g. a remapped row's projection)."""
    s, e = locus["start"], locus["end"]
    if e - s > MAX_LOCUS:
        return None
    r0 = max(0, s - SEG_FLANK)
    r_seq = ref.fetch(locus["chrom"], r0, min(ref.get_reference_length(locus["chrom"]), e + SEG_FLANK))
    ol, orr = locus.get("o_left"), locus.get("o_right")
    if ol is None or orr is None:
        return r_seq, r0, ""
    o0 = max(0, ol - SEG_FLANK - O_PAD)
    o1 = min(other.get_reference_length(locus["o_chrom"]), orr + 1 + SEG_FLANK + O_PAD)
    o_seq = other.fetch(locus["o_chrom"], o0, o1) if o1 > o0 else ""
    return r_seq, r0, revcomp(o_seq) if locus["o_strand"] == "-" else o_seq


def locus_true_ends(locus, ref, other):
    """A copy of the locus plus true_start true_end (genome 0-based half-open, None unless te_flag is ok) dup_len
    dup_seq window te_flag, from true_ends on the locus's segments.

    The locus span (the members' union) is passed to true_ends as the hint that identifies the element gap, so an
    unrelated long indel elsewhere in the flanks does not make the locus NA (only the near rule applies to it).

    te_flag: ok; too_long (locus > MAX_LOCUS, not realigned); no_clean_gap (true_ends found no clean element gap,
    or there is no other-haplotype segment); short_dup (dup_len < MIN_DUP: the ends are not pinned by a TSD, so
    true ends are None but dup_len, dup_seq and window are kept). window = max(0, dup_len - TSD_ASSUMED), the
    microhomology window. Only ok and short_dup have a dup_len; dup_seq is None when dup_len is 0."""
    na = dict(true_start=None, true_end=None, dup_len=None, dup_seq=None, window=None)
    seg = segments(locus, ref, other)
    if seg is None:
        return dict(locus, **na, te_flag="too_long")
    r_seq, r_off, o_seq = seg
    hint = (locus["start"] - r_off, locus["end"] - r_off)
    te = true_ends(r_seq, o_seq, hint) if r_seq and o_seq else None
    if te is None:
        return dict(locus, **na, te_flag="no_clean_gap")
    dup = dict(dup_len=te.dup_len, dup_seq=te.dup_seq or None, window=max(0, te.dup_len - TSD_ASSUMED))
    if te.dup_len < MIN_DUP:
        return dict(locus, true_start=None, true_end=None, **dup, te_flag="short_dup")
    return dict(locus, true_start=r_off + te.start, true_end=r_off + te.end, **dup, te_flag="ok")


def write_true_elements(loci, out_gff):
    """GFF3 of the true elements for hap_truth.py: one per te_flag == ok locus, ID = locus_id, at the true ends.
    Returns the number written.

    LTRs: the median member's LTR lengths clipped to the true span; pseudo-LTRs (write_gff, mode truth) when either
    is unknown or together they leave < MIN_INTERNAL bp of internal region. They only affect hap_truth's solo
    window. tsd = dup_seq only when 3 <= dup_len <= 8 (ACGT), the range hap_truth.py reads (else it assumes 5)."""
    calls = []
    for i, l in enumerate(loci):
        if l["te_flag"] != "ok":
            continue
        s, e = l["true_start"], l["true_end"]
        l5, l3 = (min(x, e - s) if x else 0 for x in (l.get("median_ltr5"), l.get("median_ltr3")))
        le, rs = (s + l5, e - l3) if l5 > 0 and l3 > 0 and (e - s) - l5 - l3 >= MIN_INTERNAL else (None, None)
        seq = l["dup_seq"] or ""
        tsd = seq if 3 <= l["dup_len"] <= 8 and set(seq) <= set("ACGT") else None
        calls.append(dict(gff_id=l["locus_id"], idx=i, chrom=l["chrom"], start=s, end=e, strand=".",
                          lltr_end=le, rltr_start=rs, tsd=tsd))
    return write_gff(calls, out_gff, "true_ends", "truth")["written"]


def median_member(members):
    """Median member (tool, element_id, start, end) by the key (length, TOOLS index, start), index (n-1)//2."""
    by_len = sorted(members, key=lambda m: (m[3] - m[2], TOOLS.index(m[0]), m[2]))
    return by_len[(len(by_len) - 1) // 2]


def _confirmed_e(row):
    return row is not None and row.get("state") == "E" and row.get("truth_source") == "confirmed"


def final_pool(loci, true_rows, member_rows):
    """Pooled loci as POOL_COLS dicts, in input order (members/flags as written in the TSV).

    basis true_confirmed: te_flag ok and the true element (true_rows[locus_id]) is state E, truth_source
    confirmed. basis member_confirmed: te_flag != ok and some member's own row (member_rows[(tool, element_id)])
    is confirmed E. An ok locus never pools through its members. match_start/match_end = the true span, else
    the median member's span."""
    out = []
    for l in loci:
        if l["te_flag"] == "ok":
            basis = "true_confirmed" if _confirmed_e(true_rows.get(l["locus_id"])) else None
            ms, me = l["true_start"], l["true_end"]
        else:
            hit = any(_confirmed_e(member_rows.get((t, e))) for t, e, _, _ in l["members"])
            basis = "member_confirmed" if hit else None
            _, _, ms, me = median_member(l["members"])
        if basis is None:
            continue
        out.append(dict(locus_id=l["locus_id"], chrom=l["chrom"], match_start=ms, match_end=me,
                        true_start=l["true_start"], true_end=l["true_end"], dup_len=l["dup_len"],
                        dup_seq=l["dup_seq"], window=l["window"], te_flag=l["te_flag"], basis=basis,
                        n_members=len(l["members"]), n_pass2=l["n_pass2"], members=_fmt_members(l["members"]),
                        flags=",".join(l["flags"]) or None))
    return out


def write_pool(rows, path):
    """Pool TSV with POOL_COLS; None -> '.'."""
    write_tsv(rows, path, POOL_COLS)


def read_pool(path):
    """Pool TSV -> list of dicts: '.' -> None, POOL_INT columns parsed as int; members ('tool:element:start:end',
    comma-separated) and flags (comma-separated) stay strings as written; extra columns are kept. Fails fast on a
    missing file or column or a non-integer value."""
    if not os.path.isfile(path):
        die(f"pool file not found: {path}")
    out = []
    with open(path, newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        missing = [c for c in POOL_COLS if c not in (rd.fieldnames or [])]
        if missing:
            die(f"{path}: missing column(s): {', '.join(missing)}")
        for i, d in enumerate(rd, start=2):
            row = {k: None if v in ("", ".") else v for k, v in d.items()}
            for c in POOL_INT:
                if row[c] is not None:
                    try:
                        row[c] = int(row[c])
                    except ValueError:
                        die(f"{path}: row {i}: column {c!r} is not an integer: {row[c]!r}")
            out.append(row)
    return out


def cmd_candidates(a):
    if not os.path.isdir(a.truth_dir):
        die(f"truth directory not found: {a.truth_dir}")
    rows_by_tool, missing = {}, []
    for tool in TOOLS:
        p = os.path.join(a.truth_dir, f"{tool}.tsv")
        if os.path.isfile(p):
            rows_by_tool[tool] = read_truth(p, tool)
            if a.verbose:
                print(f"{tool}: {len(rows_by_tool[tool])} rows", file=sys.stderr)
        else:
            missing.append(tool)
    if missing:
        warn(f"missing truth file for {', '.join(missing)} (tool did not finish?); continuing without them")
    if not rows_by_tool:
        die(f"no truth files for any tool in {a.truth_dir}")
    loci = candidate_loci(rows_by_tool, a.ds)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    write_candidates(loci, a.out)
    nom = {t: sum(nominates(r) for r in rows) for t, rows in rows_by_tool.items()}
    n_pass2 = sum(r["note"] == "remapped" and nominates(r) for rows in rows_by_tool.values() for r in rows)
    n_spread = sum("spread" in l["flags"] for l in loci)
    print(f"{a.ds}: {len(loci)} loci; nominating calls " + " ".join(f"{t}={n}" for t, n in nom.items())
          + f"; pass2={n_pass2}; spread={n_spread}; missing tools={','.join(missing) or 'none'}", file=sys.stderr)
    if a.verbose:
        sizes = defaultdict(int)
        for l in loci:
            sizes[len(l["members"])] += 1
        print("loci by n_members: " + " ".join(f"{k}:{v}" for k, v in sorted(sizes.items())), file=sys.stderr)


def _open_fasta(path):
    """pysam FastaFile; fails fast if the FASTA or its .fai is missing (pysam would write an index beside it)."""
    if not os.path.isfile(path):
        die(f"FASTA not found: {path}")
    if not os.path.isfile(path + ".fai"):
        die(f"FASTA index not found: {path}.fai (index it with samtools faidx)")
    return pysam.FastaFile(path)


def cmd_true_ends(a):
    loci = read_candidates(a.candidates)
    bad = [l["locus_id"] for l in loci if not l["locus_id"].startswith(f"{a.ds}_")]
    if bad:
        die(f"{a.candidates}: {len(bad)} locus id(s) are not from dataset {a.ds!r}, e.g. {bad[0]}")
    ref, other = _open_fasta(a.ref), _open_fasta(a.other)
    miss = sorted({l["chrom"] for l in loci} - set(ref.references))
    if miss:
        die(f"{len(miss)} candidate chrom(s) not in --ref {a.ref}, e.g. {miss[0]} (wrong reference?)")
    miss = sorted({l["o_chrom"] for l in loci if l["o_left"] is not None and l["o_right"] is not None}
                  - set(other.references))
    if miss:
        die(f"{len(miss)} other-haplotype contig(s) not in --other {a.other}, e.g. {miss[0]} (wrong FASTA?)")
    out = []
    for i, l in enumerate(loci, start=1):
        out.append(locus_true_ends(l, ref, other))
        if a.verbose and i % 200 == 0:
            print(f"true-ends: {i}/{len(loci)} loci", file=sys.stderr)
    for p in (a.out_tsv, a.out_gff):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    write_candidates(out, a.out_tsv, te=True)
    n_el = write_true_elements(out, a.out_gff)
    n = Counter(l["te_flag"] for l in out)
    print(f"{a.ds}: {len(out)} loci; te_flag " + " ".join(f"{k}={n[k]}" for k in TE_FLAGS)
          + f"; true elements written={n_el}", file=sys.stderr)
    if n_el != n["ok"]:
        warn(f"{n['ok'] - n_el} ok locus/loci not written to {a.out_gff} (see warnings above)")
    if a.verbose:
        dups = Counter(l["dup_len"] for l in out if l["te_flag"] in ("ok", "short_dup"))
        print("dup_len (ok + short_dup): " + " ".join(f"{k}:{v}" for k, v in sorted(dups.items()))
              + f"; ok loci with window > 0: {sum(l['window'] > 0 for l in out if l['te_flag'] == 'ok')}"
              + f"; loci without other-haplotype anchors: {sum(l['o_left'] is None or l['o_right'] is None for l in loci)}",
              file=sys.stderr)


def cmd_final(a):
    loci = read_candidates(a.candidates_te)
    if any("te_flag" not in l for l in loci):
        die(f"{a.candidates_te}: no true-end columns ({' '.join(TE_COLS)}); make it with `pool.py true-ends`")
    true_rows = {r["element_id"]: r for r in read_truth(a.true_truth, "true_element")}
    if not os.path.isdir(a.truth_dir):
        die(f"truth directory not found: {a.truth_dir}")
    member_rows, present = {}, []
    for tool in TOOLS:
        p = os.path.join(a.truth_dir, f"{tool}.tsv")
        if os.path.isfile(p):
            present.append(tool)
            member_rows.update(((tool, r["element_id"]), r) for r in read_truth(p, tool))
    if not present:
        die(f"no truth files for any tool in {a.truth_dir}")
    no_true = [l["locus_id"] for l in loci if l["te_flag"] == "ok" and l["locus_id"] not in true_rows]
    if no_true:
        warn(f"{len(no_true)} true element(s) missing from {a.true_truth}, e.g. {no_true[0]}; they cannot pool")
    mem = [(t, e) for l in loci for t, e, _, _ in l["members"]]
    no_mem = [m for m in mem if m not in member_rows]
    if no_mem:
        warn(f"{len(no_mem)} of {len(mem)} member call(s) missing from {a.truth_dir}, e.g. {':'.join(no_mem[0])}; "
             "they cannot confirm a locus")
    pool = final_pool(loci, true_rows, member_rows)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    write_pool(pool, a.out)
    n = Counter(p["basis"] for p in pool)
    print(f"pool {a.out}: {len(pool)} of {len(loci)} loci pooled; basis "
          + " ".join(f"{b}={n[b]}" for b in BASES), file=sys.stderr)
    if a.verbose:
        tf = Counter(p["te_flag"] for p in pool)
        st = Counter(f"{r['state']}/{r['truth_source']}" for l in loci if l["te_flag"] == "ok"
                     for r in [true_rows.get(l["locus_id"])] if r is not None)
        print("pooled by te_flag: " + " ".join(f"{k}={tf[k]}" for k in TE_FLAGS)
              + f"; pooled with pass-2 nominations: {sum(p['n_pass2'] > 0 for p in pool)}"
              + "; true-element state/truth_source: " + " ".join(f"{k}:{v}" for k, v in sorted(st.items())),
              file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Insertion-polymorphism pool: nominate candidate loci from truth tables, "
                                             "find their true ends, pool the confirmed ones.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("candidates", help="nominate calls and merge them into candidate loci")
    c.add_argument("--ds", required=True, help="dataset name, e.g. camara_hap2 (prefix of locus ids)")
    c.add_argument("--truth-dir", required=True, help="directory with <tool>.tsv from hap_truth.py")
    c.add_argument("--out", required=True, help="output candidates TSV")
    c.add_argument("-v", "--verbose", action="store_true", help="per-tool and per-locus-size detail")
    c.set_defaults(func=cmd_candidates)
    t = sub.add_parser("true-ends", help="tool-independent true element ends per candidate locus, and the true "
                                         "elements as GFF3 for hap_truth.py")
    t.add_argument("--ds", required=True, help="dataset name; every locus id must start with '<ds>_'")
    t.add_argument("--candidates", required=True, help="candidates TSV from `pool.py candidates`")
    t.add_argument("--ref", required=True, help="reference haplotype FASTA the calls are on (.fai required)")
    t.add_argument("--other", required=True, help="other haplotype FASTA, target of the lift (.fai required)")
    t.add_argument("--out-tsv", required=True, help="output candidates_te TSV: the candidates columns plus "
                                                    "true_start true_end dup_len dup_seq window te_flag")
    t.add_argument("--out-gff", required=True, help="output GFF3 with one element per te_flag ok locus")
    t.add_argument("-v", "--verbose", action="store_true", help="progress every 200 loci and dup_len counts")
    t.set_defaults(func=cmd_true_ends)
    f = sub.add_parser("final", help="pool the loci whose true element (or, without true ends, a member call) is "
                                     "confirmed E")
    f.add_argument("--candidates-te", required=True, help="candidates_te TSV from `pool.py true-ends`")
    f.add_argument("--true-truth", required=True, help="hap_truth.py table of the true-element GFF3")
    f.add_argument("--truth-dir", required=True, help="directory with each tool's <tool>.tsv from hap_truth.py "
                                                      "(the member calls' own rows)")
    f.add_argument("--out", required=True, help="output pool TSV")
    f.add_argument("-v", "--verbose", action="store_true", help="pooled counts by te_flag and true-element states")
    f.set_defaults(func=cmd_final)
    a = ap.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
