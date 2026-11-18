#!/usr/bin/env python3
"""Convert a calls TSV (0-based half-open) to GFF3 (1-based inclusive) for ltrkin and hap_truth.py.

Each element is an LTR_retrotransposon (ID=<tool>_<idx:06d>, Name=<call id>, tsd=<SEQ> when known) with two
long_terminal_repeat children and, when the TSD is known, two target_site_duplication features.
Mode 'ltrkin' skips calls without LTR edges (skipped_no_ltr) and calls with inconsistent edges, i.e. not
start < lltr_end <= rltr_start < end (skipped_inconsistent). Mode 'truth' gives both kinds pseudo-LTRs (first and
last PSEUDO_LTR of the span, at least 1 bp; counted in pseudo_ltr, and inconsistent ones also in inconsistent_ltr);
a span too short for two pseudo-LTRs is skipped (skipped_too_short). Each skipped or repaired call is warned about
on stderr (tool, idx, coordinates), never fatal."""
import argparse
import sys

from lg_config import PSEUDO_LTR
from lg_io import die, read_calls


def gff_id(tool, idx):
    return f"{tool}_{idx:06d}"


def _esc(s):
    """Percent-encode GFF3 attribute-reserved characters (; = , & %) and tabs/newlines."""
    return "".join(f"%{ord(ch):02X}" if ch in ";=,&%\t\n\r" else ch for ch in str(s))


def write_gff(calls, out, tool, mode):
    """Write GFF3; returns counts {written, skipped_no_ltr, pseudo_ltr, inconsistent_ltr, skipped_inconsistent,
    skipped_too_short}. Records sorted by chrom and start. Per-call problems warn on stderr and never abort."""
    if mode not in ("ltrkin", "truth"):
        die(f"unknown mode {mode!r} (use ltrkin or truth)")
    n = dict(written=0, skipped_no_ltr=0, pseudo_ltr=0, inconsistent_ltr=0, skipped_inconsistent=0, skipped_too_short=0)

    def warn(c, msg):
        print(f"warning: {tool} idx={c.get('idx')} {c['chrom']}:{c['start']}-{c['end']}: {msg}", file=sys.stderr)

    recs = []
    for c in calls:
        s, e = c["start"], c["end"]
        le, rs = c.get("lltr_end"), c.get("rltr_start")
        if le is None or rs is None:
            if mode == "ltrkin":
                n["skipped_no_ltr"] += 1
                continue
            le = rs = None
        elif not (s < le <= rs < e):
            if mode == "ltrkin":
                warn(c, f"inconsistent LTR edges {le}/{rs}; skipped")
                n["skipped_inconsistent"] += 1
                continue
            warn(c, f"inconsistent LTR edges {le}/{rs}; using pseudo-LTRs")
            n["inconsistent_ltr"] += 1
            le = rs = None
        if le is None:
            p = max(1, int((e - s) * PSEUDO_LTR))
            le, rs = s + p, e - p
            if not (s < le <= rs < e):
                warn(c, "span too short for two pseudo-LTRs; skipped")
                n["skipped_too_short"] += 1
                n["inconsistent_ltr"] -= 1 if (c.get("lltr_end") is not None and c.get("rltr_start") is not None) else 0
                continue
            n["pseudo_ltr"] += 1
        gid = c.get("gff_id") or gff_id(tool, c["idx"])
        strand = c.get("strand") if c.get("strand") in ("+", "-") else "."
        t = len(c["tsd"]) if c.get("tsd") else 0
        attrs = f"ID={gid}"
        if c.get("id") not in (None, "", "."):
            attrs += f";Name={_esc(c['id'])}"
        if t:
            attrs += f";tsd={c['tsd']}"
        ch = c["chrom"]
        f = [(ch, s, "LTR_retrotransposon", s + 1, e, strand, attrs),
             (ch, s, "long_terminal_repeat", s + 1, le, strand, f"ID={gid}_L;Parent={gid}"),
             (ch, s, "long_terminal_repeat", rs + 1, e, strand, f"ID={gid}_R;Parent={gid}")]
        if t:
            f += [(ch, s, "target_site_duplication", s - t + 1, s, strand, f"ID={gid}_T5;Parent={gid}"),
                  (ch, s, "target_site_duplication", e + 1, e + t, strand, f"ID={gid}_T3;Parent={gid}")]
        recs.append((ch, s, -e, gid, f))
        n["written"] += 1
    recs.sort(key=lambda r: r[:4])
    with open(out, "w") as fh:
        fh.write("##gff-version 3\n")
        for r in recs:
            for ch, _, ty, a, b, st, at in r[4]:
                fh.write(f"{ch}\t{tool}\t{ty}\t{a}\t{b}\t.\t{st}\t.\t{at}\n")
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True, help="calls TSV (calls_to_bed.py format, 0-based half-open)")
    ap.add_argument("--tool", required=True, help="tool name: GFF3 source column and ID prefix")
    ap.add_argument("--mode", required=True, choices=["ltrkin", "truth"],
                    help="ltrkin: skip calls with missing or inconsistent LTR edges; truth: give them pseudo-LTRs "
                         "(spans too short for two are skipped); all with a stderr warning per call")
    ap.add_argument("--out", required=True, help="output GFF3 path")
    ap.add_argument("-v", "--verbose", action="store_true", help="also report the input call count")
    a = ap.parse_args()
    calls = read_calls(a.calls, a.tool)
    n = write_gff(calls, a.out, a.tool, a.mode)
    print(f"calls_to_gff: {a.tool} {a.mode}: written={n['written']} skipped_no_ltr={n['skipped_no_ltr']} "
          f"pseudo_ltr={n['pseudo_ltr']} inconsistent_ltr={n['inconsistent_ltr']} "
          f"skipped_inconsistent={n['skipped_inconsistent']} skipped_too_short={n['skipped_too_short']}", file=sys.stderr)
    if a.verbose:
        print(f"calls_to_gff: {len(calls)} input calls", file=sys.stderr)


if __name__ == "__main__":
    main()
