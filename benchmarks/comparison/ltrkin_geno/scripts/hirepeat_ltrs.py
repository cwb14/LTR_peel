#!/usr/bin/env python3
"""Add LTR edges to HiRepeat calls (which report none) from a Kmer2LTR annot.tsv.

Kmer2LTR seq_id is <normalized chrom>:<start1>-<end>. Only status 'pass' or 'weak_pair' counts, with
0 < ltr5_end < ltr3_start - 1 < end - start. Then lltr_end = start + ltr5_end and rltr_start = start + ltr3_start - 1
(0-based, as in adjudicate.py). Calls without a paired entry keep None edges and are counted as unpaired."""
import argparse
import sys

from lg_io import die, read_calls, read_names, read_tsv, write_calls


def add_kmer2ltr_ltrs(calls, annot_path, names):
    """Set lltr_end/rltr_start on calls from Kmer2LTR annot.tsv. names: normalized -> original chrom.
    Returns (calls, number of calls without paired LTR edges)."""
    rev = {v: k for k, v in names.items()}
    rows = read_tsv(annot_path)
    if rows and not {"seq_id", "status", "ltr5_end", "ltr3_start"} <= set(rows[0]):
        die(f"{annot_path}: needs columns seq_id, status, ltr5_end, ltr3_start")
    inner = {}
    for r in rows:
        if r["status"] not in ("pass", "weak_pair"):
            continue
        try:
            c, rng = r["seq_id"].rsplit(":", 1)
            s1, e = rng.split("-")
            key = (c, int(s1) - 1, int(e))
        except ValueError:
            die(f"{annot_path}: malformed seq_id {r['seq_id']!r} (expected chrom:start-end)")
        try:
            a, b = int(r["ltr5_end"]), int(r["ltr3_start"]) - 1
        except ValueError:
            continue
        if 0 < a < b < key[2] - key[1]:
            inner[key] = (a, b)
    n_unpaired = 0
    for c in calls:
        hit = inner.get((rev.get(c["chrom"], c["chrom"]), c["start"], c["end"]))
        if hit:
            c["lltr_end"], c["rltr_start"] = c["start"] + hit[0], c["start"] + hit[1]
        else:
            n_unpaired += 1
    return calls, n_unpaired


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True, help="HiRepeat calls TSV (calls_to_bed.py format)")
    ap.add_argument("--annot", required=True, help="Kmer2LTR annot.tsv for the same calls")
    ap.add_argument("--names", required=True, help="dataset directory holding names.tsv (new_id, original_id)")
    ap.add_argument("--out", required=True, help="output calls TSV with lltr_end/rltr_start filled")
    ap.add_argument("-v", "--verbose", action="store_true", help="also report paired counts and unpaired examples")
    a = ap.parse_args()
    calls = read_calls(a.calls, "hirepeat")
    out, n = add_kmer2ltr_ltrs(calls, a.annot, read_names(a.names))
    write_calls(out, a.out)
    print(f"hirepeat_ltrs: {n} of {len(out)} calls unpaired (no LTR edges)", file=sys.stderr)
    if a.verbose:
        print(f"hirepeat_ltrs: {len(out) - n} paired; first unpaired: "
              f"{[(c['chrom'], c['start'], c['end']) for c in out if c['lltr_end'] is None][:5]}", file=sys.stderr)


if __name__ == "__main__":
    main()
