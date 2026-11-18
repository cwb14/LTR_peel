#!/usr/bin/env python3
"""Background set B and per-tool control call sets for ltrkin.

ltrkin learns its genotype prior, read error, depth calibration and boundary tolerance from the whole annotation it
is given, so a tool's genotypes at real loci depend on its other calls. The control holds those parameters nearly
equal across tools: every tool's control annotation is the same background B plus that tool's own calls at the
pooled loci.

Subcommands:
  build    B = calls of all primary tools with consistent LTR edges that overlap no candidate locus (its candidate
           span, widened to its true span when it has true ends; input candidates_te.tsv), clustered at
           reciprocal overlap >= RO (union-find); one representative per cluster (the median by length, then tool
           order, then start). Written as background.tsv (calls columns + tool idx gff_id; ids bg_000000.. in
           genome order).
  control  GFF3 (ltrkin mode) of B plus the tool's calls with consistent LTR edges that overlap a pooled locus's
           match span by >= 1 bp. Tool calls keep their own ID <tool>_<idx:06d>.
Coordinates are 0-based half-open (GFF3 is written 1-based by write_gff)."""
import argparse
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# Appended, never prepended: the main benchmark's scripts/ has a queue.py that would shadow stdlib queue.
sys.path.append(os.path.join(os.environ["BENCH_ROOT"], "scripts"))
from score_elements import Index  # noqa: E402

from calls_to_gff import gff_id, write_gff  # noqa: E402
from lg_config import RO, TOOLS  # noqa: E402
from lg_io import BG_COLS, die, final_calls_path, read_background, read_calls, write_tsv  # noqa: E402
from pool import read_candidates, read_pool  # noqa: E402


def has_ltrs(c):
    """True when both LTR edges are present and consistent: start < lltr_end <= rltr_start < end."""
    le, rs = c.get("lltr_end"), c.get("rltr_start")
    return le is not None and rs is not None and c["start"] < le <= rs < c["end"]


def _find(p, i):
    while p[i] != i:
        p[i] = p[p[i]]
        i = p[i]
    return i


def exclusion_spans(candidates):
    """Span each candidate locus excludes from B: the members' union [start, end), widened to
    [min(start, true_start), max(end, true_end)) when the locus has true ends (the true element can extend past the
    union)."""
    out = []
    for l in candidates:
        s, e = l["start"], l["end"]
        ts, te = l.get("true_start"), l.get("true_end")
        if ts is not None and te is not None:
            s, e = min(s, ts), max(e, te)
        out.append(dict(chrom=l["chrom"], start=s, end=e))
    return out


def build_background(calls_by_tool, candidates, stats=None):
    """Background B (list of call dicts + gff_id=bg_<k:06d>, k in genome order: chrom, start, end).

    Uses calls of TOOLS only, with usable LTR edges, overlapping no candidate locus by >= 1 bp: neither its
    candidate span nor, when known, its true span (exclusion_spans). If stats (dict) is given it gets n_input,
    n_excluded, n_no_ltr, n_clusters."""
    allc = [c for t in TOOLS for c in calls_by_tool.get(t, [])]
    no_ltr = [c for c in allc if not has_ltrs(c)]
    ok = [c for c in allc if has_ltrs(c)]
    cidx = Index(exclusion_spans(candidates)) if candidates else None
    kept = [c for c in ok if cidx is None or not cidx.overlaps(c["chrom"], c["start"], c["end"])]
    idx = Index(kept) if kept else None
    par = list(range(len(kept)))
    for i, c in enumerate(kept):
        li = c["end"] - c["start"]
        for j in idx.overlaps(c["chrom"], c["start"], c["end"]):
            if j <= i:
                continue
            d = kept[j]
            ov = min(c["end"], d["end"]) - max(c["start"], d["start"])
            if min(ov / li, ov / (d["end"] - d["start"])) >= RO:
                par[_find(par, i)] = _find(par, j)
    groups = {}
    for i in range(len(kept)):
        groups.setdefault(_find(par, i), []).append(kept[i])
    reps = []
    for g in groups.values():
        g.sort(key=lambda c: (c["end"] - c["start"], TOOLS.index(c["tool"]), c["start"]))
        reps.append(g[(len(g) - 1) // 2])
    reps.sort(key=lambda c: (c["chrom"], c["start"], c["end"], TOOLS.index(c["tool"])))
    bg = [dict(c, gff_id=f"bg_{k:06d}") for k, c in enumerate(reps)]
    if stats is not None:
        stats.update(n_input=len(allc), n_excluded=len(ok) - len(kept), n_no_ltr=len(no_ltr),
                     n_clusters=len(groups))
    return bg


def control_calls(bg, tool_calls, pool):
    """B plus the tool's calls with usable LTR edges that overlap any pooled locus's match span by >= 1 bp.

    Tool calls get gff_id = <tool>_<idx:06d> (their own ID)."""
    spans = [dict(chrom=p["chrom"], start=p["match_start"], end=p["match_end"]) for p in pool]
    pidx = Index(spans) if spans else None
    out = list(bg)
    for c in tool_calls:
        if has_ltrs(c) and pidx is not None and pidx.overlaps(c["chrom"], c["start"], c["end"]):
            out.append(dict(c, gff_id=gff_id(c["tool"], c["idx"])))
    return out


def cmd_build(a):
    for p in (a.calls_dir, a.candidates):
        if not os.path.exists(p):
            die(f"not found: {p}")
    by_tool = {}
    for t in TOOLS:
        p = final_calls_path(a.calls_dir, t)
        if not os.path.isfile(p):
            print(f"warning: {a.ds}: no final calls for {t} ({p}); not in the background", file=sys.stderr)
            continue
        by_tool[t] = read_calls(p, t)
        if a.verbose:
            print(f"background: {t}: {len(by_tool[t])} calls", file=sys.stderr)
    if not by_tool:
        die(f"no tool calls found under {a.calls_dir}")
    loci = read_candidates(a.candidates)
    if loci and "te_flag" not in loci[0]:
        die(f"{a.candidates}: no true-end columns; give candidates_te.tsv from `pool.py true-ends`, so that calls "
            "at a locus's true span are excluded too")
    st = {}
    bg = build_background(by_tool, loci, st)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    write_tsv([dict(b, **{"class": b.get("cls")}) for b in bg], a.out, BG_COLS)
    print(f"background: {a.ds}: input_calls={st['n_input']} excluded_at_candidates={st['n_excluded']} "
          f"without_usable_ltrs={st['n_no_ltr']} clusters={st['n_clusters']} B={len(bg)}", file=sys.stderr)


def cmd_control(a):
    bg = read_background(a.background)
    calls = read_calls(a.calls, a.tool)
    pool = read_pool(a.pool)
    out = control_calls(bg, calls, pool)
    n_tool = len(out) - len(bg)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    n = write_gff(out, a.out, a.tool, "ltrkin")
    skipped = n["skipped_no_ltr"] + n["skipped_inconsistent"]
    print(f"background: control {a.tool}: B={len(bg)} tool_calls_at_pool={n_tool} written={n['written']} "
          f"skipped_no_usable_ltrs={skipped}", file=sys.stderr)
    if a.verbose:
        pooled = Index([dict(chrom=p["chrom"], start=p["match_start"], end=p["match_end"]) for p in pool]) if pool else None
        n_ov = sum(1 for c in calls if pooled and pooled.overlaps(c["chrom"], c["start"], c["end"]))
        print(f"background: {len(calls)} {a.tool} calls, {n_ov} overlap the pool, {n_ov - n_tool} of those "
              f"skipped for missing or inconsistent LTR edges", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="build background.tsv for one dataset")
    b.add_argument("--ds", required=True, help="dataset name (labels messages)")
    b.add_argument("--calls-dir", required=True, help="directory with final calls of every tool "
                   "(<tool>.tsv; hirepeat.ltr.tsv for hirepeat)")
    b.add_argument("--candidates", required=True, help="candidate loci with true ends (candidates_te.tsv from "
                   "`pool.py true-ends`); calls overlapping a locus's candidate span or true span are excluded")
    b.add_argument("--out", required=True, help="output background.tsv")
    b.add_argument("-v", "--verbose", action="store_true", help="also report per-tool call counts")
    c = sub.add_parser("control", help="write a tool's control GFF3 (B + its calls at pooled loci)")
    c.add_argument("--background", required=True, help="background.tsv from build")
    c.add_argument("--calls", required=True, help="the tool's final calls TSV")
    c.add_argument("--tool", required=True, help="tool name (GFF3 source column, ID prefix)")
    c.add_argument("--pool", required=True, help="pool.tsv of confirmed insertion loci")
    c.add_argument("--out", required=True, help="output GFF3")
    c.add_argument("-v", "--verbose", action="store_true", help="also report how many pooled-locus calls were skipped")
    a = ap.parse_args()
    {"build": cmd_build, "control": cmd_control}[a.cmd](a)


if __name__ == "__main__":
    main()
