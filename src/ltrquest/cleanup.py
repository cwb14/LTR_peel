"""Clear a finished run's working files out of the run folder, and say where its results are.

A run leaves its results -- per genome the GFF3, the depth tables and FASTAs with
the false positives removed, and the plots -- beside everything it built on the
way: round tables and work directories, unfiltered depth tables, pooled cluster
tables, logs and sidecars, the reuse record, and the helper tools it fetched.
The driver calls this once every result is in place, and never after a failure:
it removes the working files unless --keep is given, then prints the results.

Removed: only the names this run gives its working files, directly inside the
run folder, and the tools directory the driver names when it lies inside that
folder and holds nothing but the helpers LTRquest fetches there (anything else
in it is kept, with a warning). A symlink is removed, never what it points to.
Never removed: the
results; an input passed with --protect, or any directory holding one; the
FP-correction staging directory, which the driver owns; and anything at all
when a genome has no results, since its unfiltered tables are then all it has.

Usage:
  python -m ltrquest.cleanup --indir D --prefix P [P ...] --run-prefix R
                             [--protect PATH ...] [--tools-dir T] [--keep]
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import sys

from . import record
from .helper_tools import HELPERS

DOCS = "https://github.com/cwb14/LTRquest#outputs"

# After a genome's '<prefix>_': what detection leaves for the pooled stages (the
# names a reuse record lists), then what the pooled stages and the finishing
# passes leave beside the results.
_WORKING = (record._OUTPUTS
            + r"|reboundary\.tsv|pre_reboundary|fpcheck\..*|all_ltr\..*"
            + r"|all_depth_protein_LTR_cleaned\.gff3|FP_masked\.fa(\.fai)?")
# After the run prefix's '<prefix>_': the pooled clustering and false-positive call.
_POOLED = r"all_ltr\..*|fpcheck\..*"
# After any of the run's prefixes: the results, kept whatever else a name matches.
_RESULTS = r"all_depth_LTR_cleaned\.gff3|depth\d+_clean_ltr\.(tsv|fa)|plots"
_STAGING = "_FPstaging"


def _named(prefix: str, alternatives: str) -> re.Pattern[str]:
    return re.compile(rf"{re.escape(prefix)}_({alternatives})")


def _depth_files(indir: str, prefix: str) -> list[tuple[int, str, str]]:
    """(depth, 'tsv' | 'fa', name) of each cleaned depth table and FASTA, by depth."""
    pattern = re.compile(rf"{re.escape(prefix)}_depth(\d+)_clean_ltr\.(tsv|fa)")
    found = []
    for name in os.listdir(indir):
        m = pattern.fullmatch(name)
        if m and os.path.exists(os.path.join(indir, name)):
            found.append((int(m.group(1)), m.group(2), name))
    return sorted(found, key=lambda f: (f[0], f[1] != "tsv"))


def primary_outputs(indir: str, prefix: str) -> list[str]:
    """One genome's results as they stand: the GFF3, the cleaned depth tables and
    FASTAs by depth (each table before its FASTA), then the plots folder."""
    indir = os.path.abspath(indir)
    found = []
    gff3 = os.path.join(indir, f"{prefix}_all_depth_LTR_cleaned.gff3")
    if os.path.isfile(gff3):
        found.append(gff3)
    found += [os.path.join(indir, name) for _depth, _ext, name in _depth_files(indir, prefix)]
    plots = os.path.join(indir, f"{prefix}_plots")
    if os.path.isdir(plots):
        found.append(plots)
    return found


def missing_results(indir: str, prefixes: list[str]) -> list[str]:
    """The prefixes with no GFF3 or no cleaned depth table in indir."""
    return [p for p in prefixes
            if not os.path.isfile(os.path.join(indir, f"{p}_all_depth_LTR_cleaned.gff3"))
            or not any(ext == "tsv" for _depth, ext, _name in _depth_files(indir, p))]


def _holds(path: str, protected: str) -> bool:
    """True when removing path would remove protected (a real path): the same file,
    or a directory above it."""
    real = os.path.realpath(path)
    return os.path.commonpath([real, protected]) == real


def _inside(path: str, directory: str) -> bool:
    """True when path's real location lies below directory's, and is not directory itself."""
    real, top = os.path.realpath(path), os.path.realpath(directory)
    return real != top and os.path.commonpath([real, top]) == top


def _strangers(tools_dir: str | None, indir: str) -> list[str] | None:
    """For a tools dir inside indir, what its top holds besides the helpers LTRquest
    fetches there; None when there is no such tools dir. A folder that merely shares
    the default name is the user's, and anything else in it says so."""
    if not tools_dir or not os.path.isdir(tools_dir) or not _inside(tools_dir, indir):
        return None
    try:
        return sorted(set(os.listdir(tools_dir)) - HELPERS)
    except OSError as exc:
        return [f"(unreadable: {exc.strerror})"]


def removable(indir: str, prefixes: list[str], run_prefix: str, protect: list[str],
              tools_dir: str | None) -> list[str]:
    """Sorted absolute paths this run may delete. Empty when a genome has no results."""
    indir = os.path.abspath(indir)
    if missing_results(indir, prefixes):
        return []
    working = [_named(p, _WORKING) for p in prefixes] + [_named(run_prefix, _POOLED)]
    # The reuse record, and the pass-2 library with the index TEsorter2 builds beside it.
    exact = ({p + record.SUFFIX for p in prefixes}
             | {f"temp_ltr_2pass_lib_{p}.fa{ext}" for p in prefixes for ext in ("", ".fxi")})
    candidates = {os.path.join(indir, name) for name in os.listdir(indir)
                  if name in exact or any(w.fullmatch(name) for w in working)}
    if _strangers(tools_dir, indir) == []:
        candidates.add(os.path.abspath(tools_dir))

    results = [_named(p, _RESULTS) for p in {*prefixes, run_prefix}]
    shields = [os.path.realpath(p) for p in protect]
    return sorted(path for path in candidates
                  if not os.path.basename(path).endswith(_STAGING)
                  and not any(r.fullmatch(os.path.basename(path)) for r in results)
                  and not any(_holds(path, s) for s in shields))


def _remove(path: str) -> None:
    """Remove a file, a directory tree, or a symlink itself (never what it points to)."""
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    else:
        os.remove(path)


def _nested(depth: int) -> str:
    if depth == 0:
        return "elements with nothing nested inside"
    return f"elements with {depth} layer{'s' if depth > 1 else ''} nested inside"


def results_summary(indir: str, prefixes: list[str]) -> str:
    """The "Your results" block: per genome, one line for each result and what it holds."""
    indir = os.path.abspath(indir)
    blocks = []
    for prefix in prefixes:
        rows = []
        gff3 = f"{prefix}_all_depth_LTR_cleaned.gff3"
        if os.path.isfile(os.path.join(indir, gff3)):
            rows.append((gff3, "all elements, for a genome browser"))
        kinds: dict[int, list[str]] = {}
        for depth, ext, _name in _depth_files(indir, prefix):
            kinds.setdefault(depth, []).append(ext)
        for depth, exts in sorted(kinds.items()):
            stem = f"{prefix}_depth{depth}_clean_ltr"
            rows.append((f"{stem}.tsv / .fa" if len(exts) == 2 else f"{stem}.{exts[0]}",
                         _nested(depth)))
        tegv = f"{prefix}_plots/{prefix}_TEGV.html"
        if os.path.isfile(os.path.join(indir, tegv)):
            rows.append((tegv, "open in a web browser"))
        elif os.path.isdir(os.path.join(indir, f"{prefix}_plots")):
            rows.append((f"{prefix}_plots/", "plots"))
        blocks.append((prefix, rows))

    width = max((len(name) for _prefix, rows in blocks for name, _ in rows), default=0) + 2
    lines = []
    for prefix, rows in blocks:
        lines.append(f"Your results for {prefix} (in {indir}):")
        lines += [f"  {name:<{width}}{meaning}" for name, meaning in rows] \
            or ["  (no result files found)"]
    lines.append(f"How to read them: {DOCS}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="ltrquest-cleanup",
        description="Remove a finished run's working files, keeping each genome's GFF3, "
                    "cleaned depth tables and FASTAs, and plots; then say where those are.")
    parser.add_argument("--indir", required=True, help="the run folder")
    parser.add_argument("--prefix", nargs="+", required=True,
                        help="every genome's output prefix, e.g. Athal_LTRs")
    parser.add_argument("--run-prefix", required=True,
                        help="the prefix naming the pooled cluster tables and FP files")
    parser.add_argument("--protect", nargs="+", action="extend", default=[], metavar="PATH",
                        help="an input file: never removed, nor any directory holding it")
    parser.add_argument("--tools-dir", default=None,
                        help="the helper-tools directory the run made; removed only when "
                             "it lies inside --indir and holds nothing but the helpers "
                             "LTRquest fetches")
    parser.add_argument("--keep", action="store_true",
                        help="remove nothing; only say where the results are")
    args = parser.parse_args(argv)
    if not os.path.isdir(args.indir):
        parser.error(f"--indir is not a directory: {args.indir}")
    # An empty prefix would turn '<prefix>_all_ltr.*' into any '_all_ltr.*'.
    if not all(args.prefix) or not args.run_prefix:
        parser.error("--prefix and --run-prefix need non-empty names")

    if not args.keep:
        missing = missing_results(args.indir, args.prefix)
        strangers = _strangers(args.tools_dir, args.indir)
        if missing:
            print(f"WARNING: no GFF3 or cleaned depth tables for {', '.join(missing)}; "
                  f"every working file stays in {os.path.abspath(args.indir)}",
                  file=sys.stderr)
        elif strangers:
            shown = ", ".join(strangers[:3]) + (f" and {len(strangers) - 3} more"
                                                if len(strangers) > 3 else "")
            print(f"WARNING: kept {os.path.abspath(args.tools_dir)}: besides the helper "
                  f"tools LTRquest fetches it holds {shown}", file=sys.stderr)
        for path in removable(args.indir, args.prefix, args.run_prefix, args.protect,
                              args.tools_dir):
            try:
                _remove(path)
            except FileNotFoundError:
                pass
            except OSError as exc:
                print(f"WARNING: could not remove {path}: {exc}", file=sys.stderr)
    print(results_summary(args.indir, args.prefix))
    return 0


if __name__ == "__main__":
    sys.exit(main())
