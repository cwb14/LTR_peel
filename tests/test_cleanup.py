"""The end-of-run clean-up: what goes, what stays, and what is never touched.

The run folders below are empty files named the way a finished run names them
(one genome as in the demo run; four genomes plus their pooled `merged_*`
tables as in a multi-genome run), next to files of the user's own.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from ltrquest import cleanup
from ltrquest.helper_tools import HELPERS

P = "Athal_tair10_chr2_LTRs"
TOOLS = "Athal_tair10_chr2_tools"      # named after the genome, without the _LTRs
INPUTS = ("Athal_tair10_chr2.fa.gz", "Athal.pep.gz")


def touch(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("")


def make(run: Path, names) -> set[str]:
    """Create each path under run; return the names they put directly inside it."""
    for name in names:
        touch(run / name)
    return {name.split("/")[0] for name in names}


def results(run: Path, p: str, depths=(0, 1), plots: bool = True) -> set[str]:
    """One genome's results: what a finished run must leave."""
    names = [f"{p}_all_depth_LTR_cleaned.gff3"]
    names += [f"{p}_depth{d}_clean_ltr.{ext}" for d in depths for ext in ("tsv", "fa")]
    if plots:
        names += [f"{p}_plots/{p}_TEGV.html", f"{p}_plots/{p}_summary.pdf",
                  f"{p}_plots/{p}.chrom.sizes", f"{p}_plots/{p}_combined_ltr.tsv",
                  f"{p}_plots/struct/fam00001_individual_expanded.pdf"]
    return make(run, names)


def working(run: Path, p: str, rounds=(1, 2), depths=(0, 1)) -> set[str]:
    """Everything else a finished run leaves for one genome."""
    names = [f"{p}_r{r}_ltr.{ext}" for r in rounds for ext in ("tsv", "fa")]
    names += [f"{p}_r{r}.work/TEsorter2/{p}_r{r}.cls.tsv" for r in rounds]
    names += [f"{p}_depth{d}_ltr.{ext}" for d in depths for ext in ("tsv", "fa")]
    names += [f"{p}_rejoined.tsv", f"{p}_strand_recovery.tsv", f"{p}_reboundary.tsv",
              f"{p}_fpcheck.log", f"{p}_fpcheck.tandem.tsv",
              f"{p}_all_depth_protein_LTR_cleaned.gff3", f"{p}.detect.json",
              f"{p}_pre_reboundary/{p}_depth0_clean_ltr.tsv",
              f"{p}_FP_masked.fa", f"{p}_FP_masked.fa.fai", f"temp_ltr_2pass_lib_{p}.fa",
              f"temp_ltr_2pass_lib_{p}.fa.fxi"]
    return make(run, names)


def pooled(run: Path, r: str) -> set[str]:
    """The pooled clustering and false-positive call, named after the run prefix."""
    names = [f"{r}_all_ltr.{s}" for s in ("fa", "fa.fai", "tsv", "consensus.fa",
                                         "consensus.fa.fai", "consensus_id0.75_cluster.tsv",
                                         "internal_id0.75_cluster.tsv")]
    names += [f"{r}_fpcheck.{s}" for s in ("log", "tandem.tsv", "family_scores.tsv",
                                          "fp_families.tsv", "fp_LTRs.fa", "fp_structure.pdf")]
    return make(run, names)


def tree(root: Path) -> list[str]:
    return sorted(str(p.relative_to(root)) for p in root.rglob("*"))


@pytest.fixture
def single(tmp_path):
    """A finished one-genome run, everything kept, its inputs in the run folder."""
    run = tmp_path / "run"
    inputs = make(run, [*INPUTS, "run.sh", "run.log"])
    kept = results(run, P)
    gone = (working(run, P) | pooled(run, P)
            | make(run, [f"{TOOLS}/Kmer2LTR/src/kmer2ltr/cli.py"]))
    return run, inputs, kept, gone


def clean(run: Path, *extra: str, protect=INPUTS) -> int:
    argv = ["--indir", str(run), "--prefix", P, "--run-prefix", P]
    for name in protect:
        argv += ["--protect", str(run / name)]
    return cleanup.main([*argv, *extra])


def test_single_genome_keeps_only_primary_outputs(single):
    run, inputs, kept, _gone = single
    assert clean(run, "--tools-dir", str(run / TOOLS)) == 0
    assert set(os.listdir(run)) == inputs | kept
    assert (run / f"{P}_plots" / "struct" / "fam00001_individual_expanded.pdf").is_file()


def test_removable_names_exactly_the_working_files(single):
    run, _inputs, _kept, gone = single
    found = cleanup.removable(str(run), [P], P, [], str(run / TOOLS))
    assert found == sorted(str(run / name) for name in gone)


def test_keep_removes_nothing(single, capsys):
    run, *_ = single
    before = tree(run)
    assert clean(run, "--tools-dir", str(run / TOOLS), "--keep") == 0
    assert tree(run) == before
    assert "Your results" in capsys.readouterr().out


def test_unrelated_files_survive(single):
    run, inputs, kept, _gone = single
    mine = make(run, ["notes.txt", "mydata/reads.bam",
                      # another run's prefix, not passed
                      "other_LTRs_r1_ltr.tsv", "other_LTRs_depth0_ltr.fa", "other_LTRs.detect.json",
                      "other_LTRs_fpcheck.log", "temp_ltr_2pass_lib_other_LTRs.fa",
                      # this prefix, but not a name the run gives its working files
                      f"{P}_r1_ltr.tsv.bak", f"{P}_notes.txt", f"{P}2_r1_ltr.tsv"])
    clean(run, "--tools-dir", str(run / TOOLS))
    assert set(os.listdir(run)) == inputs | kept | mine
    assert (run / "mydata" / "reads.bam").is_file()


def test_protected_inputs_survive(single, tmp_path, monkeypatch):
    run, *_ = single
    elsewhere = tmp_path / "inputs"
    elsewhere.mkdir()
    # through a symlink that names a file in the run folder
    (elsewhere / "genome.fa").symlink_to(run / f"{P}_r2_ltr.fa")
    # a run-named symlink pointing at an input that lives elsewhere
    touch(elsewhere / "proteins.fa")
    (run / f"{P}_depth0_ltr.fa").unlink()
    (run / f"{P}_depth0_ltr.fa").symlink_to(elsewhere / "proteins.fa")

    monkeypatch.chdir(run)                     # a relative path, from the run folder
    cleanup.main(["--indir", str(run), "--prefix", P, "--run-prefix", P,
                  "--protect", f"{P}_r1_ltr.fa", str(elsewhere / "genome.fa"),
                  str(elsewhere / "proteins.fa")])

    assert (run / f"{P}_r1_ltr.fa").is_file()
    assert (run / f"{P}_r2_ltr.fa").is_file()
    assert (run / f"{P}_depth0_ltr.fa").is_symlink()
    assert (elsewhere / "proteins.fa").is_file()
    assert not (run / f"{P}_r1_ltr.tsv").exists()          # their neighbors still go


def test_a_directory_holding_an_input_survives(single):
    run, *_ = single
    touch(run / f"{P}_r1.work" / "genome.fa")
    touch(run / TOOLS / "proteins.fa")
    clean(run, "--tools-dir", str(run / TOOLS),
          protect=(f"{P}_r1.work/genome.fa", f"{TOOLS}/proteins.fa"))
    assert (run / f"{P}_r1.work" / "genome.fa").is_file()
    assert (run / TOOLS / "proteins.fa").is_file()
    assert not (run / f"{P}_r2.work").exists()


def test_a_symlink_is_removed_not_followed(single, tmp_path):
    run, *_ = single
    elsewhere = tmp_path / "elsewhere"
    touch(elsewhere / "TEsorter2" / "x.cls.tsv")
    touch(elsewhere / "table.tsv")
    shutil.rmtree(run / f"{P}_r1.work")
    (run / f"{P}_r1.work").symlink_to(elsewhere, target_is_directory=True)
    (run / f"{P}_r1_ltr.tsv").unlink()
    (run / f"{P}_r1_ltr.tsv").symlink_to(elsewhere / "table.tsv")
    clean(run)
    assert not os.path.lexists(run / f"{P}_r1.work")
    assert not os.path.lexists(run / f"{P}_r1_ltr.tsv")
    assert (elsewhere / "TEsorter2" / "x.cls.tsv").is_file()
    assert (elsewhere / "table.tsv").is_file()


def test_the_staging_dir_is_never_touched(single):
    run, *_ = single
    make(run, [f"{P}_FPstaging/round1/{P}_r1_ltr.tsv"])
    clean(run)
    assert (run / f"{P}_FPstaging" / "round1" / f"{P}_r1_ltr.tsv").is_file()


def test_nothing_is_removed_when_a_genome_has_no_results(single, capsys):
    """With no cleaned tables the unfiltered ones are all the user has."""
    run, *_ = single
    for path in run.glob(f"{P}_depth*_clean_ltr.*"):
        path.unlink()
    before = tree(run)
    assert clean(run, "--tools-dir", str(run / TOOLS)) == 0
    assert tree(run) == before
    assert P in capsys.readouterr().err


def test_multi_genome_cleanup(tmp_path):
    run = tmp_path / "run"
    inputs = make(run, ["A.fa", "A.fa.fai", "B.fa", "B.fa.fai", "memo.md", "ltrquest.log",
                        "ltr_biology/notes.md"])
    kept = results(run, "A_LTRs", depths=(0, 1, 2)) | results(run, "B_LTRs", depths=(0,))
    working(run, "A_LTRs", rounds=(1, 2, 3), depths=(0, 1, 2))
    working(run, "B_LTRs", rounds=(1,), depths=(0,))
    pooled(run, "merged")
    make(run, ["merged_tools/Kmer2LTR/src/kmer2ltr/cli.py"])
    assert cleanup.main(["--indir", str(run), "--prefix", "A_LTRs", "B_LTRs",
                         "--run-prefix", "merged", "--protect", str(run / "A.fa"),
                         str(run / "B.fa"), "--tools-dir", str(run / "merged_tools")]) == 0
    assert set(os.listdir(run)) == inputs | kept


def test_tools_dir_removed_only_when_given_and_inside(single, tmp_path):
    run, *_ = single
    clean(run)                                               # not given
    assert (run / TOOLS).is_dir()

    outside = tmp_path / "shared_tools"
    touch(outside / "TEsorter2" / "pipeline.py")
    clean(run, "--tools-dir", str(outside))                  # outside the run folder
    assert (outside / "TEsorter2" / "pipeline.py").is_file()

    link = run / "linked_tools"                              # inside by name only
    link.symlink_to(outside, target_is_directory=True)
    clean(run, "--tools-dir", str(link))
    assert link.is_symlink() and (outside / "TEsorter2" / "pipeline.py").is_file()

    clean(run, "--tools-dir", str(run))                      # the run folder itself
    assert (run / f"{P}_all_depth_LTR_cleaned.gff3").is_file()

    clean(run, "--tools-dir", str(run / TOOLS))              # given and inside
    assert not (run / TOOLS).exists()


def test_a_tools_dir_holding_anything_else_is_kept(single, capsys):
    """A folder that merely shares the default name -- maize_tools beside maize.fa -- is
    the user's, whatever else is in it."""
    run, *_ = single
    touch(run / TOOLS / "my_scripts" / "plot_me.R")
    clean(run, "--tools-dir", str(run / TOOLS))
    assert (run / TOOLS / "my_scripts" / "plot_me.R").is_file()
    assert (run / TOOLS / "Kmer2LTR" / "src" / "kmer2ltr" / "cli.py").is_file()
    warnings = [line for line in capsys.readouterr().err.splitlines() if TOOLS in line]
    assert len(warnings) == 1 and "my_scripts" in warnings[0]
    assert not (run / f"{P}_r1.work").exists()               # the rest is still cleaned


def test_a_tools_dir_holding_only_helpers_is_removed(single, capsys):
    run, *_ = single
    for helper in ("TEsorter2", "TRF-mod", "sdust", "FASTAN", "minimap2", "miniprot"):
        touch(run / TOOLS / helper / "Makefile")
    clean(run, "--tools-dir", str(run / TOOLS))
    assert not (run / TOOLS).exists()
    assert TOOLS not in capsys.readouterr().err


def test_a_relative_tools_dir_is_read_from_the_current_directory(single, monkeypatch):
    run, *_ = single
    monkeypatch.chdir(run)
    clean(run, "--tools-dir", f"./{TOOLS}")
    assert not (run / TOOLS).exists()


def test_every_helper_fetched_into_a_tools_dir_is_known():
    """The names the fetchers give their helpers at the top of a tools dir, read off their
    code, so a new helper cannot make every tools dir look like the user's."""
    found = set()
    for path in Path(cleanup.__file__).parent.rglob("*.py"):
        text = path.read_text()
        found |= set(re.findall(r'tools_dir\)?\s*/\s*"([^"/]+)"', text))
        found |= set(re.findall(r'_resolve_helper\(\s*"([^"]+)"', text))
    assert len(found) >= 7 and found <= HELPERS, found - HELPERS


@pytest.mark.parametrize("names", [
    ["--prefix", "", "--run-prefix", P],
    ["--prefix", P, "", "--run-prefix", P],
    ["--prefix", P, "--run-prefix", ""],
])
def test_empty_prefixes_are_refused(single, names):
    """An empty run prefix would make the pooled patterns match any *_all_ltr.* file."""
    run, *_ = single
    before = tree(run)
    with pytest.raises(SystemExit) as excinfo:
        cleanup.main(["--indir", str(run), *names])
    assert excinfo.value.code == 2
    assert tree(run) == before


def test_summary_lists_each_genome_and_tegv(tmp_path):
    run = tmp_path / "run"
    results(run, "A_LTRs", depths=(0, 1, 2))
    results(run, "B_LTRs", depths=(0,), plots=False)
    text = cleanup.results_summary(str(run), ["A_LTRs", "B_LTRs"])
    lines = text.splitlines()
    assert f"Your results for A_LTRs (in {run}):" in lines
    assert f"Your results for B_LTRs (in {run}):" in lines
    wanted = {
        "A_LTRs_all_depth_LTR_cleaned.gff3": "all elements, for a genome browser",
        "A_LTRs_depth0_clean_ltr.tsv / .fa": "elements with nothing nested inside",
        "A_LTRs_depth1_clean_ltr.tsv / .fa": "elements with 1 layer nested inside",
        "A_LTRs_depth2_clean_ltr.tsv / .fa": "elements with 2 layers nested inside",
        "A_LTRs_plots/A_LTRs_TEGV.html": "open in a web browser",
        "B_LTRs_all_depth_LTR_cleaned.gff3": "all elements, for a genome browser",
        "B_LTRs_depth0_clean_ltr.tsv / .fa": "elements with nothing nested inside",
    }
    for name, meaning in wanted.items():
        hits = [line for line in lines if name in line]
        assert len(hits) == 1, name
        assert re.fullmatch(rf"  {re.escape(name)} +{re.escape(meaning)}", hits[0]), hits[0]
    assert "plots" not in text[text.index("Your results for B_LTRs"):]
    assert lines[-1] == "How to read them: https://github.com/cwb14/LTRquest#outputs"


def test_primary_outputs_in_reading_order(tmp_path):
    run = tmp_path / "run"
    results(run, "A_LTRs", depths=(0, 1, 2, 10))
    got = cleanup.primary_outputs(str(run), "A_LTRs")
    assert all(os.path.isabs(path) for path in got)
    assert [os.path.basename(path) for path in got] == [
        "A_LTRs_all_depth_LTR_cleaned.gff3",
        *[f"A_LTRs_depth{d}_clean_ltr.{ext}" for d in (0, 1, 2, 10) for ext in ("tsv", "fa")],
        "A_LTRs_plots"]


def test_runs_as_a_module(single):
    run, *_ = single
    result = subprocess.run([sys.executable, "-m", "ltrquest.cleanup", "--indir", str(run),
                             "--prefix", P, "--run-prefix", P, "--keep"],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert f"Your results for {P} (in {run}):" in result.stdout


# ---------------------------------------------------------------------------
# The driver's side: flags, records, and where the call sits
# ---------------------------------------------------------------------------

def run_driver(driver, *args, **kw):
    return subprocess.run(["bash", str(driver), *args], capture_output=True, text=True, **kw)


def orchestrator(driver) -> str:
    text = driver.read_text()
    return text[text.index("run_fp_orchestrator() {"):text.index("# Parse args")]


def cleanup_block(driver) -> str:
    """Where the orchestrator builds the clean-up's arguments, up to the call."""
    body = orchestrator(driver)
    return body[body.index("local -a cleanup_args=("):
                body.index('"${CLEANUP[@]}" "${cleanup_args[@]}"')]


def test_driver_calls_cleanup_after_plots(driver):
    body = orchestrator(driver)
    plots = body.index('bash "${SELF_DIR}/plots.sh"')
    call = body.index('"${CLEANUP[@]}" "${cleanup_args[@]}"')
    assert plots < call < body.index('echo "Done."')
    for arg in ('--indir "$final_dir"', '--prefix "${OUT_PREFIXES[@]}"',
                '--run-prefix "$RUN_PREFIX"', '--protect "${abs_genomes[@]}"',
                '--protect "$abs_proteins"', '--tools-dir "$TOOLS_DIR"', "--keep"):
        assert arg in body[plots:call], arg
    # nothing that builds the arguments can stop a finished run under set -e
    assert "$(" not in cleanup_block(driver)
    # a clean-up that fails warns; the run it follows has already succeeded
    assert re.search(r'"\$\{CLEANUP\[@\]\}" "\$\{cleanup_args\[@\]\}" \\\n\s+\|\| echo "WARNING',
                     body)
    assert re.search(r'^CLEANUP=\(\s+"\$PY" -m ltrquest\.cleanup\s+\)$', driver.read_text(),
                     re.MULTILINE)


def test_a_tools_dir_named_or_already_there_is_kept(driver):
    """Whether the default tools dir was already there is noted before anything can make
    it (the Kmer2LTR check at start-up may clone into it). One the user named, or one
    already there, is shielded at the end; only one this run made is passed for removal."""
    text = driver.read_text()
    derived = text.index('TOOLS_DIR="${LTRQUEST_TOOLS_DIR:-./${RUN_PREFIX}_tools}"')
    noted = text.index('[[ -e "$TOOLS_DIR" || -L "$TOOLS_DIR" ]] && TOOLS_DIR_PREEXISTED=true')
    assert derived < noted < text.index("kmer2ltr.api(sys.argv[1])")
    block = cleanup_block(driver)
    branch = block[block.index('if [[ "$USER_TOOLS_DIR" == true || '
                               '"$TOOLS_DIR_PREEXISTED" == true ]]; then'):]
    assert branch.index('--protect "$TOOLS_DIR"') < branch.index("else") \
        < branch.index('--tools-dir "$TOOLS_DIR"')


def test_a_failed_step_keeps_the_working_files(driver):
    """A finishing step that failed and let the run go on can only be redone from the
    working files, so then the clean-up removes nothing."""
    text = driver.read_text()
    for warning in ("strand recovery failed for", "could not re-orient the _clean_ depth FASTAs",
                    "re-boundarying failed; the calls are kept as detected"):
        at = text.index(warning)
        assert "note_failed_step" in text[at:at + 400], warning
    body = orchestrator(driver)
    at = body.index("plotting stage reported failures")
    assert "failed_steps+=(" in body[at:at + 300]
    # noted in the staging dir, which promotion never empties into the results, and read
    # back before that dir is removed
    assert 'FAILED_STEPS="${staging}/' in body
    assert body.index('mapfile -t failed_steps < "$FAILED_STEPS"') \
        < body.index('rm -rf "$staging"\n')
    block = cleanup_block(driver)
    failed = block.index('elif (( ${#failed_steps[@]} > 0 )); then')
    assert block.index("cleanup_args+=( --keep )", failed) < block.index("else", failed)


def test_records_only_with_keep_intermediates(driver):
    body = orchestrator(driver)
    write = body.index('"${RECORD[@]}" write')
    guard = body.rindex('if [[ "$KEEP_INTERMEDIATES" == true ]]; then', 0, write)
    assert write < body.index("\n  fi\n", guard)


def test_the_worker_removes_its_pass2_library_and_index(driver):
    """TEsorter2 indexes the pass-2 library beside it (pyfastx's .fxi); left there it
    would be promoted into the results."""
    assert 'rm -f "$temp_lib" "$temp_lib".* 2>/dev/null || true' in driver.read_text()


def test_the_driver_makes_its_scratch_files_in_the_run_folder(driver):
    """Not in the system temp folder, which is small on many cluster nodes."""
    calls = re.findall(r"\$\(mktemp\b[^)]*\)", driver.read_text())
    assert calls
    for call in calls:
        assert '"./' in call, call


def test_flags_parse(driver, tmp_path):
    tools = tmp_path / "x"
    result = run_driver(driver, "--keep-intermediates", "--tools-dir", str(tools), "--help")
    assert result.returncode == 0
    assert not tools.exists()                      # asking for help creates nothing


# Both stop at the package check (no interpreter), after the tools dir is settled.
STOP = {"PATH": "/usr/bin:/bin", "LTRQUEST_PYTHON": "/nonexistent/python"}
QUICK = ("--strand-recovery", "off", "--no-reboundary")


def test_the_orchestrator_resolves_and_creates_the_tools_dir(driver, toy_genome, tmp_path):
    result = run_driver(driver, "--genome", str(toy_genome), *QUICK, "--tools-dir", "my_tools",
                        env=STOP, cwd=tmp_path)
    assert result.returncode != 0 and "not importable" in result.stderr
    assert f"Tools directory: {(tmp_path / 'my_tools').resolve()}\n" in result.stdout
    assert (tmp_path / "my_tools").is_dir()


def test_a_worker_ignores_the_tools_dir_value(driver, toy_genome, tmp_path):
    """Workers replay the argv from a staging dir; the orchestrator's export is what counts."""
    result = run_driver(driver, "--genome", str(toy_genome), *QUICK, "--tools-dir", "my_tools",
                        env={**STOP, "_LTRQUEST_WORKER": "1"}, cwd=tmp_path)
    assert "Tools directory: ./toy_tools\n" in result.stdout
    assert not (tmp_path / "my_tools").exists()
