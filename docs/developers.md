# Developer notes

For people changing LTRquest. Users want the [README](../README.md); the output
formats are in [outputs.md](outputs.md).

## How it fits together

```
ltrquest (driver, src/ltrquest/scripts/ltrquest.sh)
  |
  |-- round 1..N --+- ltrquest.detect   LTRharvest + LTR_finder (+ FastLTR on capped chunks)
  |                |                    -> sdust filter -> Kmer2LTR -> TEsorter2
  |                +- ltrquest.mask     mask this round's hits -> next round's genome
  |
  |-- ltrquest.reconcile                pool rounds, resolve containment -> depth buckets
  |-- Kmer2LTR + ltrquest.flag_fp       cluster into families, purge false-positive families
  |                                     and tandem-array calls (rDNA, satellites)
  |-- ltrquest.recover_strand           give strand to unstranded elements   (off by default)
  |-- ltrquest.annotate                 add strand, family, confidence columns
  |-- ltrquest.recover_strand           re-orient the _clean_ FASTAs to match (off by default)
  |-- ltrquest.reboundary               re-bound calls that lack a TSD       (on by default)
  |-- ltrquest.gff3                     pooled GFF3 (after annotate when re-boundarying is off)
  |-- ltrquest-plots                    structure PDFs, summary PDF, TEGV browser
  +-- ltrquest.cleanup                  remove working files (unless --keep-intermediates)
```

The stages that run by default and their off switches are in
[outputs.md](outputs.md#1-what-runs-by-default-and-how-to-turn-it-off). Every
stage is also its own command (`ltrquest-reconcile --help`, `ltrquest-gff3
--help`, and so on), so any of them can be re-run alone. The Nextflow pipeline
runs the same stages as modules; see [nextflow.md](nextflow.md).

`ltrquest.cleanup` runs last, once every result is in place, and never after a
failed step: when plotting, strand recovery or re-boundarying failed, everything is
kept so the step can be redone. It removes only names this run gives its working
files, directly inside the run folder, and never anything passed with `--protect`.

## Where things live

- `src/ltrquest/`: the Python stages. `scripts/ltrquest.sh` is the driver and
  holds the `--help` and `--help-all` text; `scripts/plots.sh` drives the plots.
- `src/ltrquest/fastltr/`: LTRquest's adaptation of FastLTR (`FastLTR.c`, with
  FASTAN's license beside it).
- `tests/`: pytest. Small synthetic fixtures; `tests/data/` holds the shipped
  *A. thaliana* chr2 and protein set used by the end-to-end test.
  The non-slow suite needs none of the external tools and finishes in seconds. Tests marked `slow` run the full toolchain and never run in CI.
- `dev/check.sh`: runs what CI runs (below).
- `main.nf`, `workflows/`, `modules/`, `subworkflows/`, `conf/`, `nextflow.config`,
  `nextflow_schema.json`: the Nextflow pipeline.
- `recipe/meta.yaml`: the conda recipe. `Dockerfile`, `bin/entrypoint.sh`,
  `bin/ltrquest-container`: the image and its launcher.
- `benchmarks/`: the benchmark pages and the code to rerun them.
- `docs/images/`: the README figures. `workflow_alt.png` (one genome, step by step) is a spare version of the
  workflow figure, not shown in the README.

## A note on Kmer2LTR boundaries

Kmer2LTR re-classifies LTR regions with a more sensitive approach. It usually
adjusts only a few bases, but for LTR-RTs with internal repeats it can adjust
them substantially, and the internal region then becomes quite small. Ideas not
yet taken: drop candidates whose internal region is shorter than 80 bp, or
whose boundaries Kmer2LTR moved by more than 20%; or fall back to the LTRharvest
and LTR_FINDER boundaries when Kmer2LTR moves them a lot.

## Before every push: run what CI runs

Each push to `main` runs the jobs in `.github/workflows/ci.yml` and
`docker.yml`, and one failure marks the commit "Some checks were not
successful". A job that is already red hides the next failure in it: from
24 Sep 2026 every push failed ruff, so shellcheck, the step after it, never ran;
the newest-Nextflow stub run failed from 12 Oct, the day Nextflow 26.10 shipped;
and the image build failed from 26 Oct, so no image was published after that.

```bash
bash dev/check.sh
```

runs ruff, shellcheck, the non-slow tests and the Nextflow stub runs, and exits
non-zero if one fails. A check whose tool is missing is listed as skipped, not
passed. Point it at both Nextflow versions CI uses:

```bash
mamba create -n ltrquest_nf -c conda-forge -c bioconda nextflow shellcheck  # newest
mamba create -n ltrquest_nf2404 -c conda-forge -c bioconda nextflow=24.04.4  # oldest allowed
NXF_OLD=<envs>/ltrquest_nf2404/bin/nextflow NXF_NEW=<envs>/ltrquest_nf/bin/nextflow \
    PATH=$PATH:<envs>/ltrquest_nf/bin bash dev/check.sh
```

Run it from an environment with the package's dev extras (`pip install -e
".[dev]"`) and the tools in `environment.yml`, or the tests that need blastn or
Kmer2LTR skip.

### The jobs, and what has turned each red

- **lint: ruff** (`ruff check .`; `ruff check --fix .` sorts imports). Covers
  `tests/` too, at a Python 3.9 floor; CI installs the newest ruff.
- **lint: shellcheck** (`--severity=warning` on `src/ltrquest/scripts/*.sh`,
  `bin/entrypoint.sh`, `bin/ltrquest-container`). Runs only once ruff passes.
- **pytest** on Python 3.9 (on `ubuntu-24.04`; `ubuntu-latest` has no 3.9 builds),
  3.11 and 3.13 (`pytest -m "not slow"`). Locally,
  put `src/` first on `PYTHONPATH`, as `dev/check.sh` does: in an environment
  with an older LTRquest installed, pytest otherwise tests that copy.
- **nextflow stub** under 24.04.4 and `latest-stable`: the test profile, then
  `--max_rounds 1`, then the off switches (`--reboundary false --strand_recovery
  off`, which must leave both steps out of `trace.txt`; `dev/check.sh` skips this
  third run). `latest-stable` picks up each
  Nextflow release the day it ships, and releases keep tightening the syntax:
  26.10 refused `return` in a workflow body (`main.nf` now uses `exit 0`).
  `nextflow lint .` under the newest Nextflow lists the errors; its warnings
  (`Channel` for `channel`, `when:` blocks, implicit `it`) are the next ones.
- **pip install** builds the wheel and runs each entry point's `--help`. Not in
  `dev/check.sh`; packaging changes need it.
- **container: build and smoke-test**. Not in `dev/check.sh`; needs Docker
  (`docker build .`). Dockerfile `RUN` lines execute under `bash -o pipefail`
  (micromamba's `SHELL`), so `tool | grep -q x` fails when the tool exits
  non-zero even though grep matched: write `{ tool 2>&1 || true; } | grep -q x`.
  Every `ARG *_REF` must be a commit that exists on GitHub, not only locally.
- **Slow tests** (`pytest -m slow`) never run in CI. Run them after touching the
  driver or a stage's outputs; `LTRQUEST_TOOLS_DIR` set to a built tools
  directory skips the clones.

### After pushing

Check the Actions tab (or `gh run list --limit 4`) and fix a red check before
the next push.

## Release checklist

1. Move the `[Unreleased]` entries in [CHANGELOG.md](../CHANGELOG.md) under the new
   version and date.
2. Set the version to `X.Y.Z` in each of these (they must agree):
   `pyproject.toml`, `src/ltrquest/__init__.py`, `CITATION.cff`,
   `recipe/meta.yaml`, the `org.opencontainers.image.version` label in
   `Dockerfile`, and `manifest.version` in `nextflow.config`. In the Nextflow
   modules (`modules/local/ltrquest/*/main.nf`), also set the
   `ghcr.io/cwb14/ltrquest:X.Y.Z` container tags
   (`test_every_module_pins_the_current_image` checks them) and the
   `ltrquest: X.Y.Z` lines in their `stub:` blocks;
   `grep -rnE 'ltrquest: ?[0-9]' modules` lists both.
3. Update the image tag in the README and [install.md](install.md) if they name one.
   Other places that name the version in examples: the `Dockerfile` comment
   lines 11-12, the docstring at `tests/test_end_to_end.py` line 10 (and
   `bin/ltrquest-container` if any version remains there).
4. Run `bash dev/check.sh`, and the slow tests after any change to the driver.
5. Commit, push, and wait for CI and the container build to go green.
6. Tag `vX.Y.Z` and push the tag; the conda recipe downloads the tag's tarball.
