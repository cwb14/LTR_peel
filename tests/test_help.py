"""The two help pages: a short one for a first run, the full one on request."""

from __future__ import annotations

import re
import subprocess

LEGACY = {"--run-sdust", "--reboundary", "--strand-recovery-ppt", "--run-trf"}


def run(driver, *args):
    return subprocess.run(["bash", str(driver), *args], capture_output=True, text=True)


def test_short_help_is_short(driver):
    result = run(driver, "--help")
    assert result.returncode == 0
    lines = result.stdout.splitlines()
    assert len(lines) <= 35
    assert max(map(len, lines)) <= 80
    for flag in ("--genome", "--proteins", "--threads", "--mutation-rate", "--out_prefix",
                 "--keep-intermediates", "--help-all"):
        assert flag in result.stdout


def test_help_is_only_the_help(driver):
    for args in (("--help",), ("-h",), ("--help-all",), ()):
        out = run(driver, *args)
        assert "Command:" not in out.stdout + out.stderr, args
    assert run(driver, "-h").stdout == run(driver, "--help").stdout


def test_no_args_prints_short_help(driver):
    bare, short = run(driver), run(driver, "--help")
    assert bare.returncode == 1
    assert bare.stdout + bare.stderr == short.stdout


def test_help_all_lists_every_option(driver):
    result = run(driver, "--help-all")
    assert result.returncode == 0
    flags = set(re.findall(r"^\s+(--[\w-]+)(?:\|[^)]*)?\)", driver.read_text(), re.M))
    flags = {f for f in flags if not f.startswith("--dev-")} - LEGACY - {"--help", "--help-all"}
    assert flags
    for flag in sorted(flags):
        assert flag in result.stdout, flag
    assert max(map(len, result.stdout.splitlines())) <= 80


def test_help_all_states_new_defaults(driver):
    text = run(driver, "--help-all").stdout
    for word in ("300", "sensitive", "--no-sdust", "--no-reboundary",
                 "--no-strand-recovery-ppt", "--strand-recovery MODE", "off"):
        assert word in text, word


def test_help_all_explains_what_keeping_and_the_tools_dir_do(driver):
    text = " ".join(run(driver, "--help-all").stdout.split())
    assert "re-plotting with ltr_peel-plots" in text
    assert "never removed" in text and "only if this run created it" in text
