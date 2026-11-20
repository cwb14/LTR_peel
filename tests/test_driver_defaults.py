"""The driver's defaults (sdust, strand recovery, PPT, re-boundarying, LTR timeout) and their off switches."""

from __future__ import annotations

import subprocess
import sys

import pytest


def run(*args, **kw):
    return subprocess.run(args, capture_output=True, text=True, **kw)


def test_defaults_block(driver):
    text = driver.read_text()
    for line in ("RUN_SDUST=true", 'STRAND_RECOVERY=""', "STRAND_RECOVERY_PPT=true",
                 "REBOUNDARY=true", "LTR_TIMEOUT=300"):
        assert any(raw.startswith(line) for raw in text.splitlines()), line


@pytest.mark.parametrize("flags", [
    ["--no-sdust"], ["--no-strand-recovery-ppt"], ["--no-reboundary"],
    ["--strand-recovery", "off"], ["--ltr-timeout", "0"],
])
def test_off_switches_parse(driver, flags):
    assert run("bash", str(driver), *flags, "--help").returncode == 0


def test_legacy_switches_parse(driver):
    result = run("bash", str(driver), "--run-sdust", "--reboundary", "--strand-recovery-ppt",
                 "--run-trf", "--help")
    assert result.returncode == 0


@pytest.mark.parametrize("value", ["-5", "abc"])
def test_bad_timeout_rejected(driver, toy_genome, tmp_path, value):
    result = run("bash", str(driver), "--genome", str(toy_genome), "--ltr-timeout", value,
                 cwd=tmp_path)
    assert result.returncode != 0 and "ltr-timeout" in result.stderr


def test_missing_blast_names_the_off_switches(driver, toy_genome, tmp_path):
    result = run("bash", str(driver), "--genome", str(toy_genome),
                 env={"PATH": "/usr/bin:/bin", "LTR_PEEL_PYTHON": "/nonexistent/python"},
                 cwd=tmp_path)
    assert result.returncode != 0
    assert "--strand-recovery off" in result.stderr and "--no-reboundary" in result.stderr


def test_detect_gets_the_timeout_before_extra_args(driver):
    text = driver.read_text()
    call = text[text.index('"${DETECT[@]}" \\'):]
    assert call.index('--ltr-timeout "$LTR_TIMEOUT"') < call.index('"${extra_round_args[@]}"')


def test_detect_default_is_300(driver):
    """ltr_peel-detect's own default matches the driver's, so every route caps chunks at 300 s."""
    help_text = subprocess.run([sys.executable, "-m", "ltr_peel.detect", "--help"],
                               capture_output=True, text=True).stdout
    assert "default 300" in " ".join(help_text.split())


def test_timeout_is_a_detection_setting(driver):
    text = driver.read_text()
    block = text[text.index("detection_settings() {"):text.index("plan_reuse() {")]
    assert 'ltr_timeout=${LTR_TIMEOUT}' in block


def test_a_relative_tools_dir_from_the_environment_is_made_absolute(driver, toy_genome, tmp_path):
    """Workers run inside round<N>/ and the merged stage in its own folder: a relative
    LTR_PEEL_TOOLS_DIR would put the tools there, and promoting the results removes a
    folder of that name from the run folder."""
    env = {"PATH": "/usr/bin:/bin", "LTR_PEEL_PYTHON": "/nonexistent/python",
           "LTR_PEEL_TOOLS_DIR": "./rel_tools"}
    result = run("bash", str(driver), "--genome", str(toy_genome), "--strand-recovery", "off",
                 "--no-reboundary", env=env, cwd=tmp_path)
    assert f"Tools directory: {tmp_path / 'rel_tools'}\n" in result.stdout


def test_the_orchestrator_exports_the_absolute_environment_tools_dir(driver):
    text = driver.read_text()
    block = text[text.index('elif [[ -n "${LTR_PEEL_TOOLS_DIR:-}" && -z "${_LTR_PEEL_WORKER:-}" ]]'):]
    block = block[:block.index("\nfi\n")]
    assert '"${PWD}/${LTR_PEEL_TOOLS_DIR#./}"' in block and "export LTR_PEEL_TOOLS_DIR" in block
