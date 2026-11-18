"""The Nextflow pipeline follows the command line's defaults.

Capped callers, sdust and re-boundarying are on unless the user turns them off;
strand recovery is off unless the user picks a preset.
"""

import json
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG = REPO_ROOT / "nextflow.config"

pytestmark = pytest.mark.skipif(not CONFIG.is_file(), reason="not running from a source checkout")


def _param(name):
    m = re.search(rf"^\s*{name}\s*=\s*(.+?)\s*(?://.*)?$", CONFIG.read_text(), re.MULTILINE)
    assert m, name
    return m.group(1)


def _schema_props():
    schema = json.loads((REPO_ROOT / "nextflow_schema.json").read_text())
    return {k: v for d in schema["$defs"].values() for k, v in d.get("properties", {}).items()}


def test_config_declares_the_cli_defaults():
    assert _param("ltr_timeout") == "300"
    assert _param("skip_sdust") == "false"
    assert _param("strand_recovery") == "'off'"
    assert _param("strand_recovery_ppt") == "true"
    assert _param("reboundary") == "true"


def test_schema_lists_the_same_defaults():
    props = _schema_props()
    assert props["ltr_timeout"]["default"] == 300
    assert props["skip_sdust"] == {**props["skip_sdust"], "type": "boolean", "default": False}
    assert props["strand_recovery"]["default"] == "off"
    assert "off" in props["strand_recovery"]["enum"]
    assert props["strand_recovery_ppt"]["default"] is True
    assert props["reboundary"]["default"] is True


def test_detect_passes_sdust_unless_skipped():
    text = (REPO_ROOT / "modules/local/ltrquest/detect/main.nf").read_text()
    script = text[text.index("script:"):text.index("stub:")]
    assert "params.skip_sdust" in script and "'--sdust'" in script
    assert "${sdust_opt}" in script
    # detect.py's own defaults already equal the driver's, so no value is passed
    assert "--sdust-args" not in script and "--max-dust-frac" not in script
    env = (REPO_ROOT / "modules/local/ltrquest/detect/environment.yml").read_text()
    assert re.search(r"^\s*-\s*sdust\b", env, re.MULTILINE)


def test_off_turns_strand_recovery_off_everywhere():
    wf = (REPO_ROOT / "workflows/ltrquest.nf").read_text()
    assert "strand_recovery.toString() in ['off', 'false']" in wf
    # Nextflow 26.10 passes `--reboundary false` as the String "false" (truthy)
    assert "params.reboundary.toString().toBoolean()" in wf
    assert "if (params.reboundary)" not in wf
    assert "if (params.strand_recovery)" not in wf
    assert "strand_recovery_ppt && !" not in wf          # default-on ppt must not trip the check
    mc = (REPO_ROOT / "conf/modules.config").read_text()
    assert "!params.strand_recovery &&" not in mc
    assert "'off'" in mc


def test_help_text_shows_the_new_defaults_and_switches():
    main = (REPO_ROOT / "main.nf").read_text()
    for flag in ("--skip_sdust", "--reboundary", "--strand_recovery_ppt"):
        assert flag in main
    assert "off by default" not in main


def test_docs_say_how_to_turn_each_stage_off():
    doc = (REPO_ROOT / "docs/nextflow.md").read_text()
    for flag in ("--skip_sdust", "--strand_recovery off", "--strand_recovery_ppt false",
                 "--reboundary false", "--ltr_timeout 0"):
        assert flag in doc
