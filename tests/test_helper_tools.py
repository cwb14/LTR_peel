"""Resolution of the helper binaries LTRquest shells out to.

`ensure_tools` used to clone and compile minimap2 and miniprot from source
before it would consider the copies already on PATH. That made the container
unusable -- it ships no compiler and no git, on purpose -- and made every
conda user pay a source build for tools their environment already had. These
tests pin the order so it cannot regress.
"""

from __future__ import annotations

import stat
import subprocess
from pathlib import Path

import pytest

from ltrquest import detect, fastltr, helper_tools, kmer2ltr
from ltrquest.detect import _resolve_helper, ensure_sdust, ensure_tebinsorter, ensure_trfmod

URL = "https://example.invalid/never-cloned"
OFFLINE = ("Offline? Run once with internet using --tools-dir DIR, then pass the same "
           "--tools-dir; or use the container.")


def fake_binary(path: Path) -> Path:
    """An executable that satisfies a `usable` predicate."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\nexit 0\n")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def is_present(path) -> bool:
    """A realistic `usable` predicate: a binary is usable if it is really there."""
    return Path(path).exists()


def never_usable(_p) -> bool:
    return False


class TestResolutionOrder:
    def test_a_prebuilt_copy_in_tools_dir_wins(self, tmp_path, monkeypatch):
        # This is the container's case: /opt/ltrquest/tools is populated at
        # build time and must be used as-is.
        prebuilt = fake_binary(tmp_path / "minimap2" / "minimap2")
        monkeypatch.setattr("shutil.which", lambda n: "/usr/bin/minimap2")
        assert _resolve_helper("minimap2", tmp_path, URL, is_present) == prebuilt

    def test_falls_back_to_the_copy_on_path(self, tmp_path, monkeypatch):
        # The conda case: nothing prebuilt, but the environment has the tool.
        on_path = fake_binary(tmp_path / "bin" / "miniprot")
        monkeypatch.setattr("shutil.which", lambda n: str(on_path) if n == "miniprot" else None)
        assert _resolve_helper("miniprot", tmp_path, URL, is_present) == on_path

    def test_path_is_consulted_before_any_clone(self, tmp_path, monkeypatch):
        # The regression itself: a usable binary on PATH must short-circuit,
        # so neither git nor a compiler is ever needed.
        on_path = fake_binary(tmp_path / "bin" / "minimap2")
        monkeypatch.setattr("shutil.which", lambda n: str(on_path) if n == "minimap2" else None)

        def explode(*a, **k):
            raise AssertionError("cloned or built despite a usable binary on PATH")

        monkeypatch.setattr("ltrquest.detect.run", explode)
        assert _resolve_helper("minimap2", tmp_path, URL, is_present) == on_path

    def test_an_unusable_binary_on_path_is_not_accepted(self, tmp_path, monkeypatch):
        fake_binary(tmp_path / "bin" / "minimap2")
        monkeypatch.setattr("shutil.which", lambda n: str(tmp_path / "bin" / "minimap2"))
        with pytest.raises(RuntimeError):
            _resolve_helper("minimap2", tmp_path, URL, never_usable)


class TestErrorsAreActionable:
    def test_missing_git_names_the_problem_and_the_fix(self, tmp_path, monkeypatch):
        # What the container hit. The old code raised a bare
        # FileNotFoundError('git') from six frames down inside subprocess.
        monkeypatch.setattr("shutil.which", lambda n: None)
        with pytest.raises(RuntimeError) as excinfo:
            _resolve_helper("minimap2", tmp_path, URL, never_usable)
        message = str(excinfo.value)
        assert "git is not installed" in message
        assert "conda install" in message
        assert "--tools-dir" in message

    def test_a_failed_build_names_the_tool(self, tmp_path, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda n: "/usr/bin/git" if n == "git" else None)
        monkeypatch.setattr("ltrquest.detect.run", lambda *a, **k: type("R", (), {"returncode": 1})())
        (tmp_path / "miniprot").mkdir()
        with pytest.raises(RuntimeError, match="miniprot build failed") as excinfo:
            _resolve_helper("miniprot", tmp_path, URL, never_usable)
        assert "Offline?" not in str(excinfo.value)      # a compiler problem, not a network one

    @pytest.mark.parametrize("git", [None, "/usr/bin/git"], ids=["no-git", "failed-build"])
    def test_an_sdust_failure_names_its_off_switch(self, tmp_path, monkeypatch, git):
        # sdust runs by default, so the way around a failure is its off switch, as
        # TRF-mod's names --no-trf and FastLTR's --no-fastltr.
        monkeypatch.setattr("shutil.which", lambda n: git if n == "git" else None)
        monkeypatch.setattr("ltrquest.detect.run", lambda *a, **k: type("R", (), {"returncode": 1})())
        if git:
            (tmp_path / "sdust").mkdir()                 # cloned; the build then fails
        with pytest.raises(RuntimeError) as excinfo:
            ensure_sdust(tmp_path)
        assert "--no-sdust" in str(excinfo.value)

    @pytest.mark.parametrize("fetch", [
        lambda d: _resolve_helper("miniprot", d, URL, never_usable),
        ensure_tebinsorter,
        ensure_trfmod,
        ensure_sdust,
    ], ids=["minimap2-miniprot", "TEsorter2", "TRF-mod", "sdust"])
    def test_a_failed_clone_says_what_to_do_offline(self, tmp_path, monkeypatch, fetch):
        # A compute node with no internet: git is there, every clone fails.
        def clone_fails(cmd, cwd=None, check=False, **kw):
            result = type("R", (), {"returncode": 128, "stdout": "",
                                    "stderr": "fatal: unable to access"})()
            if check and cmd[:2] == ["git", "clone"]:
                raise RuntimeError(f"Command failed (128): {' '.join(cmd)}\n"
                                   f"STDERR:\n{result.stderr}")
            return result

        monkeypatch.setattr("shutil.which", lambda n: "/usr/bin/git" if n == "git" else None)
        monkeypatch.setattr("ltrquest.detect.run", clone_fails)
        monkeypatch.delenv("TESORTER2_SRC", raising=False)
        monkeypatch.delenv("TEBINSORTER_SRC", raising=False)
        with pytest.raises(RuntimeError) as excinfo:
            fetch(tmp_path)
        assert str(excinfo.value).endswith(OFFLINE)


def _no_internet(cmd, *args, **kwargs):
    """subprocess.run on a compute node with no internet: every git clone fails."""
    if list(cmd[:2]) == ["git", "clone"]:
        raise subprocess.CalledProcessError(128, cmd)
    return subprocess.CompletedProcess(cmd, 1, "", "")


class TestOfflineBeyondDetect:
    """FastLTR's FASTAN and Kmer2LTR are fetched outside detect's helpers."""

    def test_a_failed_fastan_clone_says_what_to_do_offline(self, tmp_path, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda n: "/usr/bin/git" if n == "git" else None)
        monkeypatch.setattr(fastltr.subprocess, "run", _no_internet)
        with pytest.raises(RuntimeError) as excinfo:
            fastltr.ensure(tmp_path)
        assert str(excinfo.value).endswith(OFFLINE)

    def test_a_failed_kmer2ltr_clone_says_what_to_do_offline(self, tmp_path, monkeypatch):
        monkeypatch.setattr(kmer2ltr.subprocess, "run", _no_internet)
        with pytest.raises(RuntimeError) as excinfo:
            kmer2ltr._clone(tmp_path)
        assert str(excinfo.value).endswith(OFFLINE)

    def test_every_fetcher_says_the_same_sentence(self):
        assert helper_tools.OFFLINE_HINT == OFFLINE == detect.OFFLINE_HINT
