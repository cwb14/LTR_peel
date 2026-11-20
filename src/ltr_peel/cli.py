"""Console entry point for the ``ltr_peel`` pipeline driver.

The driver itself is ``scripts/ltr_peel.sh``: the round loop is a long-running
sequence of external tools, and bash is the right language for that. This module
exists to make it a first-class installed command, and to hand it two things it
cannot work out for itself:

* ``LTR_PEEL_PYTHON`` - the interpreter that owns this package. Without it the
  driver would fall back to ``python3`` on ``PATH``, which in a conda/venv layout
  is not necessarily the one ``ltr_peel`` was installed into.
* the absolute path of the packaged script, wherever pip put it.

``os.execv`` replaces this process rather than spawning a child, so signals,
exit codes and job control behave exactly as if bash had been invoked directly.
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parent / "scripts"


def _script(name: str) -> str:
    path = _SCRIPTS / name
    if not path.is_file():
        sys.exit(
            f"ltr_peel: packaged script is missing: {path}\n"
            "The installation is incomplete; try reinstalling the package."
        )
    return str(path)


def _exec(name: str, argv: list[str]) -> "None":
    bash = shutil.which("bash")
    if bash is None:
        sys.exit("ltr_peel: bash is required but was not found on PATH.")

    env = dict(os.environ)
    env.setdefault("LTR_PEEL_PYTHON", sys.executable)
    env["LTR_PEEL_ARGV0"] = Path(sys.argv[0]).name or name
    os.execve(bash, [bash, _script(name), *argv], env)


def main() -> "None":
    """Run the full pipeline: ``ltr_peel --genome ... [options]``."""
    _exec("ltr_peel.sh", sys.argv[1:])


def plots() -> "None":
    """Run the plotting stage alone: ``ltr_peel-plots --prefix ... --genome ...``."""
    _exec("plots.sh", sys.argv[1:])


if __name__ == "__main__":  # pragma: no cover - exercised via the console script
    main()
