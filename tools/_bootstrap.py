"""
Helpers shared by the custom build/install scripts.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]


def _run(cmd: Sequence[str]) -> None:
    display = " ".join(cmd)
    print(f"[srfm] {display}")
    subprocess.run(cmd, check=True)


def ensure_python_packages(packages: Iterable[tuple[str, str]]) -> None:
    """
    Ensure that the supplied Python modules are importable, installing the
    matching pip package name when needed.
    """
    missing: set[str] = set()
    for module_name, pip_name in packages:
        try:
            __import__(module_name)
        except ModuleNotFoundError:
            missing.add(pip_name)
    if missing:
        _run([sys.executable, "-m", "pip", "install", "--upgrade", *sorted(missing)])


def ensure_meson_and_ninja() -> None:
    """
    Meson and Ninja are required when numpy.f2py uses the Meson backend.
    Install both if either the module or executable is missing.
    """
    requirements: list[tuple[str, str]] = [("mesonbuild", "meson")]
    ninja_in_path = shutil.which("ninja") is not None
    if not ninja_in_path:
        requirements.append(("ninja", "ninja"))
    ensure_python_packages(requirements)


def reraise_from_process(callable_args: Sequence[str]) -> None:
    """
    Run a subprocess while mirroring the command prior to execution. This is
    primarily a convenience wrapper to keep the higher level scripts tidy.
    """
    _run(list(callable_args))

