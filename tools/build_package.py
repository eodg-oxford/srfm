#!/usr/bin/env python3
"""
Cross-platform entry point that mirrors the developer Makefile targets.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Sequence

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parent.parent))
    import tools._bootstrap as _bootstrap  # type: ignore
else:
    from . import _bootstrap


REPO_ROOT = _bootstrap.REPO_ROOT
DIST_DIR = REPO_ROOT / "dist"
BUILD_EXTENSIONS = REPO_ROOT / "build_extensions.py"
INSTALL_SCRIPT = REPO_ROOT / "tools" / "install_from_dist.py"


def _run(cmd: Sequence[str]) -> None:
    printable = " ".join(cmd)
    print(f"[srfm] {printable}")
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def build_native(components: Sequence[str]) -> None:
    args = [sys.executable, str(BUILD_EXTENSIONS)]
    if components:
        args.extend(components)
    _bootstrap.ensure_meson_and_ninja()
    _run(args)


def build_artifact(kind: str) -> None:
    _bootstrap.ensure_python_packages([("build", "build")])
    _run(
        [
            sys.executable,
            "-m",
            "build",
            f"--{kind}",
            "--no-isolation",
            "--skip-dependency-check",
        ]
    )


def clean_dist() -> None:
    if DIST_DIR.exists():
        print(f"[srfm] removing '{DIST_DIR}'")
        shutil.rmtree(DIST_DIR)


def install_from_dist(pip_args: str) -> None:
    args = [sys.executable, str(INSTALL_SCRIPT)]
    if pip_args:
        args.extend(["--pip-args", pip_args])
    _run(args)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="SRFM build helper.")
    parser.add_argument(
        "action",
        choices=["native", "wheel", "sdist", "dist", "install", "clean-dist"],
        help="Task to run. 'dist' builds both wheel and sdist.",
    )
    parser.add_argument(
        "--components",
        nargs="*",
        default=(),
        help="Subset of native components to build (passed to build_extensions.py).",
    )
    parser.add_argument(
        "--pip-args",
        default="",
        help="Quoted string of extra arguments for pip when running the install action.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if args.action == "native":
        build_native(args.components)
    elif args.action == "wheel":
        build_artifact("wheel")
    elif args.action == "sdist":
        build_artifact("sdist")
    elif args.action == "dist":
        clean_dist()
        build_artifact("wheel")
        build_artifact("sdist")
    elif args.action == "clean-dist":
        clean_dist()
    elif args.action == "install":
        install_from_dist(args.pip_args)
    else:
        raise SystemExit(f"Unknown action: {args.action}")


if __name__ == "__main__":
    main()
