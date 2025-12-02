#!/usr/bin/env python3
"""
Install SRFM from the pre-built artifacts stored under dist/.
"""

from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path
from typing import Sequence

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parent.parent))
    import tools._bootstrap as _bootstrap  # type: ignore
else:
    from . import _bootstrap


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Install SRFM from the local dist/ directory."
    )
    parser.add_argument(
        "--dist",
        type=Path,
        default=_bootstrap.REPO_ROOT / "dist",
        help="Path to the dist/ folder containing wheels and/or an sdist.",
    )
    parser.add_argument(
        "--package",
        default="SRFM",
        help="Package name to install. Defaults to SRFM.",
    )
    parser.add_argument(
        "--pip-args",
        default="",
        help="Quoted string of additional arguments passed to pip before the package name.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    dist = args.dist
    if not dist.is_dir():
        raise SystemExit(f"dist folder '{dist}' does not exist.")
    artifacts = list(dist.glob("*.whl")) + list(dist.glob("*.tar.gz"))
    if not artifacts:
        raise SystemExit(f"No wheels or source archives found under '{dist}'.")
    _bootstrap.ensure_meson_and_ninja()
    cmd = [
        sys.executable,
        "-m",
        "pip",
        "install",
        "--upgrade",
        "--no-index",
        "--find-links",
        str(dist),
    ]
    cmd.extend(shlex.split(args.pip_args))
    cmd.append(args.package)
    _bootstrap.reraise_from_process(cmd)


if __name__ == "__main__":
    main()
