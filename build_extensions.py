#!/usr/bin/env python3
"""Cross-platform helper for building the SRFM Fortran extensions."""

from __future__ import annotations

import argparse
import datetime as _dt
import heapq
import os
import re
import shlex
import shutil
import subprocess
import sys
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set

try:
    from tools import _bootstrap as _build_bootstrap
except ModuleNotFoundError:
    _build_bootstrap = None


REPO_ROOT = Path(__file__).resolve().parent
SRC_ROOT = REPO_ROOT / "src" / "srfm"
LOG_FILE = SRC_ROOT / "build.log"
DEFAULT_FFLAGS = "-O3 -g -Wall"
F2PY_CMD = [sys.executable, "-m", "numpy.f2py"]
_CUSTOM_FFLAGS = False
_MESON_FLAG_WARNING = False
MODULE_DEF_RE = re.compile(r"^\s*module\s+(?!procedure\b)(?!function\b)(\w+)", re.I)
USE_RE = re.compile(r"^\s*use\s+(?:(?:,\s*intrinsic\s*)?::\s*)?(\w+)", re.I)
INTRINSIC_MODULES = {
    "iso_c_binding",
    "iso_fortran_env",
    "omp_lib",
    "omp_lib_kinds",
}


def _ensure_meson_and_ninja() -> None:
    missing = []
    try:
        import mesonbuild  # type: ignore
    except ModuleNotFoundError:
        missing.append("meson")
    if shutil.which("ninja") is None:
        missing.append("ninja")
    if missing:
        print("[build] Installing required build tools:", ", ".join(missing))
        subprocess.run(
            [sys.executable, "-m", "pip", "install", "--upgrade", *missing],
            check=True,
        )


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the SRFM Fortran extensions with numpy.f2py."
    )
    parser.add_argument(
        "components",
        nargs="*",
        choices=["all", "mie", "disort_single", "disort_double", "rfm"],
        help=(
            "Subset of components to build. "
            "If omitted or if 'all' is supplied, every component is compiled."
        ),
    )
    parser.add_argument(
        "--fflags",
        default=None,
        help=(
            "Custom Fortran compiler flags. "
            "Defaults to the value of SRFM_FFLAGS or '-O3 -g -Wall'."
        ),
    )
    parser.add_argument(
        "--backend",
        choices=["auto", "meson", "numpy.distutils"],
        default="auto",
        help=(
            "f2py backend to use. "
            "'auto' selects Meson when it is available, "
            "and falls back to the default numpy.distutils backend otherwise."
        ),
    )
    return parser.parse_args(argv)


def _init_log(fflags: str, backend: str) -> None:
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    LOG_FILE.write_text(
        "SRFM native build log\n"
        f"Timestamp: {_dt.datetime.now().isoformat()}\n"
        f"Python: {sys.version.split()[0]}\n"
        f"fflags: {fflags}\n"
        f"backend: {backend}\n"
    )


def _log_command(cmd: Iterable[str]) -> None:
    human = " ".join(shlex.quote(str(part)) for part in cmd)
    with LOG_FILE.open("a", encoding="utf-8") as log:
        log.write(f"\n$ {human}\n")


def _run_command(
    cmd: List[str], cwd: Path | None = None, env: dict[str, str] | None = None
) -> None:
    str_cmd = [str(part) for part in cmd]
    human = " ".join(shlex.quote(part) for part in str_cmd)
    location = cwd if cwd is not None else Path.cwd()
    print(f"[build] {human} (cwd={location})")
    _log_command(str_cmd)
    with LOG_FILE.open("a", encoding="utf-8") as log:
        subprocess.run(
            str_cmd,
            cwd=cwd,
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
        )


def _resolve_backend(user_choice: str) -> list[str]:
    if user_choice == "auto":
        backend = "meson" if shutil.which("meson") else None
    else:
        backend = user_choice
    if backend is None:
        return []
    return ["--backend", backend]


def _ensure_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required source file '{path}' was not found.")


def _combine_sources(target_dir: Path, file_names: Iterable[str]) -> Path:
    combined = target_dir / "code.f"
    with combined.open("wb") as destination:
        for name in file_names:
            source = target_dir / name
            _ensure_exists(source)
            with source.open("rb") as src_file:
                shutil.copyfileobj(src_file, destination)
            destination.write(b"\n")
    return combined


def _combine_absolute(paths: Iterable[Path], output: Path) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as destination:
        for path in paths:
            _ensure_exists(path)
            with path.open("rb") as src_file:
                shutil.copyfileobj(src_file, destination)
            destination.write(b"\n")
    return output


def _parse_fortran_file(path: Path) -> tuple[Set[str], Set[str]]:
    modules: Set[str] = set()
    uses: Set[str] = set()
    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        # Strip comments
        line = raw_line.split("!", 1)[0].strip()
        if not line:
            continue
        module_match = MODULE_DEF_RE.match(line)
        if module_match:
            modules.add(module_match.group(1).lower())
            continue
        use_match = USE_RE.match(line)
        if use_match:
            name = use_match.group(1).lower()
            if name not in INTRINSIC_MODULES:
                uses.add(name)
    return modules, uses


def _order_sources_with_dependencies(source_paths: list[Path]) -> list[Path]:
    if not source_paths:
        return []

    file_modules: Dict[Path, Set[str]] = {}
    file_uses: Dict[Path, Set[str]] = {}
    module_map: Dict[str, Path] = {}

    for path in source_paths:
        modules, uses = _parse_fortran_file(path)
        file_modules[path] = modules
        file_uses[path] = uses
        for module_name in modules:
            module_map.setdefault(module_name, path)

    dependencies: Dict[Path, Set[Path]] = {path: set() for path in source_paths}
    reverse_deps: Dict[Path, Set[Path]] = defaultdict(set)

    for path in source_paths:
        for module_name in file_uses[path]:
            provider = module_map.get(module_name)
            if provider and provider != path:
                dependencies[path].add(provider)
                reverse_deps[provider].add(path)

    in_degree: Dict[Path, int] = {
        path: len(deps) for path, deps in dependencies.items()
    }
    queue: List[tuple[str, Path]] = [
        (path.name.lower(), path) for path, degree in in_degree.items() if degree == 0
    ]
    heapq.heapify(queue)
    ordered: List[Path] = []

    while queue:
        _, path = heapq.heappop(queue)
        ordered.append(path)
        for dependent in reverse_deps.get(path, ()):
            in_degree[dependent] -= 1
            if in_degree[dependent] == 0:
                heapq.heappush(queue, (dependent.name.lower(), dependent))

    if len(ordered) != len(source_paths):
        # Cycle detected; append remaining files in alphabetical order to avoid stalls.
        remaining = [path for path in source_paths if path not in ordered]
        ordered.extend(sorted(remaining, key=lambda item: item.name.lower()))

    return ordered


def _prepare_command(
    base_cmd: list[str], fflags: str, backend_args: list[str]
) -> tuple[list[str], dict[str, str] | None]:
    global _MESON_FLAG_WARNING
    cmd = list(base_cmd)
    env = None
    uses_meson = backend_args and backend_args[-1] == "meson"
    if not uses_meson:
        cmd.extend(["--f90flags", fflags])
    elif _CUSTOM_FFLAGS:
        if not _MESON_FLAG_WARNING:
            warnings.warn(
                "Custom Fortran flags are ignored when using the Meson backend "
                "(Python >= 3.12). Set compiler flags through your compiler "
                "environment variables if needed.",
                stacklevel=2,
            )
            _MESON_FLAG_WARNING = True
    cmd.extend(backend_args)
    return cmd, env


def _build_mie(fflags: str, backend_args: list[str]) -> None:
    print("Building mie_module...")
    source = SRC_ROOT / "mie_ewp.f90"
    _ensure_exists(source)
    base_cmd = F2PY_CMD + [
        "-c",
        str(source),
        "-m",
        "mie_module",
    ]
    cmd, env = _prepare_command(base_cmd, fflags, backend_args)
    _run_command(cmd, cwd=SRC_ROOT, env=env)


def _build_disort_single(fflags: str, backend_args: list[str]) -> None:
    print("Building disort_module_s (single precision)...")
    source_dir = SRC_ROOT / "DISORT"
    files = [
        "DISOTESTAUX.f",
        "DISORT.f",
        "BDREF.f",
        "DISOBRDF.f",
        "ERRPACK.f",
        "LINPAK.f",
        "LAPACK.f",
        "RDI1MACH.f",
    ]
    combined = _combine_sources(source_dir, files)
    base_cmd = F2PY_CMD + [
        "-c",
        str(combined),
        "-m",
        "disort_module_s",
    ]
    cmd, env = _prepare_command(base_cmd, fflags, backend_args)
    _run_command(cmd, cwd=source_dir, env=env)


def _build_disort_double(fflags: str, backend_args: list[str]) -> None:
    print("Building disort_module_d (double precision)...")
    source_dir = SRC_ROOT / "DISORT_dbl"
    files = [
        "DISOTESTAUX.f",
        "DISORT.f",
        "BDREF.f",
        "DISOBRDF.f",
        "ERRPACK.f",
        "LINPACK_D.f",
        "LAPACK.f",
        "RDI1MACH.f",
    ]
    combined = _combine_sources(source_dir, files)
    base_cmd = F2PY_CMD + [
        "-c",
        str(combined),
        "-m",
        "disort_module_d",
    ]
    cmd, env = _prepare_command(base_cmd, fflags, backend_args)
    _run_command(cmd, cwd=source_dir, env=env)


def _build_rfm(fflags: str, backend_args: list[str]) -> None:
    print("Building rfm_py...")
    rfm_dir = SRC_ROOT / "RFM"
    wrapper = rfm_dir / "rfm_wrapper.pyf"
    _ensure_exists(wrapper)
    source_dir = rfm_dir / "source"
    if not source_dir.is_dir():
        raise FileNotFoundError(f"RFM source folder '{source_dir}' was not found.")
    exclusions = {"rfm.f90", "combined.f90"}
    sources = [
        path for path in source_dir.glob("*.f90") if path.name.lower() not in exclusions
    ]
    sources = _order_sources_with_dependencies(sources)
    if not sources:
        raise FileNotFoundError("No RFM Fortran sources were discovered.")
    combined = _combine_absolute(sources, source_dir / "combined.f90")
    base_cmd = F2PY_CMD + [
        "-c",
        "-m",
        "rfm_py",
        str(wrapper),
        str(combined),
    ]
    cmd, env = _prepare_command(base_cmd, fflags, backend_args)
    try:
        _run_command(cmd, cwd=rfm_dir, env=env)
    finally:
        if combined.exists():
            combined.unlink()


BUILDERS = {
    "mie": _build_mie,
    "disort_single": _build_disort_single,
    "disort_double": _build_disort_double,
    "rfm": _build_rfm,
}


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    selection = list(BUILDERS)
    if args.components and "all" not in args.components:
        selection = args.components
    if _build_bootstrap is not None:
        _build_bootstrap.ensure_meson_and_ninja()
    else:
        _ensure_meson_and_ninja()
    fflags_source = args.fflags or os.environ.get("SRFM_FFLAGS")
    global _CUSTOM_FFLAGS
    _CUSTOM_FFLAGS = fflags_source is not None
    fflags = fflags_source or DEFAULT_FFLAGS
    backend_args = _resolve_backend(args.backend)
    backend_label = backend_args[-1] if backend_args else "numpy.distutils"
    _init_log(fflags, backend_label)
    try:
        for component in selection:
            BUILDERS[component](fflags, backend_args)
    except subprocess.CalledProcessError as exc:
        print(
            f"[error] '{' '.join(exc.cmd)}' failed with exit code {exc.returncode}. "
            f"See '{LOG_FILE}' for details."
        )
        raise SystemExit(exc.returncode) from exc
    except FileNotFoundError as exc:
        print(f"[error] {exc}")
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
