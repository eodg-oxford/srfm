"""Run every family in DISORT's authoritative Fortran validation driver."""

from __future__ import annotations

from pathlib import Path
from hashlib import sha256
import json
import os
import platform
import re
import shutil
import subprocess
import tempfile

import pytest

pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).resolve().parents[2]
DISORT_ROOT = REPO_ROOT / "src" / "srfm" / "DISORT_dbl"


def _cached_disort_executable(cache_directory, compiler, sources, compiler_flags):
    """Build the native validation driver when its compilation inputs change.

    Args:
        cache_directory (Path): Persistent pytest cache directory.
        compiler (str): Fortran compiler executable.
        sources (sequence[Path]): Ordered source files passed to the compiler.
        compiler_flags (sequence[str]): Compiler options affecting the executable.

    Returns:
        Path: Cached executable matching the sources, compiler, flags, and platform.

    Notes:
        Compilation uses a private directory and publishes by atomic rename.
        Failed or concurrent builds cannot expose a partially written executable.
        Only the executable is cached; the numerical validation runs every time.
    """
    version = subprocess.run(
        [compiler, "--version"], text=True, capture_output=True, check=True, timeout=20
    )
    build_inputs = {
        "compiler": str(Path(compiler).resolve()),
        "compiler_version": version.stdout,
        "flags": list(compiler_flags),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "sources": [
            (str(source.resolve()), sha256(source.read_bytes()).hexdigest())
            for source in sources
        ],
    }
    build_key = sha256(json.dumps(build_inputs, sort_keys=True).encode()).hexdigest()
    cache_directory.mkdir(parents=True, exist_ok=True)
    suffix = ".exe" if os.name == "nt" else ""
    cached_executable = cache_directory / f"disotest-{build_key}{suffix}"
    if cached_executable.is_file() and os.access(cached_executable, os.X_OK):
        return cached_executable

    with tempfile.TemporaryDirectory(prefix="build-", dir=cache_directory) as directory:
        executable = Path(directory) / f"disotest{suffix}"
        compiled = subprocess.run(
            [compiler, *compiler_flags, *(str(source) for source in sources),
             "-o", str(executable)],
            cwd=directory,
            text=True,
            capture_output=True,
            check=False,
            timeout=120,
        )
        assert compiled.returncode == 0, compiled.stdout + compiled.stderr
        executable.replace(cached_executable)
    return cached_executable


def test_all_disort_driver_families_match_reference_results(tmp_path, pytestconfig):
    """Execute all 17 families and 57 cases, compiling only on a cache miss.

    The driver owns the published reference arrays and comparison tolerances.
    Its process exit status does not encode comparison failures, so the final
    47-check numerical summary is asserted in addition to successful execution.
    """
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran is required to compile the DISORT validation driver")

    sources = (
        "disotest.f90",
        "DISOTESTAUX.f",
        "DISORT.f",
        "BDREF.f",
        "DISOBRDF.f",
        "ERRPACK.f",
        "LINPACK_D.f",
        "LAPACK.f",
        "RDI1MACH.f",
    )
    executable = _cached_disort_executable(
        pytestconfig.cache.mkdir("srfm-disort"),
        compiler,
        tuple(DISORT_ROOT / source for source in sources),
        ("-O2", "-fcheck=all"),
    )

    completed = subprocess.run(
        [str(executable)],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )
    output = completed.stdout + completed.stderr
    assert completed.returncode == 0, output

    cases = re.findall(r"Test Case No\.\s*(\d+)[a-r]?", output)
    assert {int(case) for case in cases} == set(range(1, 18))
    assert len(cases) == 57
    assert re.search(r"Percent of unit tests passed:\s+100\.00%", output)
    assert re.search(r"Number of unit tests:\s+47", output)
    assert re.search(r"Number of unit tests passed:\s+47", output)
    assert re.search(r"Number of unit tests failed:\s+0", output)


@pytest.mark.parametrize("changed_input", ["source", "flags", "version", "compiler"])
def test_disort_build_cache_reuses_and_invalidates_executable(
    tmp_path, monkeypatch, changed_input
):
    """An unchanged build is reused; each changed compilation input rebuilds it."""
    source = tmp_path / "driver.f90"
    source.write_text("program original\nend\n")
    compiler = str(tmp_path / "compiler")
    compiler_flags = ("-O2", "-fcheck=all")
    compiler_version = "test compiler 1"
    compilation_commands = []

    def compile_stub(command, **kwargs):
        if command[1] == "--version":
            return subprocess.CompletedProcess(command, 0, compiler_version, "")
        compilation_commands.append(command)
        executable = Path(command[command.index("-o") + 1])
        executable.write_bytes(b"compiled driver")
        executable.chmod(0o755)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(subprocess, "run", compile_stub)
    cache_directory = tmp_path / "cache"
    first = _cached_disort_executable(cache_directory, compiler, (source,), compiler_flags)
    assert _cached_disort_executable(
        cache_directory, compiler, (source,), compiler_flags
    ) == first
    assert len(compilation_commands) == 1

    if changed_input == "source":
        previous_stat = source.stat()
        source.write_text("program modified\nend\n")
        os.utime(source, ns=(previous_stat.st_atime_ns, previous_stat.st_mtime_ns))
    elif changed_input == "flags":
        compiler_flags = ("-O0", "-fcheck=all")
    elif changed_input == "version":
        compiler_version = "test compiler 2"
    else:
        compiler = str(tmp_path / "other-compiler")

    rebuilt = _cached_disort_executable(cache_directory, compiler, (source,), compiler_flags)
    assert rebuilt != first
    assert rebuilt.is_file()
    assert len(compilation_commands) == 2


def test_disort_build_cache_does_not_publish_failed_compilation(tmp_path, monkeypatch):
    """A failed compiler cannot leave a cache entry that a later test would reuse."""
    source = tmp_path / "driver.f90"
    source.write_text("program driver\nend\n")
    compilation_attempts = []

    def compile_stub(command, **kwargs):
        if command[1] == "--version":
            return subprocess.CompletedProcess(command, 0, "test compiler 1", "")
        compilation_attempts.append(command)
        executable = Path(command[command.index("-o") + 1])
        executable.write_bytes(b"driver")
        executable.chmod(0o755)
        status = 1 if len(compilation_attempts) == 1 else 0
        return subprocess.CompletedProcess(command, status, "", "compiler failed" if status else "")

    monkeypatch.setattr(subprocess, "run", compile_stub)
    cache_directory = tmp_path / "cache"
    with pytest.raises(AssertionError, match="compiler failed"):
        _cached_disort_executable(cache_directory, "compiler", (source,), ("-O2",))
    assert not list(cache_directory.iterdir())
    executable = _cached_disort_executable(cache_directory, "compiler", (source,), ("-O2",))
    assert executable.is_file()
    assert len(compilation_attempts) == 2
