"""Run every family in DISORT's authoritative Fortran validation driver."""

from __future__ import annotations

from pathlib import Path
import re
import shutil
import subprocess

import pytest

pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).resolve().parents[2]
DISORT_ROOT = REPO_ROOT / "src" / "srfm" / "DISORT_dbl"


def test_all_disort_driver_families_match_reference_results(tmp_path):
    """Compile and execute all 17 families and 57 cases in ``disotest``.

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
    executable = tmp_path / "disotest"
    compile_result = subprocess.run(
        [
            compiler,
            "-O2",
            "-fcheck=all",
            *(str(DISORT_ROOT / source) for source in sources),
            "-o",
            str(executable),
        ],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )
    assert compile_result.returncode == 0, compile_result.stderr

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
