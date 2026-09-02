"""Shared pytest configuration and intentionally small scientific fixtures."""

from __future__ import annotations

import os
from pathlib import Path
import pickle
import subprocess
import sys
import uuid

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))
_pythonpath = os.environ.get("PYTHONPATH")
os.environ["PYTHONPATH"] = (
    str(SRC_ROOT) if not _pythonpath else os.pathsep.join((str(SRC_ROOT), _pythonpath))
)

# Keep imports and plotting deterministic in restricted local environments.
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/srfm-matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/srfm-numba")


@pytest.fixture
def tiny_ri_file(tmp_path: Path) -> Path:
    """Three-point refractive-index spectrum with an analytically linear trend."""
    path = tmp_path / "synthetic.ri"
    path.write_text(
        "# FORMAT = wavl n k\n"
        "# DESCRIPTION = pytest synthetic refractive indices\n"
        "1.0 1.40 0.010\n"
        "2.0 1.50 0.020\n"
        "4.0 1.70 0.040\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def descending_ri_file(tmp_path: Path) -> Path:
    path = tmp_path / "descending.ri"
    path.write_text(
        "# FORMAT = wavn n k\n"
        "10000 1.40 0.010\n"
        "5000 1.50 0.020\n"
        "2500 1.70 0.040\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def tiny_iasi_ils(tmp_path: Path) -> Path:
    """Small symmetric RFM-format instrument line shape for E2E convolution."""
    path = tmp_path / "tiny.ils"
    path.write_text(
        "! Synthetic pytest IASI-like ILS\n"
        "! No external instrument files required\n"
        "! Npts, first offset, spacing\n"
        "5 -0.250000 0.125000\n"
        "0.25 0.5 1.0 0.5 0.25\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def tiny_atmosphere(tmp_path: Path) -> Path:
    """Three-level atmosphere suitable for parser and native RFM tests."""
    from srfm.rfm_functions import write_atm_file

    path = tmp_path / "tiny.atm"
    write_atm_file(
        {
            "HGT [km]": np.array([0.0, 1.0, 2.0]),
            "PRE [mb]": np.array([1013.0, 900.0, 800.0]),
            "TEM [K]": np.array([290.0, 284.0, 278.0]),
            "F11 [ppmv]": np.array([1.0, 1.0, 1.0]),
        },
        path,
        header="Synthetic pytest atmosphere",
    )
    return path


@pytest.fixture
def tiny_altitude_grid(tmp_path: Path) -> Path:
    from srfm.rfm_functions import write_atm_file

    path = tmp_path / "heights.atm"
    write_atm_file({"HGT [km]": [0.0, 1.0, 2.0]}, path)
    return path


@pytest.fixture
def tiny_xsc_file(tmp_path: Path) -> Path:
    """Tiny F11 HITRAN-format cross section, requiring no external database."""
    from srfm.rfm_functions import write_xsc_file

    path = tmp_path / "f11.xsc"
    write_xsc_file(
        {
            "surface": {
                "molec": "F11",
                "low_spc": 999.0,
                "upp_spc": 1001.0,
                "npts": 3,
                "temperature": 290.0,
                "pressure": 760.0,
                "beta_ext": [1.0e-20, 2.0e-20, 1.0e-20],
            },
            "aloft": {
                "molec": "F11",
                "low_spc": 999.0,
                "upp_spc": 1001.0,
                "npts": 3,
                "temperature": 278.0,
                "pressure": 600.0,
                "beta_ext": [0.8e-20, 1.6e-20, 0.8e-20],
            },
        },
        path,
    )
    return path


@pytest.fixture
def require_native():
    """Return a helper that skips with an actionable local build command."""

    def _require(module, name: str) -> None:
        if module is None:
            pytest.skip(
                f"{name} extension is unavailable; run python build_extensions.py"
            )

    return _require


@pytest.fixture
def run_native_case(tmp_path: Path):
    """Run one native test case outside the pytest process.

    A missing result file is treated as a failure even when a Fortran ``STOP``
    exits with status zero.
    """

    def _run(case: str, payload: dict):
        token = uuid.uuid4().hex
        args_path = tmp_path / f"native-{token}-args.pkl"
        result_path = tmp_path / f"native-{token}-result.pkl"
        with args_path.open("wb") as handle:
            pickle.dump(payload, handle)

        completed = subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / "tests" / "_native_cases.py"),
                case,
                str(args_path),
                str(result_path),
            ],
            cwd=REPO_ROOT,
            env=os.environ.copy(),
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0:
            pytest.fail(
                f"Native case {case!r} exited with status {completed.returncode}.\n"
                f"stdout:\n{completed.stdout}\n"
                f"stderr:\n{completed.stderr}",
                pytrace=False,
            )
        if not result_path.is_file():
            pytest.fail(
                f"Native case {case!r} exited without a result payload; a Fortran "
                f"STOP or native early termination is likely.\n"
                f"stdout:\n{completed.stdout}\n"
                f"stderr:\n{completed.stderr}",
                pytrace=False,
            )
        with result_path.open("rb") as handle:
            result = pickle.load(handle)
        return result, completed

    return _run
