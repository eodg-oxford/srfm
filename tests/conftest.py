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
    """Create a three-point refractive-index spectrum.

    Its analytically linear trend makes interpolation results transparent.

    Args:
        tmp_path: Pytest temporary-path fixture.

    Returns:
        Path to the generated refractive-index file.
    """
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
    """Create a descending three-point refractive-index spectrum.

    The fixture exercises readers that must normalize wavenumber ordering.

    Args:
        tmp_path: Pytest temporary-path fixture.

    Returns:
        Path to the generated refractive-index file.
    """
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
    """Create a small symmetric RFM-format instrument line shape.

    The compact kernel exercises E2E convolution without external IASI data.

    Args:
        tmp_path: Pytest temporary-path fixture.

    Returns:
        Path to the generated line-shape file.
    """
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
    """Create a three-level synthetic atmosphere.

    Its F11 profile is suitable for parser and native RFM tests.

    Args:
        tmp_path: Pytest temporary-path fixture.

    Returns:
        Path to the generated atmosphere file.
    """
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
    """Create a three-level altitude-only RFM atmosphere file.

    The file supplies the explicit output grid used by native RFM tests.

    Args:
        tmp_path: Pytest temporary-path fixture.

    Returns:
        Path to the generated altitude-grid file.
    """
    from srfm.rfm_functions import write_atm_file

    path = tmp_path / "heights.atm"
    write_atm_file({"HGT [km]": [0.0, 1.0, 2.0]}, path)
    return path


@pytest.fixture
def tiny_xsc_file(tmp_path: Path) -> Path:
    """Create a tiny F11 HITRAN-format cross section.

    Two atmospheric states avoid any dependency on an external database.

    Args:
        tmp_path: Pytest temporary-path fixture.

    Returns:
        Path to the generated cross-section file.
    """
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
def tiny_iasi_atmosphere(tmp_path: Path) -> Path:
    """Create an atmosphere containing every profile read by ``iasi_main``.

    The additional trace-gas profiles let the IASI runner construct its merged
    ECMWF/RFM atmosphere without relying on external climatology files.

    Args:
        tmp_path: Pytest temporary-path fixture.

    Returns:
        Path to the generated IASI-compatible atmosphere file.
    """
    from srfm.rfm_functions import write_atm_file

    path = tmp_path / "tiny-iasi.atm"
    values = {
        "HGT [km]": np.array([0.0, 1.0, 2.0]),
        "PRE [mb]": np.array([1013.0, 900.0, 800.0]),
        "TEM [K]": np.array([290.0, 284.0, 278.0]),
        "F11 [ppmv]": np.ones(3),
        "O3 [ppmv]": np.array([0.03, 0.04, 0.05]),
        "CO [ppmv]": np.full(3, 0.1),
        "N2O [ppmv]": np.full(3, 0.3),
        "CO2 [ppmv]": np.full(3, 420.0),
        "CH4 [ppmv]": np.full(3, 1.8),
        "H2O [ppmv]": np.array([8000.0, 4000.0, 2000.0]),
    }
    write_atm_file(values, path, header="Synthetic pytest IASI atmosphere")
    return path


@pytest.fixture
def tiny_iasi_observation(tmp_path: Path) -> Path:
    """Create a one-pixel processed-IASI pickle with synthetic profiles.

    The spectrum follows IASI's hard-coded 645--2760 cm-1 grid while the
    atmospheric and angular arrays contain only the selected pixel.

    Args:
        tmp_path: Pytest temporary-path fixture.

    Returns:
        Path to the generated processed observation file.
    """
    path = tmp_path / "synthetic_20250323_A.pkl"
    n_points = int(np.floor((2760.0 - 645.0) / 0.25)) + 1
    data = {
        "spec_bbt": np.full((1, n_points), 285.0),
        "zen": np.array([0.0]),
        "sza": np.array([0.0]),
        "saa": np.array([0.0]),
        "azi": np.array([0.0]),
        "ecmwf_z": np.array([[2.0, 1.0, 0.1]]),
        "ecmwf_T": np.array([[278.0, 284.0, 290.0]]),
        "ecmwf_p": np.array([[800.0, 900.0, 1013.0]]),
        "ecmwf_O3": np.array([[0.05, 0.04, 0.03]]),
        "ecmwf_CO": np.array([[0.1, 0.1, 0.1]]),
        "ecmwf_N2O": np.array([[0.3, 0.3, 0.3]]),
        "ecmwf_CO2": np.array([[420.0, 420.0, 420.0]]),
        "ecmwf_CH4": np.array([[1.8, 1.8, 1.8]]),
        "ecmwf_H2O": np.array([[2000.0, 4000.0, 8000.0]]),
    }
    with path.open("wb") as handle:
        pickle.dump(data, handle)
    return path


@pytest.fixture
def tiny_iasi_nedt(tmp_path: Path) -> Path:
    """Create a minimal IASI noise-equivalent-temperature table.

    Three header rows mirror the format consumed by ``iasi_main`` before two
    points bracket the synthetic E2E wavenumber range.

    Args:
        tmp_path: Pytest temporary-path fixture.

    Returns:
        Path to the generated NEDT text file.
    """
    path = tmp_path / "tiny-iasi-nedt.txt"
    path.write_text(
        "# Synthetic IASI NEDT\n# wavenumber temperature\n# cm-1 K\n"
        "999.0 0.2\n1001.0 0.2\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def require_native():
    """Provide a helper for requiring native extensions.

    Unavailable modules cause a skip with an actionable local build command.

    Returns:
        Callable that checks one imported native module.
    """

    def _require(module, name: str) -> None:
        """Skip a test when one requested native module is unavailable.

        Available modules let the calling test continue unchanged.

        Args:
            module: Imported extension module or ``None`` when unavailable.
            name: Human-readable extension name for the skip message.
        """
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

    Args:
        tmp_path: Pytest temporary-path fixture.

    Returns:
        Callable that executes one named native case.
    """

    def _run(case: str, payload: dict):
        """Execute a named native case in an isolated Python subprocess.

        Native termination is diagnosed from both status and result payload.

        Args:
            case: Key selecting an entry point in ``_native_cases.py``.
            payload: Pickle-safe arguments consumed by the selected case.

        Returns:
            A pair containing the unpickled result and subprocess record.
        """
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
