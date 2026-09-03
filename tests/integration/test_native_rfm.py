import numpy as np
import pytest

from srfm import rfm_helper
from srfm.RFM import rfm_py

pytestmark = pytest.mark.integration


def test_structured_driver_runs_real_rfm_and_captures_optical_depths(
    tmp_path, tiny_atmosphere, tiny_altitude_grid, tiny_xsc_file, require_native
):
    """Verify a structured driver through the compiled RFM boundary.

    The result is checked for parsed differential optical depths on the tiny
    synthetic atmosphere and cross-section grid.

    Args:
        tmp_path: Pytest temporary-path fixture.
        tiny_atmosphere: Synthetic atmospheric-profile fixture.
        tiny_altitude_grid: Synthetic altitude-grid fixture.
        tiny_xsc_file: Synthetic F11 cross-section fixture.
        require_native: Fixture helper for extension availability.
    """
    require_native(rfm_py, "RFM")
    result = rfm_helper.run_rfm_with_parameters(
        header="pytest self-contained RFM",
        flags=("OPT", "NAD", "SFC", "PRF", "LEV", "DBL"),
        spectral=(rfm_helper.SpectralRange(999.5, 1000.5, 0.5),),
        gases=("F11",),
        atmosphere=(tiny_altitude_grid, tiny_atmosphere),
        tangent=("1.0",),
        lev=("0", "1", "2"),
        sfc=("TEMSFC=290",),
        xsc=(tiny_xsc_file,),
        run_id="pytest",
        directory=tmp_path,
        output_mode="capture",
        enable_capture=True,
        optical_levels=(0, 1, 2),
    )
    assert result.ok
    assert result.output_df.shape == (2, 16)
    differential = result.output_df.filter(regex=r"^dOD_").to_numpy()
    assert differential.shape == (2, 3)
    assert np.all(np.isfinite(differential))
    assert np.all(differential >= 0)
    np.testing.assert_allclose(
        differential,
        [[0.026936, 0.035914, 0.026936], [0.034026, 0.045368, 0.034026]],
        rtol=3e-4,
        atol=2e-7,
    )
    integrated = result.output_df.filter(regex=r"^iOD_").to_numpy()
    np.testing.assert_allclose(integrated[-1], differential.sum(axis=0), rtol=2e-6)


def test_real_rfm_failure_emits_log_to_stderr_and_removes_file(
    tmp_path, require_native, capsys
):
    """Verify native RFM failures emit and clean their transient log.

    This protects actionable scheduler diagnostics without leaving stale files
    that could contaminate subsequent runs.

    Args:
        tmp_path: Pytest temporary-path fixture.
        require_native: Fixture helper for extension availability.
        capsys: Pytest fixture capturing standard streams.
    """
    require_native(rfm_py, "RFM")
    missing_atmosphere = tmp_path / "deliberately-missing.atm"
    driver_lines = (
        "*HDR",
        "  pytest deliberate RFM failure",
        "*FLG",
        "  OPT NAD SFC PRF LEV DBL",
        "*SPC",
        "  999.5 1000.5 0.5",
        "*GAS",
        "  F11",
        "*ATM",
        f"  {missing_atmosphere}",
        "*TAN",
        "  1.0",
        "*LEV",
        "  0.0 1.0",
        "*SFC",
        "  TEMSFC=290",
        "*END",
    )

    with pytest.raises(RuntimeError, match="RFM"):
        rfm_helper.run_rfm(
            directory=tmp_path,
            driver_lines=driver_lines,
            output_mode="capture",
            enable_capture=True,
            optical_levels=(0.0, 1.0),
        )

    streams = capsys.readouterr()
    assert "--- rfm.log ---" in streams.err
    assert "pytest deliberate RFM failure" in streams.err
    assert not list(tmp_path.glob("rfm.log*"))
