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


def test_reconciled_collision_runs_native_rfm_and_returns_effective_levels(
    tmp_path, tiny_xsc_file, require_native
):
    """Prevent the 3/3.001 km ADDLEV collision through real RFM capture."""
    from srfm.rfm_functions import write_atm_file
    from srfm.rfm_levels import read_matching_grid, reconcile_levels

    require_native(rfm_py, "RFM")
    atmosphere = tmp_path / "collision.atm"
    write_atm_file({"HGT [km]": [2, 3, 4], "PRE [mb]": [1013, 900, 800],
                    "TEM [K]": [290, 284, 278], "F11 [ppmv]": [1, 1, 1]}, atmosphere)
    flags = ("OPT", "NAD", "SFC", "PRF", "LEV", "DBL")
    grid = read_matching_grid({"atmosphere": (atmosphere,), "flags": flags}, {}, tmp_path)
    levels, tracking, boundaries = reconcile_levels(grid, [2, 3, 4], {"aerosol": (3.001, 3.101)})
    result = rfm_helper.run_rfm_with_parameters(
        header="pytest reconciled collision", flags=flags,
        spectral=(rfm_helper.SpectralRange(999.5, 1000.5, 0.5),),
        gases=("F11",), atmosphere=(atmosphere,), tangent=("1",),
        lev=tuple(map(str, levels)), sfc=("TEMSFC=290",), xsc=(tiny_xsc_file,),
        directory=tmp_path, output_mode="capture", optical_levels=levels,
        optical_depth_format="compact",
    )
    captured = result.optical_depth_grid
    assert result.ok
    np.testing.assert_array_equal(captured.altitude_lower[::-1], levels[:-1])
    np.testing.assert_array_equal(captured.altitude_upper[::-1], levels[1:])
    assert tracking == [None, "aerosol", None]
    assert boundaries["aerosol"] == (3, float(np.float32(3.101)))
    assert np.all(np.isfinite(captured.differential_tau))
    assert np.all(captured.differential_tau >= 0)


@pytest.mark.parametrize("subdivide", [False, True])
def test_existing_driver_prelev_pressure_observer_matches_bundled_rfm(
    tmp_path, tiny_xsc_file, require_native, capsys, subdivide
):
    """Check the inspected pressure path and in-memory LEV replacement natively."""
    from srfm.rfm_levels import read_matching_grid, reconcile_levels

    require_native(rfm_py, "RFM")
    path = tmp_path / "existing.drv"
    # A pressure computational grid initializes PRESFC before *OBS. With a
    # height grid the bundled RFM leaves it unset until later initialization.
    pressure_atmosphere = tmp_path / "pressure.atm"
    pressure_atmosphere.write_text(
        "3\n*PRE [mb]\n1013 900 800\n*HGT [km]\n0 1 2\n"
        "*TEM [K]\n290 284 278\n*F11 [ppmv]\n1 1 1\n*END\n"
    )
    content = (
        "*HDR\npytest pressure observer\n*FLG\nOPT NAD SFC PRF LEV DBL OBS"
        + (" LAY" if subdivide else "") + "\n"
        "*SPC\n999.5 1000.5 0.5\n*GAS\nF11\n"
        f"*ATM\n{pressure_atmosphere}\n*TAN\n1\n*OBS\nPREOBS=810\n"
        "*LEV\n0 1 1.0008 2\n*SFC\nTEMSFC=290\n"
        f"*XSC\n{tiny_xsc_file}\n*END\n"
    )
    path.write_text(content)
    configuration = {"driver_path": path, "verbose": False, "capture_files_content": False}
    grid = read_matching_grid({}, configuration, tmp_path)
    levels, _, _ = reconcile_levels(grid, [*grid, 1.0008], {})
    result = rfm_helper.rfm_main(
        configuration=configuration,
        driver_inputs={"header": "unused", "flags": ("OPT", "NAD", "LEV"),
                       "spectral": (), "gases": (), "atmosphere": (), "tangent": ("1",)},
        levels=levels, rfm_out_fldr=tmp_path,
    )
    assert result.ok
    np.testing.assert_array_equal(result.optical_depth_grid.altitude_lower[::-1], levels[:-1])
    assert "Change pressure" in capsys.readouterr().out
    assert path.read_text() == content
