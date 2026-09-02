from pprint import pformat

import numpy as np
import pytest
from netCDF4 import Dataset

import srfm
from srfm import rfm_helper
from srfm.DISORT import disort_module_s
from srfm.DISORT_dbl import disort_module_d
from srfm.RFM import rfm_py
from srfm.inputs import Inputs


pytestmark = pytest.mark.e2e


# This documents every input that the current run_srfm implementation reads.
# A complete driver-table test below asserts that none is accidentally omitted.
COMPLETE_RUN_INPUT_KEYS = {
    "adjust_maxcmu",
    "albedo",
    "azi",
    "base_plots",
    "bbt",
    "bbt_out_fname",
    "btemp",
    "convolve_iasi",
    "date",
    "deltamplus",
    "disort_precision",
    "do_pseudo_sphere",
    "driver_inputs",
    "earth_radius",
    "fin_res",
    "fin_wvnmhi",
    "fin_wvnmlo",
    "fisot",
    "header",
    "iasi_ils",
    "ibcnd",
    "lamber",
    "levels",
    "maxcmu",
    "maxphi",
    "maxulv",
    "maxumu",
    "nmom",
    "onlyfl",
    "out_mode",
    "planck",
    "plot_type",
    "prnt",
    "rad",
    "rad_out_fname",
    "results_fldr",
    "rfm_config",
    "saa",
    "scat_lyrs_inputs",
    "show_plots",
    "spc_res",
    "spc_units",
    "spc_wvnmhi",
    "spc_wvnmlo",
    "sun",
    "sza",
    "temis",
    "ttemp",
    "usrang",
    "usrtau",
    "utau",
    "zen",
}


def _complete_input_values(
    results, tiny_atmosphere, tiny_altitude_grid, tiny_xsc_file
):
    """Return a complete, tiny equivalent of the basic-example driver inputs."""
    return {
        "fwd_model": "SRFM",
        "instrument": "synthetic",
        "plot_profiles": False,
        "results_fldr": str(results),
        "date": (2025, 3, 23),
        "fin_wvnmlo": 999.5,
        "fin_wvnmhi": 1000.5,
        "fin_res": 0.5,
        "spc_wvnmlo": 999.5,
        "spc_wvnmhi": 1000.5,
        "spc_res": 0.5,
        "spc_units": "cm-1",
        "rfm_config": {
            "output_mode": "capture",
            "driver_path": None,
            "generate_driver": False,
            "verbose": False,
            "capture_files_content": False,
        },
        "driver_inputs": {
            "header": "pytest self-contained SRFM E2E",
            "flags": ("OPT", "NAD", "SFC", "PRF", "LEV", "DBL"),
            "spectral": (rfm_helper.SpectralRange(999.5, 1000.5, 0.5),),
            "gases": ("F11",),
            "atmosphere": (str(tiny_altitude_grid), str(tiny_atmosphere)),
            "xsc": (str(tiny_xsc_file),),
        },
        "levels": [0.0, 1.0, 2.0],
        "scat_lyrs_inputs": {
            "synthetic_aerosol": {
                "name": "synthetic_aerosol",
                "low_spc": 999.5,
                "upp_spc": 1000.5,
                "res": 0.5,
                "spec_units": "cm-1",
                "mass_loading": 0.005,
                "n": None,
                "r": 0.2,
                "s": 1.4,
                "rho": 1800.0,
                "s_a_den": None,
                "v_den": None,
                "dist_type": "log_normal",
                "comp": "ri",
                "refractive_index": 1.50 - 0.005j,
                "center_alt": 0.5,
                "thick": 0.2,
                "alt_upp": None,
                "alt_low": None,
                "radii": 10,
                "eta": 1e-5,
                "phase_quad_N": 16,
                "phase_quad_type": "L",
                "radii_quad_type": "T",
                "leg_coeffs": True,
                "leg_coeffs_type": "normalised",
                "multiprocess": False,
            }
        },
        "nmom": 8,
        "maxcmu": 4,
        "maxumu": 1,
        "maxphi": 1,
        "maxulv": 1,
        "usrang": True,
        "usrtau": True,
        "ibcnd": 0,
        "onlyfl": False,
        "prnt": [False] * 5,
        "planck": True,
        "lamber": True,
        "deltamplus": False,
        "do_pseudo_sphere": False,
        "utau": [0.0],
        "fisot": 0.0,
        "albedo": 0.0,
        "temis": 1.0,
        "earth_radius": 6371.0,
        "header": "NO HEADER",
        "btemp": 290.0,
        "ttemp": 278.0,
        "disort_precision": "double",
        "adjust_maxcmu": False,
        "sun": False,
        "sza": 0.0,
        "saa": 0.0,
        "zen": 0.0,
        "azi": 0.0,
        "convolve_iasi": False,
        "iasi_ils": None,
        "out_mode": None,
        "rad": False,
        "bbt": False,
        "rad_out_fname": None,
        "bbt_out_fname": None,
        "base_plots": False,
        "show_plots": False,
        "plot_type": "bbt",
    }


def _require_all_native_extensions(require_native, precision="double"):
    require_native(rfm_py, "RFM")
    require_native(srfm.mie_module, "Mie")
    if precision == "single":
        require_native(disort_module_s, "single-precision DISORT")
    else:
        require_native(disort_module_d, "double-precision DISORT")


def _assert_complete_result(result, results, values):
    assert result is not None
    n_points = (
        int(
            np.floor(
                (values["fin_wvnmhi"] - values["fin_wvnmlo"])
                / values["fin_res"]
            )
        )
        + 1
    )
    expected_grid = values["fin_wvnmlo"] + np.arange(n_points) * values["fin_res"]
    np.testing.assert_allclose(result["wvnm"], expected_grid)
    assert result["uu"].shape == (n_points, 1, 1, 1)
    assert result["bbt"].shape == (n_points, 1, 1, 1)
    assert np.all(np.isfinite(result["uu"]))
    assert np.all(np.isfinite(result["bbt"]))
    assert np.all(result["uu"] > 0)
    assert np.all((result["bbt"] > 150) & (result["bbt"] < 350))
    assert (results / "rfm_files" / "grid.spc").is_file()


def _write_driver(path, values):
    path.write_text(
        "from srfm.rfm_helper import SpectralRange\n\n"
        f"inputs = {pformat(values, sort_dicts=False)}\n",
        encoding="utf-8",
    )


def test_complete_srfm_pathway_with_self_contained_scientific_inputs(
    tmp_path,
    tiny_atmosphere,
    tiny_altitude_grid,
    tiny_xsc_file,
    require_native,
    run_native_case,
):
    """Exercise direct Inputs -> RFM -> Mie -> DISORT -> SRFM."""
    _require_all_native_extensions(require_native)
    results = tmp_path / "direct-results"
    values = _complete_input_values(
        results, tiny_atmosphere, tiny_altitude_grid, tiny_xsc_file
    )

    result, _ = run_native_case("e2e", {"values": values})

    _assert_complete_result(result, results, values)


def test_complete_driver_table_is_read_and_executed(
    tmp_path,
    tiny_atmosphere,
    tiny_altitude_grid,
    tiny_xsc_file,
    require_native,
    run_native_case,
):
    """Exercise a basic-example-style driver table through the full native model."""
    _require_all_native_extensions(require_native)
    results = tmp_path / "driver-results"
    values = _complete_input_values(
        results, tiny_atmosphere, tiny_altitude_grid, tiny_xsc_file
    )
    driver = tmp_path / "driver_table.py"
    _write_driver(driver, values)

    loaded = Inputs()
    loaded.read_srfm_drv(driver)

    assert COMPLETE_RUN_INPUT_KEYS <= loaded.values.keys()
    assert loaded.values["driver_inputs"]["atmosphere"] == (
        str(tiny_altitude_grid),
        str(tiny_atmosphere),
    )
    assert loaded.values["driver_inputs"]["xsc"] == (str(tiny_xsc_file),)

    result, completed = run_native_case("e2e", {"driver_path": str(driver)})

    _assert_complete_result(result, results, values)
    assert COMPLETE_RUN_INPUT_KEYS <= set(result["input_keys"])
    assert "--- rfm.log ---" in completed.stdout
    assert "R-RFM: Running RFM" in completed.stdout
    assert not list(results.glob("rfm.log*"))


@pytest.mark.parametrize(
    ("scenario", "updates", "expected_files"),
    [
        pytest.param(
            "single_txt_rad_plot",
            {
                "disort_precision": "single",
                "out_mode": "txt",
                "bbt": True,
                "rad": True,
                "bbt_out_fname": None,
                "rad_out_fname": "single-radiance",
                "base_plots": True,
                "plot_type": "rad",
            },
            ("bbt.txt", "single-radiance.txt", "base_plot.png"),
            id="single-txt-radiance-plot",
        ),
        pytest.param(
            "solar_iasi_netcdf_bbt_plot",
            {
                "spc_res": 0.25,
                "fin_res": 0.25,
                "sun": True,
                "sza": 30.0,
                "saa": 45.0,
                "convolve_iasi": True,
                "out_mode": "netcdf",
                "bbt": True,
                "rad": True,
                "bbt_out_fname": "solar-bbt",
                "rad_out_fname": None,
                "base_plots": True,
                "plot_type": "bbt",
            },
            ("solar-bbt.nc", "rad.nc", "base_plot.png"),
            id="solar-iasi-netcdf-bbt-plot",
        ),
    ],
)
def test_driver_table_optional_execution_branches(
    scenario,
    updates,
    expected_files,
    tmp_path,
    tiny_atmosphere,
    tiny_altitude_grid,
    tiny_xsc_file,
    tiny_iasi_ils,
    require_native,
    run_native_case,
):
    """Exercise mutually exclusive output and native-option branches end to end."""
    results = tmp_path / scenario
    values = _complete_input_values(
        results, tiny_atmosphere, tiny_altitude_grid, tiny_xsc_file
    )
    values.update(updates)
    if values["convolve_iasi"]:
        values["iasi_ils"] = str(tiny_iasi_ils)
    _require_all_native_extensions(require_native, values["disort_precision"])

    driver = tmp_path / f"{scenario}_driver.py"
    _write_driver(driver, values)
    result, completed = run_native_case("e2e", {"driver_path": str(driver)})

    _assert_complete_result(result, results, values)
    assert "--- rfm.log ---" in completed.stdout
    for filename in expected_files:
        path = results / filename
        assert path.is_file()
        assert path.stat().st_size > 0

    if values["out_mode"] == "txt":
        bbt = np.loadtxt(results / "bbt.txt")
        radiance = np.loadtxt(results / "single-radiance.txt")
        assert bbt.shape == radiance.shape == (3, 2)
        np.testing.assert_allclose(bbt[:, 0], result["wvnm"])
        np.testing.assert_allclose(radiance[:, 0], result["wvnm"])
    else:
        with Dataset(results / "solar-bbt.nc") as dataset:
            assert dataset.variables["bbt"].shape == (5,)
            assert np.all(np.isfinite(dataset.variables["bbt"][:]))
        with Dataset(results / "rad.nc") as dataset:
            assert dataset.variables["rad"].shape == (5,)
            assert np.all(np.isfinite(dataset.variables["rad"][:]))
