"""Tests for output geometry and multidimensional serialization helpers."""

import numpy as np
import pandas as pd
import pytest
from netCDF4 import Dataset
from pathlib import Path
from unittest.mock import patch

from srfm.forward_model import SRFM
from srfm.main import (
    _calculate_output_utau,
    _create_netcdf_spectral_dimensions,
    _plot_retained_spectral_outputs,
    _plot_spectral_outputs,
    _resolve_output_geometry,
    _resolve_retained_outputs,
    _write_retained_netcdf_outputs,
    _write_srfm_netcdf,
    _write_spectral_text,
)
from srfm.iasi_main import _format_iasi_plot_title

pytestmark = pytest.mark.unit


@pytest.fixture
def optical_profile():
    """Return a top-to-bottom three-layer altitude profile."""
    return pd.DataFrame(
        {
            "h_upper (km)": [3.0, 2.0, 1.0],
            "h_lower (km)": [2.0, 1.0, 0.0],
        }
    )


def test_altitudes_resolve_to_lower_boundaries_and_total_optical_depth(optical_profile):
    """Verify altitude output uses total layer optical depth and adds TOA once."""
    values, rows = _resolve_output_geometry(
        "altitude", [0.9, 3.0], True, optical_profile
    )
    utau = _calculate_output_utau(
        "altitude", values, rows, [0.1, 0.2, 0.3], first_retained_layer=0
    )

    np.testing.assert_allclose(values, [1.0, 3.0])
    assert rows == [1, None]
    np.testing.assert_allclose(utau, [0.3, 0.0])


def test_altitude_optical_depth_accounts_for_truncated_top_layers(optical_profile):
    """Verify discarded negligible top layers do not offset DISORT coordinates."""
    values, rows = _resolve_output_geometry(
        "altitude", [1.0], True, optical_profile
    )
    utau = _calculate_output_utau(
        "altitude", values, rows, [1e-10, 0.2, 0.3], first_retained_layer=1
    )

    np.testing.assert_allclose(utau, [0.2, 0.0])


def test_altitude_output_rejects_values_outside_atmosphere(optical_profile):
    """Verify nearest-level matching cannot fold out-of-range requests inward."""
    with pytest.raises(ValueError, match="within the atmospheric grid"):
        _resolve_output_geometry("altitude", [4.0], False, optical_profile)


def test_tau_output_adds_toa_and_rejects_depth_below_surface(optical_profile):
    """Verify direct optical-depth output adds zero and enforces atmosphere bounds."""
    values, rows = _resolve_output_geometry("tau", [0.2], True, optical_profile)

    assert rows is None
    np.testing.assert_allclose(values, [0.2, 0.0])
    with pytest.raises(ValueError, match="exceeds"):
        _calculate_output_utau("tau", values, rows, [0.1, 0.05], 0)


def test_text_output_writes_correct_multidimensional_indices(tmp_path):
    """Verify flattened text rows retain polar, level, and azimuthal indices."""
    values = np.arange(16.0).reshape(2, 2, 2, 2)
    filename = tmp_path / "spectra.txt"

    _write_spectral_text(filename, np.array([900.0, 901.0]), values, "Radiance")

    written = np.loadtxt(filename)
    assert written.shape == (16, 5)
    np.testing.assert_array_equal(written[:8, 1:4], np.indices((2, 2, 2)).reshape(3, -1).T)
    np.testing.assert_array_equal(written[:, 4], values.reshape(-1))


def test_netcdf_dimensions_follow_native_disort_axis_order(tmp_path):
    """Verify NetCDF labels match SRFM's spectral, polar, level, azimuth order."""
    model = SRFM()
    model.wvnm = np.array([900.0, 901.0])
    model.uu = np.zeros((2, 2, 3, 4))
    model.output_format = "altitude"
    model.output_values = np.array([1.0, 2.0, 3.0])
    model.output_polar_angles = np.array([10.0, 20.0])
    model.output_azimuthal_angles = np.array([0.0, 90.0, 180.0, 270.0])
    filename = tmp_path / "spectra.nc"

    with Dataset(filename, "w") as dataset:
        dimensions = _create_netcdf_spectral_dimensions(dataset, model)
        variable = dataset.createVariable("rad", "f8", dimensions)
        variable[:] = model.uu

    with Dataset(filename) as dataset:
        assert dataset.variables["rad"].dimensions == (
            "wavenumber",
            "output_polar_angle",
            "output_level",
            "output_azimuthal_angle",
        )
        assert dataset.variables["rad"].shape == (2, 2, 3, 4)
        np.testing.assert_allclose(dataset.variables["output_level"][:], [1.0, 2.0, 3.0])
        assert dataset.variables["output_level"].units == "km"


def test_retained_outputs_are_written_with_native_netcdf_dimensions(tmp_path):
    """Explicitly retained DISORT fields use their matching spectral axes."""
    model = SRFM()
    model.wvnm = np.array([900.0, 901.0])
    model.output_format = "altitude"
    model.output_values = np.array([1.0, 2.0, 3.0])
    model.output_polar_angles = np.array([10.0, 20.0])
    model.output_azimuthal_angles = np.array([0.0, 90.0, 180.0, 270.0])
    model.uu = np.arange(48.0).reshape(2, 2, 3, 4)
    model.flup = np.arange(6.0).reshape(2, 3)
    model.albmed = np.arange(4.0).reshape(2, 2)
    filename = tmp_path / "retained.nc"

    with Dataset(filename, "w") as dataset:
        _create_netcdf_spectral_dimensions(dataset, model)
        _write_retained_netcdf_outputs(
            dataset, model, {"flup", "albmed", "radiance"}
        )

    with Dataset(filename) as dataset:
        assert dataset.variables["flup"].dimensions == (
            "wavenumber",
            "output_level",
        )
        assert dataset.variables["albmed"].dimensions == (
            "wavenumber",
            "output_polar_angle",
        )
        assert dataset.variables["uu"].dimensions == (
            "wavenumber",
            "output_polar_angle",
            "output_level",
            "output_azimuthal_angle",
        )
        np.testing.assert_allclose(dataset.variables["flup"][:], model.flup)
        np.testing.assert_allclose(dataset.variables["albmed"][:], model.albmed)
        np.testing.assert_allclose(dataset.variables["uu"][:], model.uu)


def test_unspecified_retained_outputs_do_not_change_netcdf_contents(tmp_path):
    """The historical default retention does not silently enlarge output files."""
    model = SRFM()
    model.wvnm = np.array([900.0])
    model.uu = np.zeros((1, 1, 1, 1))
    model.flup = np.ones((1, 1))
    model.output_format = "tau"
    model.output_values = np.array([0.0])
    model.output_polar_angles = np.array([0.0])
    model.output_azimuthal_angles = np.array([0.0])
    filename = tmp_path / "default.nc"

    with Dataset(filename, "w") as dataset:
        _create_netcdf_spectral_dimensions(dataset, model)
        _write_retained_netcdf_outputs(dataset, model, None)

    with Dataset(filename) as dataset:
        assert "flup" not in dataset.variables


def test_complete_netcdf_writer_supports_raw_output_without_radiance_or_bbt(
    tmp_path,
):
    """A raw retained DISORT field can be the only result in the generic file."""
    model = SRFM()
    model.wvnm = np.array([900.0, 901.0])
    model.output_format = "altitude"
    model.output_values = np.array([1.0, 2.0, 3.0])
    model.output_polar_angles = np.array([0.0])
    model.output_azimuthal_angles = np.array([0.0])
    model.rfldn = np.arange(6.0).reshape(2, 3)
    filename = tmp_path / "raw-only.nc"

    _write_srfm_netcdf(
        filename,
        model,
        ("rfldn",),
        {"driver_inputs": {"spectral": (900.0, 901.0)}},
        np.array([10.0, 11.0]),
        {},
    )

    with Dataset(filename) as dataset:
        result_names = {
            "rfldir",
            "rfldn",
            "flup",
            "dfdt",
            "uavg",
            "uu",
            "albmed",
            "trnmed",
            "bbt",
        }
        assert result_names.intersection(dataset.variables) == {"rfldn"}
        assert dataset.variables["rfldn"].dimensions == (
            "wavenumber",
            "output_level",
        )
        np.testing.assert_allclose(dataset.variables["rfldn"][:], model.rfldn)


def test_retained_output_selection_accepts_every_radiance_alias():
    """The removed rad boolean is replaced by equivalent retention aliases."""
    for alias in ("rad", "radiance", "uu"):
        requested, runtime = _resolve_retained_outputs(
            {
                "retain_outputs": ("bbt", alias, "flup"),
                "convolve_iasi": False,
            }
        )
        assert requested == {"bbt", "uu", "flup"}
        assert runtime == {"uu", "flup"}


def test_base_plots_include_every_retained_primary_output(tmp_path):
    """Retaining both BBT and radiance creates two non-overwriting plot sets."""
    model = SRFM()
    with patch("srfm.main._plot_spectral_outputs") as plot:
        _plot_retained_spectral_outputs(
            model, tmp_path, ("bbt", "radiance"), False
        )

    assert [call.args[2] for call in plot.call_args_list] == ["bbt", "rad"]
    assert [call.kwargs["filename_prefix"] for call in plot.call_args_list] == [
        "base_plot_bbt",
        "base_plot_rad",
    ]


def test_base_plots_warn_when_no_primary_output_is_retained(tmp_path):
    """Flux-only retention cannot create a BBT or radiance base plot."""
    with pytest.warns(UserWarning, match="no base plots were created"):
        _plot_retained_spectral_outputs(SRFM(), tmp_path, ("flup",), False)


def test_spectral_plot_titles_use_angles_but_filenames_use_indices(tmp_path):
    """Plot titles use configured angles while filenames retain array indices."""
    model = SRFM()
    model.wvnm = np.array([900.0, 901.0])
    model.bbt = np.zeros((2, 2, 1, 2))
    model.output_format = "altitude"
    model.output_values = np.array([3.0])
    model.output_polar_angles = np.array([15.0, 42.5])
    model.output_azimuthal_angles = np.array([0.0, 135.0])

    with patch("srfm.main.plt.title") as title, patch(
        "srfm.main.plt.savefig"
    ) as savefig, patch("srfm.main.plt.close"):
        _plot_spectral_outputs(model, tmp_path, "bbt", False)

    assert title.call_count == 4
    assert title.call_args_list[0].args[0] == (
        "Output polar angle: 15°\n"
        "Output altitude (km): 3\n"
        "Output azimuthal angle: 0°"
    )
    assert title.call_args_list[-1].args[0] == (
        "Output polar angle: 42.5°\n"
        "Output altitude (km): 3\n"
        "Output azimuthal angle: 135°"
    )
    filenames = [Path(call.args[0]).name for call in savefig.call_args_list]
    assert filenames == [
        "base_plot_0_0_0.png",
        "base_plot_0_0_1.png",
        "base_plot_1_0_0.png",
        "base_plot_1_0_1.png",
    ]


def test_iasi_plot_title_uses_angle_values():
    """The legacy IASI title reports physical angles, not array indices."""
    model = SRFM()
    model.output_format = "altitude"
    model.output_values = np.array([3.0])
    model.output_polar_angles = np.array([27.5])
    model.output_azimuthal_angles = np.array([142.0])

    assert _format_iasi_plot_title("iasi", model, 0, 0, 0) == (
        "iasi\nOutput polar angle: 27.5°; "
        "altitude (km): 3; azimuthal angle: 142°"
    )
