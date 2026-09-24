"""Tests for common in-memory and file-backed spectral fields."""

from __future__ import annotations

import numpy as np
import pytest
from unittest.mock import Mock, patch

from srfm.main import (
    _prepare_boundary_spectral_fields,
    _prepare_solar_spectral_irradiance,
    _set_spectral_boundary_inputs,
)
from srfm import utilities
from srfm.spectral_fields import SpectralField, convert_solar_spectral_density

pytestmark = pytest.mark.unit


def test_spectral_field_interpolates_in_canonical_wavenumber_order():
    """Descending wavelength inputs are converted, sorted, and interpolated."""
    field = SpectralField.from_specification(
        {
            "grid": np.array([10.0, 8.0, 5.0]),
            "values": np.array([0.1, 0.3, 0.7]),
            "grid_units": "um",
        },
        "albedo",
        minimum=0,
        maximum=1,
    )

    np.testing.assert_allclose(field.wavenumber_cm_inverse, [1000.0, 1250.0, 2000.0])
    np.testing.assert_allclose(field.interpolate([1125.0, 1625.0]), [0.2, 0.5])


@pytest.mark.parametrize(
    ("units", "values", "expected"),
    [
        ("W m-2 (cm-1)-1", [2.0, 3.0], [2.0, 3.0]),
        ("W m-2 um-1", [200.0, 300.0], [2.0, 0.75]),
        ("W m-2 nm-1", [0.2, 0.3], [2.0, 0.75]),
        ("mW cm-2 um-1", [20.0, 30.0], [2.0, 0.75]),
    ],
)
def test_solar_density_unit_conversions(units, values, expected):
    """Every documented solar unit spelling includes the spectral Jacobian."""
    converted = convert_solar_spectral_density(
        np.asarray(values), units, np.array([1000.0, 2000.0])
    )
    np.testing.assert_allclose(converted, expected)


def test_file_backed_scalar_and_matrix_fields(tmp_path):
    """Text columns support scalar properties and multiple phase moments."""
    filename = tmp_path / "properties.txt"
    filename.write_text(
        "# wn tau beta0 beta1\n1000 0.1 1.0 0.2\n2000 0.3 1.0 0.4\n",
        encoding="utf-8",
    )
    optical_depth = SpectralField.from_specification(
        {
            "file": filename,
            "grid_column": 0,
            "value_column": 1,
            "skiprows": 1,
            "grid_units": "cm-1",
        },
        "tau",
        minimum=0,
    )
    moments = SpectralField.from_specification(
        {
            "file": filename,
            "grid_column": 0,
            "value_columns": [2, 3],
            "skiprows": 1,
            "grid_units": "cm-1",
        },
        "moments",
        matrix_values=True,
    )

    np.testing.assert_allclose(optical_depth.interpolate([1500]), [0.2])
    np.testing.assert_allclose(moments.interpolate([1500]), [[1.0, 0.3]])


def test_relative_spectral_file_paths_use_current_working_directory(
    tmp_path, monkeypatch
):
    """File paths do not depend on whether inputs came from a driver module."""
    filename = tmp_path / "albedo.txt"
    filename.write_text("1000 0.1\n2000 0.3\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    field = SpectralField.from_specification(
        {
            "file": "albedo.txt",
            "grid_column": 0,
            "value_column": 1,
            "grid_units": "cm-1",
        },
        "albedo",
    )

    np.testing.assert_allclose(field.interpolate([1500]), [0.2])


def test_spectral_field_rejects_invalid_values_and_incomplete_coverage():
    """Physical bounds and complete-grid coverage are enforced explicitly."""
    with pytest.raises(ValueError, match="less than or equal to 1"):
        SpectralField.from_specification(
            {"grid": [1000, 2000], "values": [0.2, 1.1], "grid_units": "cm-1"},
            "albedo",
            minimum=0,
            maximum=1,
        )

    field = SpectralField.from_specification(
        {"grid": [1000, 2000], "values": [0.2, 0.3], "grid_units": "cm-1"},
        "albedo",
    )
    with pytest.raises(ValueError, match="computational grid spans"):
        field.validate_coverage([999, 1500])


def test_boundary_preparation_preserves_scalar_and_converts_custom_solar():
    """Scalar albedo remains a fast path and custom FBEAM is not normalized."""
    albedo, solar = _prepare_boundary_spectral_fields(
        {
            "albedo": 0.25,
            "solar_spectrum": {
                "grid": [1000.0, 2000.0],
                "values": [4.0, 8.0],
                "grid_units": "cm-1",
                "value_units": "W m-2 (cm-1)-1",
            },
        },
        np.array([1000.0, 1500.0, 2000.0]),
    )

    assert albedo == 0.25
    np.testing.assert_allclose(solar.interpolate([1000, 1500, 2000]), [4, 6, 8])


def test_custom_solar_is_not_date_scaled_or_amplitude_normalized():
    """Custom FBEAM depends only on unit conversion and linear interpolation."""
    _, custom = _prepare_boundary_spectral_fields(
        {
            "albedo": 0.0,
            "solar_spectrum": {
                "grid": [1000.0, 2000.0],
                "values": [4.0, 8.0],
                "grid_units": "cm-1",
                "value_units": "W m-2 (cm-1)-1",
            },
        },
        [1000.0, 1500.0, 2000.0],
    )
    with patch(
        "srfm.main.utilities.load_solar_spectrum_Gueymard20018"
    ) as load, patch("srfm.main.utilities.scale_solar_spectrum") as scale:
        january = _prepare_solar_spectral_irradiance(
            custom, [1000.0, 1500.0, 2000.0], 1, sun=True
        )
        july = _prepare_solar_spectral_irradiance(
            custom, [1000.0, 1500.0, 2000.0], 182, sun=True
        )

    np.testing.assert_allclose(january, [4.0, 6.0, 8.0])
    np.testing.assert_array_equal(january, july)
    load.assert_not_called()
    scale.assert_not_called()


def test_omitted_custom_solar_preserves_legacy_calculation_exactly():
    """The helper retains the established Gueymard interpolation/date pathway."""
    target = np.array([1000.0, 1001.0])
    source_values, source_grid = utilities.load_solar_spectrum_Gueymard20018()
    expected = utilities.scale_solar_spectrum(
        np.interp(target, source_grid[::-1], source_values[::-1]), 91
    )

    actual = _prepare_solar_spectral_irradiance(None, target, 91, sun=True)

    np.testing.assert_array_equal(actual, expected)


def test_disabled_sun_returns_zero_without_loading_any_spectrum():
    """Sun-off runs pass zero FBEAM even when a custom spectrum was supplied."""
    _, custom = _prepare_boundary_spectral_fields(
        {
            "albedo": 0.0,
            "solar_spectrum": {
                "grid": [1000.0, 2000.0],
                "values": [4.0, 8.0],
                "grid_units": "cm-1",
                "value_units": "W m-2 (cm-1)-1",
            },
        },
        [1000.0, 2000.0],
    )
    with patch("srfm.spectral_fields.SpectralField.interpolate") as interpolate, patch(
        "srfm.main.utilities.load_solar_spectrum_Gueymard20018"
    ) as load:
        irradiance = _prepare_solar_spectral_irradiance(
            custom, [1000.0, 1500.0, 2000.0], 1, sun=False
        )

    np.testing.assert_array_equal(irradiance, 0.0)
    interpolate.assert_not_called()
    load.assert_not_called()


def test_disort_receives_per_wavenumber_albedo_and_unscaled_fbeam():
    """The spectral loop boundary helper forwards exact prepared values."""
    disort = Mock()
    albedo = np.array([0.1, 0.3])
    irradiance = np.array([4.0, 8.0])

    _set_spectral_boundary_inputs(
        disort, 1, albedo, irradiance, sun=True
    )

    disort.set_albedo.assert_called_once_with(0.3)
    disort.set_fbeam.assert_called_once_with(8.0)

    disort.reset_mock()
    _set_spectral_boundary_inputs(
        disort, 0, albedo, irradiance, sun=False
    )
    disort.set_albedo.assert_called_once_with(0.1)
    disort.set_fbeam.assert_called_once_with(0)


def test_solar_zenith_changes_umu0_without_rescaling_custom_fbeam():
    """Beam-normal custom irradiance is independent of the incidence cosine."""
    received = []
    for solar_zenith_degrees in (0.0, 60.0):
        disort = Mock()
        disort.set_umu0(np.cos(np.deg2rad(solar_zenith_degrees)).item())
        _set_spectral_boundary_inputs(
            disort,
            0,
            0.0,
            np.array([7.5]),
            sun=True,
        )
        received.append(
            (
                disort.set_umu0.call_args.args[0],
                disort.set_fbeam.call_args.args[0],
            )
        )

    assert received == [(1.0, 7.5), (pytest.approx(0.5), 7.5)]
