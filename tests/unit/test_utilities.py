from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from srfm import utilities

pytestmark = pytest.mark.unit


def test_closest_uses_joint_coordinate_distance():
    """Verify closest-point selection uses joint coordinate distance.

    Independent axis minima must not be combined into a nonexistent point.
    """
    index = utilities.closest([0, 10, 11], 9, [100, 0, 8], 2)
    assert index == 1


def test_brightness_temperature_inverts_planck_radiance():
    """Verify brightness temperature inverts Planck radiance.

    Representative temperatures and wavenumbers should round-trip accurately.
    """
    h = 6.62607015e-30
    c = 2.99792458e10
    kb = 1.380649e-19
    wavenumber = 1000.0
    temperature = 280.0
    radiance = (
        2 * h * c**2 * wavenumber**3 / np.expm1(h * c * wavenumber / (kb * temperature))
    )
    assert utilities.convert_spectral_radiance_to_bbt(
        radiance, wavenumber
    ) == pytest.approx(temperature, rel=2e-14)


def test_rayleigh_layer_depths_are_nonnegative_and_conserve_total():
    """Verify Rayleigh layer depths are nonnegative and conservative.

    Summed layer optical depth must equal the corresponding total column value.
    """
    lower = pd.Series([100.0, 500.0, 1013.0])
    upper = pd.Series([0.0, 100.0, 500.0])
    result = utilities.calc_Rayleigh_opt_depths(1013.0, lower, upper, 1000.0)
    assert result.shape == (3,)
    assert np.all(result >= 0)
    assert result.sum() == pytest.approx(
        utilities.calc_tot_Rayleigh_opt_depth(1013, 10.0)
    )


def test_rayleigh_numpy_arrays_are_supported():
    """Verify Rayleigh calculations accept NumPy profile arrays.

    Vector inputs should follow the same path as ordinary sequences.
    """
    result = utilities.calc_Rayleigh_opt_depths(
        1000, np.array([500, 1000]), np.array([0, 500]), 500
    )
    assert isinstance(result, np.ndarray)


def test_line_break_respects_limit_and_preserves_tokens():
    """Verify line wrapping respects limits without splitting tokens.

    Rejoining the wrapped text should recover the original token sequence.
    """
    original = "one two three four five six"
    broken = utilities.line_break_str(original, chars=10, delim=" ", indent=2)
    lines = broken.splitlines()
    assert all(len(line) <= 12 for line in lines)
    assert broken.replace("\n", " ").replace("  ", " ").split() == original.split()


def test_line_break_warns_when_no_delimiter_exists():
    """Verify line wrapping warns when text has no delimiter.

    The warning explains why the requested limit cannot be honored safely.
    """
    with pytest.warns(UserWarning, match="Could not split"):
        assert utilities.line_break_str("abcdefghijkl", 5, " ") == "abcdefghijkl"


def test_memory_safe_array_uses_constraints_and_percentage(monkeypatch):
    """Verify memory-safe allocation applies limits and percentages.

    Simulated available memory makes the allocation decision deterministic.

    Args:
        monkeypatch: Pytest fixture replacing available-memory reporting.
    """
    monkeypatch.setattr(
        utilities.psutil, "virtual_memory", lambda: SimpleNamespace(available=8000)
    )
    assert utilities.memory_safe_np_zeros_2d([10], pct=50, max_sec_dim=20).shape == (
        10,
        20,
    )
    assert utilities.memory_safe_np_zeros_2d([10, 50], pct=50).shape == (10, 50)
    with pytest.raises(RuntimeError, match="not enough memory"):
        utilities.memory_safe_np_zeros_2d([10, 51], pct=50)


@pytest.mark.parametrize(
    ("value", "factors"), [(2, [2]), (12, [2, 2, 3]), (17, [17]), (1, [])]
)
def test_prime_factors(value, factors):
    """Verify prime factorization for representative integers.

    Ordered factors are checked for primes, composites, and repeated factors.

    Args:
        value: Parameterized integer to factor.
        factors: Expected ordered prime factors.
    """
    assert utilities.find_prime_factors(value) == factors


def test_layer_extent_and_bounds_round_trip():
    """Verify layer centers and bounds round-trip consistently.

    Converting either representation to the other should recover the inputs.
    """
    upper, lower = utilities.calc_layer_extent(10, 3)
    assert (upper, lower) == (11.5, 8.5)
    assert utilities.calc_layer_bounds(upper, lower) == (10.0, 3.0)


def test_layer_bounds_reject_inverted_extent():
    """Verify inverted layer bounds are rejected.

    An upper boundary below the lower boundary has no physical thickness.
    """
    with pytest.raises(ValueError, match="must"):
        utilities.calc_layer_bounds(1, 2)


def test_add_layer_and_tracking_conversion():
    """Verify adding a layer updates levels and tracking arrays.

    Converting level tracking to layer tracking must preserve placement.
    """
    layer = SimpleNamespace(name="ash", alt_low=1.5, alt_upp=3.5)
    levels, tracking = utilities.add_lyr_from_Layer([0, 1, 2, 3, 4], [None] * 5, layer)
    assert levels == [0, 1, 1.5, 3.5, 4]
    assert tracking == [None, None, "ash", "ash", None]
    assert utilities.track_lev_to_track_lyr(tracking) == [None, None, "ash", None]


@pytest.mark.parametrize(
    ("values", "code"),
    [([1, 2, 3], 1), ([3, 2, 1], 2), ([1, 3, 2], 0), ([1, 1, 2], 0)],
)
def test_monotonic_classification(values, code):
    """Verify monotonic sequences receive the expected classification.

    Increasing, decreasing, constant, and mixed sequences have distinct codes.

    Args:
        values: Parameterized sequence to classify.
        code: Expected monotonicity code.
    """
    assert utilities.monotonic(values) == code


def test_monotonic_rejects_non_vector_and_short_input():
    """Verify monotonicity checks reject unusable inputs.

    Scalars, multidimensional arrays, and undersized sequences are ambiguous.
    """
    with pytest.raises(ValueError, match="one-dimensional"):
        utilities.monotonic([[1, 2]])
    with pytest.raises(ValueError, match="at least two"):
        utilities.monotonic([1])


@pytest.mark.parametrize("units", ["cm-1", "um", "nm"])
def test_spectral_grids_are_round_trip_consistent(units):
    """Verify spectral grids round-trip in every supported unit.

    Wavelength and wavenumber outputs must describe identical samples.

    Args:
        units: Parameterized wavelength or wavenumber unit.
    """
    limits = {"cm-1": (1000, 1002, 1), "um": (10, 12, 1), "nm": (10000, 12000, 1000)}
    wvnm, wavelength = utilities.calc_grids(*limits[units], units)
    np.testing.assert_allclose(wvnm * wavelength, 1e4, rtol=2e-15)
    assert len(wvnm) == 3


@pytest.mark.parametrize(
    "args",
    [
        (1000, 999, 1, "cm-1"),
        (0, 1, 1, "cm-1"),
        (1000, 1002, 0, "cm-1"),
        (1, 2, 1, "Hz"),
    ],
)
def test_invalid_spectral_grids_rejected(args):
    """Verify invalid spectral-grid requests are rejected.

    Reversed bounds, nonpositive resolution, and unsupported units are covered.

    Args:
        args: Parameterized positional arguments for grid construction.
    """
    with pytest.raises(ValueError):
        utilities.calc_grids(*args)


@pytest.mark.parametrize("dist_type", [None, "log_normal"])
def test_mass_loading_number_concentration_round_trip(dist_type):
    """Verify loading and concentration conversions round-trip.

    Both supported distribution types should preserve the original mass.

    Args:
        dist_type: Parameterized particle-size distribution type.
    """
    spread = 1.6
    mass = utilities.mass_loading_from_number_conc(
        25, 2, 2300, spread, dist_type, r=0.5
    )
    number = utilities.number_conc_from_mass_loading(
        mass, 2300, 2, spread, dist_type, r=0.5
    )
    assert number == pytest.approx(25, rel=2e-15)


def test_total_optical_depth_accepts_arrays_series_and_lists():
    """Verify total optical depth accepts common vector containers.

    NumPy arrays, pandas series, and lists should produce equivalent sums.
    """
    result = utilities.calc_tot_dtauc(
        np.array([1, 2]), pd.Series([0.1, 0.2]), [0.01, 0.02]
    )
    np.testing.assert_allclose(result, [1.11, 2.22])


def test_solar_scaling_bounds_and_known_perihelion_amplification():
    """Verify solar scaling bounds and perihelion amplification.

    Seasonal distance correction should remain physical and strongest nearby.
    """
    assert utilities.scale_solar_spectrum(1.0, 4) > 1.0
    for day in (0, 367):
        with pytest.raises(ValueError):
            utilities.scale_solar_spectrum(1.0, day)


def test_altitude_profile_is_increasing_and_isothermal_step_matches_hypsometric_equation():
    """Verify altitude profiles and an isothermal hypsometric step.

    The generated grid must increase and reproduce the analytic special case.
    """
    altitude = utilities.get_altitude_prf([100000, 90000], [280, 280])
    assert altitude[0] == 0
    assert altitude[1] > 0
    expected = -np.log(0.9) * 8.314 * 280 / (28.965 * 9.81)
    assert altitude[1] == pytest.approx(expected)


def test_convolution_preserves_constant_spectrum_away_from_edges(tmp_path):
    """Verify convolution preserves a constant spectrum away from edges.

    Normalized kernel weights should leave interior constant values unchanged.

    Args:
        tmp_path: Pytest temporary-path fixture for the synthetic kernel.
    """
    ils = tmp_path / "ils.txt"
    ils.write_text("-1 0\n0 1\n1 0\n", encoding="utf-8")
    result = utilities.convolve_spectrum(np.ones(9), np.arange(9.0), str(ils))
    np.testing.assert_allclose(result[1:-1], 1.0)
    with pytest.raises(ValueError, match="not regular"):
        utilities.convolve_spectrum(np.ones(4), [0, 1, 3, 4], str(ils))


def test_json_handler_converts_scientific_types():
    """Verify JSON handling converts scientific Python and NumPy types.

    Datetimes, arrays, scalars, and complex values need serializable forms.
    """
    assert utilities.json_handler(np.array([1, 2])) == [1, 2]
    assert utilities.json_handler(np.float64(2.5)) == 2.5
