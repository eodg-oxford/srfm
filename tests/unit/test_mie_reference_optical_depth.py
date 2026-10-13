"""Reference optical-depth inputs and the existing Mie calculation path."""

import numpy as np
import pytest

from srfm import optical_properties
from srfm.layer import MieLayer

pytestmark = pytest.mark.unit


@pytest.fixture
def optical_depth_layer_inputs():
    """Return a complete reference optical-depth layer without material density."""
    return {
        "name": "reference_optical_depth",
        "low_spc": 8.0,
        "upp_spc": 12.0,
        "spec_units": "um",
        "res": 1.0,
        "tau_ref": 0.6,
        "tau_ref_at": 10.0,
        "tau_ref_at_units": "um",
        "r": 0.4,
        "s": 1.6,
        "dist_type": "log_normal",
        "comp": "ri",
        "refractive_index": 1.5 - 0.01j,
        "alt_low": 4.0,
        "alt_upp": 5.5,
        "radii": 24,
        "eta": 1e-6,
        "radii_quad_type": "T",
        "phase_quad_N": 12,
        "phase_quad_type": "L",
        "leg_coeffs": True,
        "leg_coeffs_type": "normalised",
        "multiprocess": False,
    }


@pytest.fixture
def optical_calculation_calls(monkeypatch):
    """Record calls to a linear extinction model at the Mie solver boundary."""
    calls = []

    def calculate_optical_properties(**parameters):
        calls.append(parameters)
        wavelengths = parameters["wavelengths"]
        extinction_coefficient = (
            2e-4 * parameters["distribution"].n * wavelengths / 10.0
        )
        return {
            "beta_ext": extinction_coefficient,
            "ssalb": np.full(wavelengths.size, 0.8),
            "phase_function": np.ones((wavelengths.size, 3)),
            "legendre_coefficient": np.ones((wavelengths.size, 1)),
        }

    monkeypatch.setattr(optical_properties, "ewp_hs", calculate_optical_properties)
    return calls


@pytest.mark.parametrize(
    ("reference_coordinate", "reference_units"),
    [(1000.0, "cm-1"), (10.0, "um"), (10000.0, "nm")],
)
def test_reference_units_feed_the_existing_full_calculation(
    optical_depth_layer_inputs, optical_calculation_calls,
    reference_coordinate, reference_units,
):
    """Every reference unit yields one fixed concentration and a full Mie call."""
    optical_depth_layer_inputs.update(
        tau_ref_at=reference_coordinate, tau_ref_at_units=reference_units
    )
    layer = MieLayer()
    layer.set_input_from_dict(optical_depth_layer_inputs)

    with pytest.warns(UserWarning, match="rho or thickness"):
        layer.calculate_op()
    layer.add_op_calc_output()
    layer.calc_tau()

    reference_call, spectrum_call = optical_calculation_calls
    np.testing.assert_allclose(reference_call["wavelengths"], [10.0])
    assert reference_call["distribution"].n == 1.0
    assert reference_call["multiprocess"] is False
    assert spectrum_call["distribution"].n == pytest.approx(2.0)
    assert spectrum_call["legendre_coefficients_flag"] is True
    np.testing.assert_allclose(spectrum_call["wavelengths"], [8, 9, 10, 11, 12])
    assert layer.tau[2] == pytest.approx(0.6)
    assert layer.n == pytest.approx(2.0)
    assert layer.rho is None
    assert layer.mass_loading is None


@pytest.mark.parametrize("spectral_units", ["um", "cm-1", "nm"])
def test_reference_between_grid_points_interpolates_supplied_refractive_index(
    optical_depth_layer_inputs, optical_calculation_calls, spectral_units,
):
    """Index arrays use wavelength interpolation without rounding the reference."""
    spectral_grids = {
        "um": (8.0, 12.0, 1.0),
        "cm-1": (800.0, 1200.0, 100.0),
        "nm": (8000.0, 12000.0, 1000.0),
    }
    lower_coordinate, upper_coordinate, grid_spacing = spectral_grids[spectral_units]
    optical_depth_layer_inputs.update(
        low_spc=lower_coordinate, upp_spc=upper_coordinate, res=grid_spacing,
        spec_units=spectral_units, tau_ref_at=9.5, rho=2000.0,
    )
    layer = MieLayer(**optical_depth_layer_inputs)
    layer.calc_grids()
    layer.refractive_index = 1.4 + 0.01 * layer.wvls - 0.02j
    original_refractive_indices = layer.refractive_index.copy()

    layer.calculate_op()

    assert optical_calculation_calls[0]["refractive_index"] == pytest.approx(1.495 - 0.02j)
    np.testing.assert_array_equal(optical_calculation_calls[0]["wavelengths"], [9.5])
    np.testing.assert_array_equal(layer.refractive_index, original_refractive_indices)
    np.testing.assert_array_equal(
        optical_calculation_calls[1]["refractive_index"], original_refractive_indices
    )


def test_reference_can_be_outside_calculation_grid_for_constant_index(
    optical_depth_layer_inputs, optical_calculation_calls,
):
    """A reference outside the output band is valid when its index is known."""
    layer = MieLayer(**(optical_depth_layer_inputs | {"tau_ref_at": 0.55, "rho": 2000.0}))

    layer.calculate_op()

    np.testing.assert_array_equal(optical_calculation_calls[0]["wavelengths"], [0.55])
    np.testing.assert_array_equal(optical_calculation_calls[1]["wavelengths"], [8, 9, 10, 11, 12])


@pytest.mark.parametrize(
    ("updates", "error_message"),
    [
        ({"tau_ref": -0.1}, "tau_ref must"),
        ({"tau_ref": 0.0}, "tau_ref must"),
        ({"tau_ref": np.nan}, "tau_ref must"),
        ({"tau_ref": np.inf}, "tau_ref must"),
        ({"tau_ref": True}, "tau_ref must"),
        ({"tau_ref": [0.1]}, "tau_ref must"),
        ({"tau_ref": None}, "tau_ref is required"),
        ({"tau_ref_at": None}, "tau_ref_at must"),
        ({"tau_ref_at": 0.0}, "tau_ref_at must"),
        ({"tau_ref_at": -1.0}, "tau_ref_at must"),
        ({"tau_ref_at": np.nan}, "tau_ref_at must"),
        ({"tau_ref_at": np.inf}, "tau_ref_at must"),
        ({"tau_ref_at": True}, "tau_ref_at must"),
        ({"tau_ref_at_units": None}, "tau_ref_at_units must"),
        ({"tau_ref_at_units": "m"}, "tau_ref_at_units must"),
        ({"tau_ref_at_units": []}, "tau_ref_at_units must"),
        ({"n": 1.0}, "together with n"),
        ({"mass_loading": 0.2}, "together with mass_loading"),
        ({"s_a_den": 0.1}, "together with s_a_den"),
        ({"v_den": 0.1}, "together with v_den"),
        ({"dist_type": "multimode_log_normal", "r": [0.2, 0.5], "s": [1.5, 1.6]},
         "relative mode concentrations"),
    ],
)
def test_reference_input_errors_precede_optical_calculations(
    optical_depth_layer_inputs, optical_calculation_calls, updates, error_message,
):
    """Malformed and ambiguous inputs fail before invoking the native solver."""
    layer = MieLayer(**(optical_depth_layer_inputs | updates))
    with pytest.raises(ValueError, match=error_message):
        layer.calculate_op()
    assert optical_calculation_calls == []


@pytest.mark.parametrize("extinction_coefficient", [0.0, -1.0, np.nan, np.inf])
def test_reference_requires_positive_finite_extinction(
    optical_depth_layer_inputs, monkeypatch, extinction_coefficient,
):
    """An impossible optical-depth normalization fails before the spectrum run."""
    monkeypatch.setattr(
        optical_properties, "ewp_hs",
        lambda **parameters: {"beta_ext": np.array([extinction_coefficient])},
    )
    layer = MieLayer(**optical_depth_layer_inputs)
    with pytest.raises(ValueError, match="Reference extinction coefficient"):
        layer.calculate_op()


@pytest.mark.parametrize(
    ("reference_wavelength", "refractive_indices", "error_message"),
    [(7.0, np.full(5, 1.5 - 0.01j), "outside"),
     (13.0, np.full(5, 1.5 - 0.01j), "outside"),
     (10.0, np.full(4, 1.5 - 0.01j), "match the layer grid")],
)
def test_reference_rejects_missing_index_coverage(
    optical_depth_layer_inputs, optical_calculation_calls,
    reference_wavelength, refractive_indices, error_message,
):
    """Supplied index arrays cannot be silently clipped or reshaped."""
    layer = MieLayer(**(optical_depth_layer_inputs | {
        "tau_ref_at": reference_wavelength, "refractive_index": refractive_indices,
    }))
    with pytest.raises(ValueError, match=error_message):
        layer.calculate_op()
    assert optical_calculation_calls == []


def test_reference_recalculation_retains_input_and_refreshes_derived_quantities(
    optical_depth_layer_inputs, optical_calculation_calls,
):
    """Derived concentrations do not conflict on reruns or mask new inputs."""
    layer = MieLayer(**(optical_depth_layer_inputs | {"rho": 2000.0}))
    layer.calculate_op()
    original_mass_loading = layer.mass_loading
    original_volume_density = layer.v_den
    layer.set_reference_optical_depth(1.2, 1000.0, "cm-1")

    layer.calculate_op()

    assert layer.n == pytest.approx(4.0)
    assert layer.mass_loading == pytest.approx(2 * original_mass_loading)
    assert layer.v_den == pytest.approx(2 * original_volume_density)
    assert layer.tau_ref == 1.2
    layer.rho = None
    with pytest.warns(UserWarning, match="rho or thickness"):
        layer.calculate_op()
    assert layer.mass_loading is None
    layer.set_n(9.0)
    with pytest.raises(ValueError, match="together with n"):
        layer.calculate_op()


def test_reference_uses_resolved_boundary_thickness(
    optical_depth_layer_inputs, optical_calculation_calls,
):
    """Normalization uses the same rounded boundaries as the output optical depth."""
    optical_depth_layer_inputs.pop("alt_low")
    optical_depth_layer_inputs.pop("alt_upp")
    layer = MieLayer(**(optical_depth_layer_inputs | {
        "center_alt": 1.23456, "thick": 0.12345, "rho": 2000.0,
    }))
    layer.calculate_op()
    layer.add_op_calc_output()
    layer.calc_tau()
    assert layer.tau[2] == pytest.approx(layer.tau_ref)


@pytest.mark.parametrize("concentration_field", ["n", "s_a_den", "v_den"])
def test_existing_concentration_inputs_accept_omitted_density(
    optical_depth_layer_inputs, optical_calculation_calls, concentration_field,
):
    """Optional density does not introduce a reference calculation for old inputs."""
    for field in ("tau_ref", "tau_ref_at", "tau_ref_at_units"):
        optical_depth_layer_inputs.pop(field)
    optical_depth_layer_inputs[concentration_field] = 2.0
    layer = MieLayer()
    layer.set_input_from_dict(optical_depth_layer_inputs)
    with pytest.warns(UserWarning, match="rho or thickness"):
        layer.calculate_op()
    assert len(optical_calculation_calls) == 1
    assert layer.n > 0
    assert layer.mass_loading is None


def test_mass_loading_input_still_requires_density(optical_depth_layer_inputs):
    """Density remains mandatory for the existing mass-to-number conversion."""
    for field in ("tau_ref", "tau_ref_at", "tau_ref_at_units"):
        optical_depth_layer_inputs.pop(field)
    layer = MieLayer(**(optical_depth_layer_inputs | {"mass_loading": 0.2}))
    with pytest.raises(ValueError, match="rho is required"):
        layer.validate_inputs()


def test_legacy_concentration_precedence_is_preserved(
    optical_depth_layer_inputs, optical_calculation_calls,
):
    """Reference-input exclusivity does not change existing concentration precedence."""
    for field in ("tau_ref", "tau_ref_at", "tau_ref_at_units"):
        optical_depth_layer_inputs.pop(field)
    layer = MieLayer(**(optical_depth_layer_inputs | {
        "n": 2.0, "s_a_den": 3.0, "v_den": 4.0, "mass_loading": 5.0, "rho": 2000.0,
    }))
    layer.calculate_op()
    assert len(optical_calculation_calls) == 1
    assert layer.n == 2.0
    assert layer.mass_loading != 5.0
