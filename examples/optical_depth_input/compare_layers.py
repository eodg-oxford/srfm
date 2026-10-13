"""Recover a randomly chosen Mie layer from one reference optical depth.

The fixed random seed makes the example reproducible. Only the compiled Mie
extension is needed; no atmospheric profiles or external spectroscopy are used.
"""

import json

import numpy as np

from srfm.layer import MieLayer


def compare_layers(seed=20261012):
    """Calculate and compare concentration and optical-depth input examples.

    Args:
        seed (int): Seed used to choose physical particle and layer parameters.

    Returns:
        dict: Physical inputs, reference optical depth, recovered concentration,
        and maximum errors between the complete calculated spectra.

    Raises:
        AssertionError: If any recovered quantity differs beyond numerical tolerance.
    """
    random_generator = np.random.default_rng(seed)
    number_concentration = float(random_generator.uniform(5.0, 250.0))
    median_radius_micrometres = float(random_generator.uniform(0.2, 1.8))
    geometric_standard_deviation = float(random_generator.uniform(1.2, 1.9))
    particle_density_kg_per_cubic_metre = float(random_generator.uniform(1000.0, 2500.0))
    lower_altitude_kilometres = float(np.round(random_generator.uniform(2.0, 8.0), 3))
    thickness_kilometres = float(np.round(random_generator.uniform(0.2, 2.0), 3))
    refractive_index = complex(
        random_generator.uniform(1.3, 1.7), -random_generator.uniform(0.001, 0.05)
    )
    layer_inputs = {
        "name": "random_particle_layer",
        "low_spc": 8.0,
        "upp_spc": 13.0,
        "res": 0.25,
        "spec_units": "um",
        "n": number_concentration,
        "r": median_radius_micrometres,
        "s": geometric_standard_deviation,
        "rho": particle_density_kg_per_cubic_metre,
        "dist_type": "log_normal",
        "comp": "ri",
        "refractive_index": refractive_index,
        "alt_low": lower_altitude_kilometres,
        "alt_upp": lower_altitude_kilometres + thickness_kilometres,
        "radii": 128,
        "eta": 1e-6,
        "phase_quad_N": 96,
        "phase_quad_type": "L",
        "radii_quad_type": "T",
        "leg_coeffs": True,
        "leg_coeffs_type": "normalised",
        "multiprocess": False,
    }

    concentration_layer = MieLayer()
    concentration_layer.set_input_from_dict(layer_inputs)
    concentration_layer.calculate_op()
    concentration_layer.add_op_calc_output()
    concentration_layer.calc_tau()

    reference_index = int(random_generator.integers(concentration_layer.wvls.size))
    reference_wavelength_micrometres = float(concentration_layer.wvls[reference_index])
    reference_optical_depth = float(concentration_layer.tau[reference_index])
    optical_depth_inputs = dict(layer_inputs)
    optical_depth_inputs.pop("n")
    optical_depth_inputs.update(
        tau_ref=reference_optical_depth,
        tau_ref_at=reference_wavelength_micrometres,
        tau_ref_at_units="um",
    )
    optical_depth_layer = MieLayer()
    optical_depth_layer.set_input_from_dict(optical_depth_inputs)
    optical_depth_layer.calculate_op()
    optical_depth_layer.add_op_calc_output()
    optical_depth_layer.calc_tau()

    maximum_relative_errors = {}
    maximum_absolute_errors = {}
    for physical_quantity in (
        "n", "mass_loading", "s_a_den", "v_den", "beta_ext", "tau", "ssalb",
        "phase_function", "legendre_coefficient",
    ):
        original_values = np.asarray(getattr(concentration_layer, physical_quantity))
        recovered_values = np.asarray(getattr(optical_depth_layer, physical_quantity))
        np.testing.assert_allclose(
            recovered_values, original_values, rtol=1e-10, atol=1e-13,
            err_msg=f"Failed to recover {physical_quantity}",
        )
        absolute_errors = np.abs(recovered_values - original_values)
        maximum_absolute_errors[physical_quantity] = float(np.max(absolute_errors))
        if physical_quantity != "legendre_coefficient":
            maximum_relative_errors[physical_quantity] = float(np.max(
                absolute_errors / np.maximum(np.abs(original_values), np.finfo(float).tiny)
            ))

    return {
        "seed": seed,
        "number_concentration_cm_minus_3": number_concentration,
        "recovered_number_concentration_cm_minus_3": optical_depth_layer.n,
        "median_radius_um": median_radius_micrometres,
        "geometric_standard_deviation": geometric_standard_deviation,
        "particle_density_kg_m_minus_3": particle_density_kg_per_cubic_metre,
        "lower_altitude_km": lower_altitude_kilometres,
        "thickness_km": thickness_kilometres,
        "refractive_index": str(refractive_index),
        "reference_wavelength_um": reference_wavelength_micrometres,
        "reference_optical_depth": reference_optical_depth,
        "spectral_points": int(concentration_layer.wvls.size),
        "maximum_relative_errors": maximum_relative_errors,
        "maximum_absolute_errors": maximum_absolute_errors,
    }


if __name__ == "__main__":
    print(json.dumps(compare_layers(), indent=2))
