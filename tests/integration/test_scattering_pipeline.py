import numpy as np
import pytest

import srfm
from srfm.inputs import Inputs
from srfm.layer import MieLayer
from srfm.size_distribution import LogNormalDistribution

pytestmark = pytest.mark.integration


def test_size_distribution_quadrature_and_native_optical_properties(
    require_native, run_native_case
):
    """Verify the distribution-to-native-Mie calculation pipeline.

    The test checks physical output invariants after Python quadrature feeds the
    compiled optical-properties calculation.

    Args:
        require_native: Fixture helper for extension availability.
        run_native_case: Fixture helper that isolates native execution.
    """
    require_native(srfm.mie_module, "Mie")
    result, _ = run_native_case("mie", {})
    assert result["beta_ext"].shape == (1,)
    assert result["phase_function"].shape == (1, 3)
    assert result["beta_ext"][0] > 0
    assert 0 <= result["ssalb"][0] <= 1
    assert np.all(np.isfinite(result["phase_function"]))
    assert np.all(result["phase_function"] >= 0)


def test_synthetic_aria_reader_feeds_native_optical_calculation(
    tiny_ri_file, require_native, run_native_case
):
    """Verify synthetic ARIA indices feed the compiled Mie calculation.

    This covers the boundary from text refractive-index data to native optical
    properties without an external data catalogue.

    Args:
        tiny_ri_file: Synthetic refractive-index fixture.
        require_native: Fixture helper for extension availability.
        run_native_case: Fixture helper that isolates native execution.
    """
    require_native(srfm.mie_module, "Mie")
    result, _ = run_native_case(
        "mie", {"ri_file": str(tiny_ri_file), "radius": 0.2, "spread": 1.4}
    )
    assert result["beta_ext"][0] > 0
    assert np.isfinite(result["ssalb"][0])


def test_inputs_to_mie_layer_prepares_consistent_distribution_and_grids(tmp_path):
    """Verify input mappings prepare a consistent Python-side Mie layer.

    Distribution parameters and spectral grids are checked before any compiled
    calculation is invoked.

    Args:
        tmp_path: Pytest temporary-path fixture.
    """
    driver = tmp_path / "driver.py"
    driver.write_text(
        "inputs = {'layer': {"
        "'name':'test', 'low_spc':1000, 'upp_spc':1002, 'res':1, 'spec_units':'cm-1',"
        "'mass_loading':0.1, 'n':None, 'r':0.3, 's':1.5, 'rho':2000,"
        "'s_a_den':None, 'v_den':None, 'dist_type':'log_normal', 'comp':'ri',"
        "'center_alt':5, 'thick':1, 'alt_upp':None, 'alt_low':None,"
        "'radii':10, 'eta':1e-5, 'phase_quad_N':12, 'phase_quad_type':'L',"
        "'radii_quad_type':'T', 'leg_coeffs':True, 'leg_coeffs_type':'normalised',"
        "'multiprocess':False, 'refractive_index':(1.5-0.01j)}}\n",
        encoding="utf-8",
    )
    inputs = Inputs()
    inputs.read_srfm_drv(driver, validate=False)
    layer = MieLayer()
    layer.set_input_from_dict(inputs.values["layer"])
    layer.calc_layer_extent()
    layer.nsv_or_ml()
    layer.calc_size_distribution()
    layer.calc_grids()
    assert layer.n > 0
    assert isinstance(layer.size_distribution, LogNormalDistribution)
    np.testing.assert_allclose(layer.wvnm, [1000, 1001, 1002])
    assert (layer.alt_low, layer.alt_upp) == (4.5, 5.5)


@pytest.mark.parametrize("particle_case", ["lognormal", "gamma_array", "aria"])
def test_native_mie_optical_depth_recovers_random_layer(
    require_native, run_native_case, particle_case
):
    """Reference optical depth recovers concentration and every native spectrum.

    Args:
        require_native: Fixture checking compiled Mie availability.
        run_native_case: Fixture isolating native calculations in a subprocess.
        particle_case: Distribution and refractive-index representation to exercise.
    """
    require_native(srfm.mie_module, "Mie")
    random_generator = np.random.default_rng(20261012)
    layer_inputs = {
        "name": "random_round_trip",
        "low_spc": 8.0,
        "upp_spc": 12.0,
        "spec_units": "um",
        "res": 0.5,
        "n": float(random_generator.uniform(5.0, 250.0)),
        "r": float(random_generator.uniform(0.2, 1.8)),
        "s": float(random_generator.uniform(1.2, 1.9)),
        "rho": float(random_generator.uniform(1000.0, 2500.0)),
        "alt_low": 3.0,
        "alt_upp": float(np.round(random_generator.uniform(3.2, 5.0), 3)),
        "dist_type": "log_normal",
        "comp": "ri",
        "refractive_index": complex(
            random_generator.uniform(1.3, 1.7), -random_generator.uniform(0.001, 0.05)
        ),
        "radii": 64,
        "eta": 1e-6,
        "radii_quad_type": "T",
        "phase_quad_N": 48,
        "phase_quad_type": "L",
        "leg_coeffs": True,
        "leg_coeffs_type": "normalised",
        "multiprocess": False,
    }
    if particle_case == "gamma_array":
        layer_inputs.update(
            dist_type="gamma", s=0.1,
            refractive_index=np.linspace(1.4 - 0.01j, 1.6 - 0.03j, 9),
        )
    elif particle_case == "aria":
        layer_inputs.update(comp="sulphuric acid", refractive_index=None)
    references = [
        {"units": "um"}, {"units": "nm"}, {"units": "cm-1"},
        {"units": "um", "include_density": False},
    ]
    if particle_case == "lognormal":
        references.append({"units": "cm-1", "multiprocess": True})
    result, _ = run_native_case("mie_optical_depth_round_trip", {
        "inputs": layer_inputs, "reference_index": 3, "references": references,
    })

    baseline = result["baseline"]
    assert baseline["n"] == layer_inputs["n"]
    assert np.all(baseline["tau"] > 0)
    assert baseline["mass_loading"] > 0
    for configuration, reconstructed in zip(references, result["reconstructed"]):
        for attribute in (
            "n", "s_a_den", "v_den", "beta_ext", "tau", "ssalb",
            "phase_function", "legendre_coefficient",
        ):
            np.testing.assert_allclose(
                reconstructed[attribute], baseline[attribute], rtol=1e-10, atol=1e-13,
                err_msg=f"{particle_case}: {attribute}, {configuration}",
            )
        if configuration.get("include_density", True):
            np.testing.assert_allclose(
                reconstructed["mass_loading"], baseline["mass_loading"], rtol=1e-10
            )
        else:
            assert reconstructed["mass_loading"] is None
