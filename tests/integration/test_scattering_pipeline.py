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
    require_native(srfm.mie_module, "Mie")
    result, _ = run_native_case(
        "mie", {"ri_file": str(tiny_ri_file), "radius": 0.2, "spread": 1.4}
    )
    assert result["beta_ext"][0] > 0
    assert np.isfinite(result["ssalb"][0])


def test_inputs_to_mie_layer_prepares_consistent_distribution_and_grids(tmp_path):
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
