import numpy as np
import pytest

from srfm.layer import GreyBodyCloud, Layer, MieLayer
from srfm.size_distribution import LogNormalDistribution


pytestmark = pytest.mark.unit


def test_layer_construction_assigns_named_parameters():
    layer = Layer("base", altitude=5, enabled=True)
    assert layer.name == "base"
    assert layer.altitude == 5
    assert layer.enabled is True


def test_mie_layer_setters_and_constructor_parameters():
    layer = MieLayer("ash", custom=3)
    layer.set_spc_lim(1000, 1002)
    layer.set_spec_units("cm-1")
    layer.set_res(1)
    layer.set_leg_coeffs_type("regular")
    assert (layer.name, layer.custom) == ("ash", 3)
    assert (layer.low_spc, layer.upp_spc, layer.res) == (1000, 1002, 1)
    assert layer.leg_coeffs_type == "regular"


def test_mie_layer_extent_from_centre_and_from_bounds():
    centred = MieLayer(center_alt=10, thick=2)
    centred.calc_layer_extent()
    assert (centred.alt_low, centred.alt_upp) == (9, 11)
    bounded = MieLayer(alt_low=4, alt_upp=8)
    bounded.calc_layer_extent()
    assert (bounded.center_alt, bounded.thick) == (6, 4)


def test_mie_layer_rejects_inconsistent_extent():
    layer = MieLayer(center_alt=10, thick=2, alt_low=7, alt_upp=11)
    with pytest.raises(AssertionError, match="do not match"):
        layer.calc_layer_extent()


def test_mie_layer_number_surface_volume_calculations_are_equivalent():
    reference = MieLayer(n=20, r=0.5, s=1.6)
    reference.n_s_v()
    by_surface = MieLayer(s_a_den=reference.s_a_den, r=0.5, s=1.6)
    by_surface.n_s_v()
    by_volume = MieLayer(v_den=reference.v_den, r=0.5, s=1.6)
    by_volume.n_s_v()
    assert by_surface.n == pytest.approx(20)
    assert by_volume.n == pytest.approx(20)


def test_mie_layer_mass_loading_and_number_concentration_round_trip():
    layer = MieLayer(mass_loading=0.2, rho=2300, thick=1, r=0.4, s=1.7, dist_type="log_normal")
    layer.nsv_or_ml()
    assert layer.n > 0
    assert layer.s_a_den > 0
    assert layer.v_den > 0
    original_mass = layer.mass_loading
    layer.mass_loading = None
    layer.nsv_or_ml()
    assert layer.mass_loading == pytest.approx(original_mass)


@pytest.mark.parametrize("missing", ["r", "s", "concentration"])
def test_mie_layer_requires_distribution_parameters(missing):
    kwargs = {"r": 0.5, "s": 1.5, "n": 2}
    if missing == "concentration":
        kwargs["n"] = None
    else:
        kwargs[missing] = None
    layer = MieLayer(**kwargs)
    with pytest.raises(RuntimeError):
        layer.nsv_or_ml()


def test_mie_layer_creates_size_distribution_and_grids():
    layer = MieLayer(n=2, r=0.5, s=1.5, dist_type="log_normal", low_spc=1000, upp_spc=1002, res=1, spec_units="cm-1")
    layer.n_s_v()
    layer.calc_size_distribution()
    layer.calc_grids()
    assert isinstance(layer.size_distribution, LogNormalDistribution)
    np.testing.assert_allclose(layer.wvnm, [1000, 1001, 1002])
    np.testing.assert_allclose(layer.wvnm * layer.wvls, 1e4)


def test_mie_layer_optical_depth_uses_kilometre_to_metre_conversion():
    layer = MieLayer(beta_ext=np.array([1e-3, 2e-3]), alt_low=4, alt_upp=5.5)
    layer.calc_tau()
    np.testing.assert_allclose(layer.tau, [1.5, 3.0])


def _grey_cloud(**overrides):
    values = dict(name="grey", low_spc=1000, upp_spc=1002, spec_units="cm-1", res=1, center_alt=5, thick=1, alt_upp=None, alt_low=None, emis=0.8, inp_tau=0.4)
    values.update(overrides)
    cloud = GreyBodyCloud()
    cloud.set_input_from_dict(values)
    return cloud


def test_grey_body_cloud_complete_python_calculation():
    cloud = _grey_cloud()
    cloud.calculate_op()
    assert (cloud.alt_low, cloud.alt_upp) == (4.5, 5.5)
    np.testing.assert_allclose(cloud.wvnm, [1000, 1001, 1002])
    np.testing.assert_allclose(cloud.tau, [0.4, 0.4, 0.4])


@pytest.mark.parametrize("emissivity", [-0.1, 1.1])
def test_grey_body_cloud_validates_emissivity(emissivity):
    cloud = _grey_cloud(emis=emissivity)
    cloud.calc_layer_extent()
    with pytest.raises(ValueError, match="Emissivity"):
        cloud.test_input_values()


def test_grey_body_regrid_preserves_constant_optical_depth():
    cloud = _grey_cloud()
    cloud.calculate_op()
    cloud.regrid(np.array([10.0, 10.005, 10.01]))
    np.testing.assert_allclose(cloud.tau, 0.4)
