import numpy as np
import pytest

from srfm.size_distribution import (
    GaussianDistribution,
    LogNormalDistribution,
    SizeDistribution,
    create_distribution,
)


pytestmark = pytest.mark.unit


def test_size_distribution_is_abstract():
    with pytest.raises(TypeError):
        SizeDistribution("abstract")


def test_gaussian_construction_and_mean_are_explicitly_python_side():
    distribution = GaussianDistribution(n=10, r=2.5, s=0.3)
    assert distribution.type == "gaussian"
    assert distribution.mean() == pytest.approx(2.5)


def test_create_distribution_selects_supported_types():
    assert isinstance(create_distribution("gaussian", n=1, r=2, s=0.2), GaussianDistribution)
    assert isinstance(create_distribution("log_normal", n=1, r=2, s=1.5), LogNormalDistribution)


def test_unknown_distribution_type_rejected():
    with pytest.raises(ValueError, match="Unknown distribution type"):
        create_distribution("gamma", n=1, r=2, s=1.5)


def test_lognormal_preserves_embedded_reference_values():
    first = create_distribution("log_normal", n=1, r=np.e, s=np.e)
    second = create_distribution("log_normal", n=1, r=np.exp(2), s=np.exp(3))
    assert first.value(1) == pytest.approx(0.24197072451914337, rel=2e-15)
    assert second.value(0.5) == pytest.approx(0.17775476480655455, rel=2e-15)


def test_lognormal_mean_and_mode_properties():
    distribution = LogNormalDistribution(n=20, r=0.5, s=1.7)
    assert distribution.mean() == pytest.approx(0.5 * np.exp(0.5 * np.log(1.7) ** 2))
    radii = np.geomspace(0.01, 10, 200)
    values = distribution.value(radii)
    assert values.shape == radii.shape
    assert np.all(np.isfinite(values))
    assert np.all(values >= 0)


def test_lognormal_integrates_to_number_density():
    distribution = LogNormalDistribution(n=12.5, r=0.4, s=1.6)
    radii = np.geomspace(1e-4, 100, 30000)
    integral = np.trapezoid(distribution.value(radii), radii)
    assert integral == pytest.approx(distribution.n, rel=2e-5)


@pytest.mark.parametrize("parameterisation", ["surface", "volume"])
def test_equivalent_lognormal_parameterisations(parameterisation):
    reference = LogNormalDistribution(n=25, r=0.7, s=1.8)
    kwargs = {"r": reference.r, "s": reference.s}
    if parameterisation == "surface":
        kwargs["surface_area_density"] = reference.surface_area_density
    else:
        kwargs["volume_density"] = reference.volume_density
    equivalent = LogNormalDistribution(**kwargs)
    assert equivalent.n == pytest.approx(reference.n, rel=2e-15)
    assert equivalent.surface_area_density == pytest.approx(reference.surface_area_density)
    assert equivalent.volume_density == pytest.approx(reference.volume_density)
    np.testing.assert_allclose(equivalent.value([0.2, 0.7, 2.0]), reference.value([0.2, 0.7, 2.0]))


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"n": 1, "r": 1},
        {"n": 0, "r": 1, "s": 1.5},
        {"n": -1, "r": 1, "s": 1.5},
        {"n": 1, "r": 0, "s": 1.5},
        {"n": 1, "r": 1, "s": 0},
        {"surface_area_density": 0, "r": 1, "s": 1.5},
        {"volume_density": -1, "r": 1, "s": 1.5},
    ],
)
def test_invalid_lognormal_parameters_rejected(kwargs):
    with pytest.raises(ValueError):
        LogNormalDistribution(**kwargs)


def test_lognormal_parameter_precedence_is_deterministic():
    distribution = LogNormalDistribution(
        n=2, r=1, s=1.5, surface_area_density=999, volume_density=999
    )
    assert distribution.n == 2
    assert distribution.surface_area_density != 999
