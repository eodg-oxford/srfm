"""Independent numerical checks for analytic particle-size distributions."""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import (
    burr12,
    gamma as scipy_gamma,
    gengamma,
    invgamma,
    lognorm,
    norm,
    truncnorm,
)

from srfm.size_distribution import (
    GaussianDistribution,
    GammaDistribution,
    InverseModifiedGammaDistribution,
    LogNormalDistribution,
    ModifiedGammaDistribution,
    MultimodeLogNormalDistribution,
    RegularisedPowerLawDistribution,
)

pytestmark = pytest.mark.unit

ORACLE_RTOL = 2e-11


def _assert_matches_scipy(distribution, oracle):
    """Compare density, moments, center metrics, area, and volume to SciPy.

    Args:
        distribution: SRFM size-distribution instance.
        oracle: Frozen ``scipy.stats`` continuous distribution normalized to one.
    """
    radii = oracle.ppf([0.1, 0.25, 0.5, 0.75, 0.9])
    np.testing.assert_allclose(
        distribution.value(radii),
        distribution.n * oracle.pdf(radii),
        rtol=ORACLE_RTOL,
        atol=1e-13,
    )
    for order in range(5):
        assert distribution.moment(order) == pytest.approx(
            distribution.n * oracle.moment(order), rel=ORACLE_RTOL
        )
    assert distribution.mean() == pytest.approx(oracle.mean(), rel=ORACLE_RTOL)
    assert distribution.median_radius == pytest.approx(oracle.median(), rel=ORACLE_RTOL)
    assert distribution.surface_area_density == pytest.approx(
        4.0 * np.pi * distribution.n * oracle.moment(2), rel=ORACLE_RTOL
    )
    assert distribution.volume_density == pytest.approx(
        4.0 * np.pi * distribution.n * oracle.moment(3) / 3.0,
        rel=ORACLE_RTOL,
    )


def test_gaussian_matches_scipy_normal_and_frozen_values():
    """Check Grainger Figure-1 Gaussian parameters against SciPy and constants."""
    distribution = GaussianDistribution(n=1.0, r=5.0, s=1.5)
    _assert_matches_scipy(distribution, norm(loc=5.0, scale=1.5))

    assert distribution.value(distribution.median_radius) == pytest.approx(
        0.265961520267622
    )
    assert distribution.mean() == pytest.approx(5.0)
    assert distribution.moment(2) == pytest.approx(27.25)
    assert distribution.moment(3) == pytest.approx(158.75)


def test_positive_truncated_gaussian_matches_scipy_truncnorm():
    """Check the positive-radius Gaussian branch and its conditional moments."""
    mean, standard_deviation = 0.5, 1.0
    distribution = GaussianDistribution(
        n=2.0, r=mean, s=standard_deviation, truncate=True
    )
    oracle = truncnorm(
        a=-mean / standard_deviation,
        b=np.inf,
        loc=mean,
        scale=standard_deviation,
    )
    _assert_matches_scipy(distribution, oracle)

    assert distribution.mean() == pytest.approx(1.00916043383703)
    assert distribution.median_radius == pytest.approx(0.896871175089544)
    assert distribution.moment(2) == pytest.approx(3.00916043383703)
    assert distribution.moment(3) == pytest.approx(5.54122195226665)


def test_lognormal_matches_scipy_lognorm_and_frozen_values():
    """Check SRFM geometric-standard-deviation conversion against SciPy."""
    number, median, geometric_standard_deviation = 4.0, 0.3, 1.7
    distribution = LogNormalDistribution(
        n=number, r=median, s=geometric_standard_deviation
    )
    oracle = lognorm(s=np.log(geometric_standard_deviation), scale=median)
    _assert_matches_scipy(distribution, oracle)

    assert distribution.value(median) == pytest.approx(10.0244010655385)
    assert distribution.mean() == pytest.approx(0.345352503635596)
    assert distribution.moment(2) == pytest.approx(0.632219543702498)
    assert distribution.moment(3) == pytest.approx(0.383438698601742)


def test_multimode_lognormal_matches_sum_of_scipy_modes():
    """Check multimode density, moments, and median against independent modes."""
    mode_number = np.array([12.0, 3.5, 0.7])
    mode_radius = np.array([0.04, 0.3, 2.0])
    mode_spread = np.array([1.45, 1.7, 2.1])
    distribution = MultimodeLogNormalDistribution(
        n=mode_number, r=mode_radius, s=mode_spread
    )
    oracles = [
        lognorm(s=np.log(spread), scale=radius)
        for radius, spread in zip(mode_radius, mode_spread)
    ]

    radii = np.geomspace(0.01, 10.0, 21)
    expected_density = sum(
        number * oracle.pdf(radii) for number, oracle in zip(mode_number, oracles)
    )
    np.testing.assert_allclose(
        distribution.value(radii), expected_density, rtol=ORACLE_RTOL
    )
    for order in range(5):
        expected_moment = sum(
            number * oracle.moment(order)
            for number, oracle in zip(mode_number, oracles)
        )
        assert distribution.moment(order) == pytest.approx(
            expected_moment, rel=ORACLE_RTOL
        )

    def normalized_mixture_cdf(radius):
        return (
            sum(
                number * oracle.cdf(radius)
                for number, oracle in zip(mode_number, oracles)
            )
            / mode_number.sum()
        )

    expected_median = brentq(
        lambda radius: normalized_mixture_cdf(radius) - 0.5, 1e-8, 1e3
    )
    assert distribution.median_radius == pytest.approx(expected_median, rel=ORACLE_RTOL)
    assert distribution.mean() == pytest.approx(0.22016145155308)
    assert distribution.median_radius == pytest.approx(0.0473423543889891)
    assert distribution.moment(2) == pytest.approx(8.99809057415678)
    assert distribution.moment(3) == pytest.approx(67.0156811311721)
    assert distribution.cdf(distribution.median_radius) == pytest.approx(
        distribution.n / 2.0
    )


def test_gamma_matches_scipy_gamma_and_reference_parameters():
    """Check alpha=2, b=0.6 against SciPy and Equation 140 moments.

    The standard deviation is derived from the Equation 140 raw moments as
    ``sqrt(M2/M0 - (M1/M0)**2)``. Equation 144 supplies only the mean.
    """
    distribution = GammaDistribution(n=1.0, effective_radius=25.0 / 3.0, s=0.2)
    oracle = scipy_gamma(a=3.0, scale=1.0 / 0.6)
    _assert_matches_scipy(distribution, oracle)

    assert distribution.alpha == pytest.approx(2.0)
    assert distribution.b == pytest.approx(0.6)
    assert distribution.mean() == pytest.approx(5.0)
    assert distribution.median_radius == pytest.approx(4.45676718953927)
    assert distribution.moment(2) == pytest.approx(33.3333333333334)
    assert distribution.moment(3) == pytest.approx(277.777777777778)
    standard_deviation = np.sqrt(
        distribution.moment(2) / distribution.n - distribution.mean() ** 2
    )
    assert standard_deviation == pytest.approx(np.sqrt(3.0) / 0.6)


def test_modified_gamma_matches_scipy_gengamma_and_reference_parameters():
    """Check Grainger Figure-3 parameters against SciPy generalized gamma."""
    alpha, b, cutoff_exponent = 2.0, 8.7e-6, 6.19
    shape = (alpha + 1.0) / cutoff_exponent
    scale = b ** (-1.0 / cutoff_exponent)
    oracle = gengamma(a=shape, c=cutoff_exponent, scale=scale)
    second, third, fourth = (oracle.moment(order) for order in (2, 3, 4))
    effective_radius = third / second
    effective_variance = second * fourth / third**2 - 1.0
    distribution = ModifiedGammaDistribution(
        n=1.0,
        effective_radius=effective_radius,
        s=effective_variance,
        gamma=cutoff_exponent,
    )
    _assert_matches_scipy(distribution, oracle)

    assert distribution.alpha == pytest.approx(alpha, rel=2e-10)
    assert distribution.b == pytest.approx(b, rel=2e-10)
    assert distribution.mean() == pytest.approx(5.0033507369075)
    assert distribution.median_radius == pytest.approx(5.12238126720944)
    assert distribution.moment(2) == pytest.approx(27.2866271796096)
    assert distribution.moment(3) == pytest.approx(158.012918338124)


@pytest.mark.parametrize(
    "distribution_class, extra_kwargs",
    [
        (GammaDistribution, {}),
        (ModifiedGammaDistribution, {"gamma": 2.5}),
    ],
    ids=["gamma", "modified_gamma"],
)
def test_gamma_radius_parameterisations_are_equivalent(
    distribution_class, extra_kwargs
):
    """Check that number-median and effective-radius inputs are equivalent."""
    from_effective_radius = distribution_class(
        n=7.0,
        effective_radius=2.3,
        s=0.08,
        **extra_kwargs,
    )
    from_median_radius = distribution_class(
        n=7.0,
        r=from_effective_radius.median_radius,
        s=0.08,
        **extra_kwargs,
    )

    assert from_median_radius.effective_radius == pytest.approx(
        from_effective_radius.effective_radius, rel=ORACLE_RTOL
    )
    assert from_median_radius.alpha == pytest.approx(
        from_effective_radius.alpha, rel=ORACLE_RTOL
    )
    assert from_median_radius.b == pytest.approx(
        from_effective_radius.b, rel=ORACLE_RTOL
    )
    radii = np.geomspace(1e-3, 20.0, 100)
    np.testing.assert_allclose(
        from_median_radius.value(radii),
        from_effective_radius.value(radii),
        rtol=ORACLE_RTOL,
    )


def test_inverse_modified_gamma_special_case_matches_scipy_invgamma():
    """Check gamma=1 inverse modified gamma against SciPy inverse gamma."""
    distribution = InverseModifiedGammaDistribution(n=3.0, alpha=6.0, b=2.0, gamma=1.0)
    oracle = invgamma(a=5.0, scale=2.0)
    _assert_matches_scipy(distribution, oracle)

    assert distribution.mean() == pytest.approx(0.5)
    assert distribution.median_radius == pytest.approx(0.428182191129108)
    assert distribution.moment(2) == pytest.approx(1.0)
    assert distribution.moment(3) == pytest.approx(1.0)


def test_regularised_power_law_matches_scipy_burr12():
    """Check the regularised power law against its Burr-XII equivalent."""
    distribution = RegularisedPowerLawDistribution(n=5.0, alpha=3.0, b=0.8, gamma=4.0)
    oracle = burr12(c=3.0, d=3.0, scale=0.8)
    _assert_matches_scipy(distribution, oracle)

    assert distribution.mean() == pytest.approx(0.537422033847176)
    assert distribution.median_radius == pytest.approx(0.510548656688515)
    assert distribution.moment(2) == pytest.approx(1.71975050831096)
    assert distribution.moment(3) == pytest.approx(1.28)


@pytest.mark.parametrize(
    "distribution",
    [
        GaussianDistribution(n=1.0, r=5.0, s=1.5, truncate=True),
        LogNormalDistribution(n=2.0, r=0.3, s=1.7),
        MultimodeLogNormalDistribution(n=[2.0, 0.5], r=[0.1, 1.0], s=[1.5, 2.0]),
        GammaDistribution(n=1.0, effective_radius=2.0, s=0.1),
        ModifiedGammaDistribution(n=1.0, effective_radius=2.0, s=0.1, gamma=2.0),
        InverseModifiedGammaDistribution(n=1.0, alpha=8.0, b=0.01, gamma=2.0),
        RegularisedPowerLawDistribution(n=1.0, alpha=3.0, b=1.0, gamma=4.0),
    ],
    ids=lambda distribution: distribution.type,
)
def test_positive_radius_distributions_integrate_to_number_density(distribution):
    """Numerically integrate each positive-radius PDF in log-radius space."""

    def log_radius_integrand(log_radius):
        radius = np.exp(log_radius)
        return distribution.value(radius) * radius

    integral = quad(
        log_radius_integrand,
        -40.0,
        40.0,
        epsabs=1e-10,
        epsrel=1e-10,
        limit=300,
    )[0]
    assert integral == pytest.approx(distribution.n, rel=2e-9, abs=1e-10)


@pytest.mark.parametrize(
    "distribution, constructor_kwargs",
    [
        (
            GammaDistribution(n=10, effective_radius=2.0, s=0.1),
            {"effective_radius": 2.0, "s": 0.1},
        ),
        (
            ModifiedGammaDistribution(n=10, effective_radius=2.0, s=0.1, gamma=2.0),
            {"effective_radius": 2.0, "s": 0.1, "gamma": 2.0},
        ),
        (
            InverseModifiedGammaDistribution(n=10, alpha=8, b=0.01, gamma=2),
            {"alpha": 8, "b": 0.01, "gamma": 2},
        ),
        (
            RegularisedPowerLawDistribution(n=10, alpha=3, b=1, gamma=4),
            {"alpha": 3, "b": 1, "gamma": 4},
        ),
    ],
    ids=["gamma", "modified_gamma", "inverse_modified_gamma", "power_law"],
)
def test_equivalent_density_parameterisations(distribution, constructor_kwargs):
    """Check number, surface, and volume inputs reconstruct the same PDF."""
    distribution_class = type(distribution)
    from_surface = distribution_class(
        surface_area_density=distribution.surface_area_density,
        **constructor_kwargs,
    )
    from_volume = distribution_class(
        volume_density=distribution.volume_density,
        **constructor_kwargs,
    )
    radii = np.geomspace(1e-3, 10.0, 50)
    assert from_surface.n == pytest.approx(distribution.n)
    assert from_volume.n == pytest.approx(distribution.n)
    np.testing.assert_allclose(from_surface.value(radii), distribution.value(radii))
    np.testing.assert_allclose(from_volume.value(radii), distribution.value(radii))


def test_heavy_tail_distributions_report_divergent_moments():
    """Ensure non-existent heavy-tail moments are explicit rather than truncated."""
    inverse = InverseModifiedGammaDistribution(n=1, alpha=3, b=1, gamma=1)
    power_law = RegularisedPowerLawDistribution(n=1, alpha=2, b=1, gamma=2)

    assert np.isinf(inverse.moment(2))
    assert np.isinf(inverse.surface_area_density)
    assert np.isinf(inverse.volume_density)
    assert np.isinf(power_law.moment(2))
    assert np.isinf(power_law.surface_area_density)
    assert np.isinf(power_law.volume_density)
