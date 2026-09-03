import numpy as np
import pytest

from srfm import optical_properties as optical
from srfm import quadrature
from srfm.size_distribution import LogNormalDistribution

pytestmark = pytest.mark.unit


def test_legendre_expansion_of_constant_phase_function():
    """Verify the Legendre expansion of a constant phase function.

    Only the zeroth coefficient should remain for isotropic scattering.
    """
    nodes, weights = quadrature.quadrature101("G", 32)
    coefficients, used = optical.legendre_polynomial_expansion(
        32, nodes, weights, np.ones(32)
    )
    assert coefficients[0] == pytest.approx(1.0, rel=2e-14)
    np.testing.assert_allclose(coefficients[1:3], 0, atol=2e-14)
    assert used >= 1


def test_normalised_and_regular_legendre_coefficients_have_expected_scaling():
    """Verify regular and normalized Legendre coefficient scaling.

    The two public conventions must differ by their documented factors.
    """
    nodes, weights = quadrature.quadrature101("G", 32)
    phase = 1 + 0.6 * nodes + 0.25 * (3 * nodes**2 - 1) / 2
    regular, _ = optical.legendre_polynomial_expansion(32, nodes, weights, phase)
    normalised, _ = optical.normalised_legendre_polynomial_expansion(
        32, nodes, weights, phase
    )
    np.testing.assert_allclose(regular[:3], [1, 0.6, 0.25], atol=2e-14)
    np.testing.assert_allclose(normalised[:3], regular[:3] / [1, 3, 5], atol=2e-14)


def test_legendre_expansion_size_limit():
    """Verify Legendre expansion rejects excessive coefficient counts.

    The angular quadrature cannot support more independent moments than nodes.
    """
    with pytest.raises(ValueError, match="Too many"):
        optical.legendre_polynomial_expansion(10001, np.ones(1), np.ones(1), np.ones(1))


def test_phase_reconstruction_from_both_coefficient_conventions():
    """Verify phase reconstruction under both coefficient conventions.

    Regular and normalized inputs should recover the same phase function.
    """
    points = np.linspace(-1, 1, 11)
    regular = optical.phase_from_legendre(3, [1, 0.6, 0.25], len(points), points)
    normalised = optical.phase_from_normalised_legendre(
        3, [1, 0.2, 0.05], len(points), points
    )
    expected = 1 + 0.6 * points + 0.25 * (3 * points**2 - 1) / 2
    np.testing.assert_allclose(regular, expected)
    np.testing.assert_allclose(normalised, expected)


def test_get_quad_user_angles_and_legendre_mode():
    """Verify user angles and Legendre mode produce valid quadrature.

    Returned angles and weights must match the requested calculation mode.
    """
    values, weights, count = optical.get_quad(angle=np.array([0, 60, 180]))
    np.testing.assert_allclose(values, [1, 0.5, -1], atol=2e-16)
    np.testing.assert_array_equal(weights, 0)
    assert count == 3
    values, weights, count = optical.get_quad(True, "L", 9, angle=np.array([45]))
    assert count == 9
    assert weights.sum() == pytest.approx(2)


def test_get_radii_produces_positive_ordered_normalised_quadrature():
    """Verify radius quadrature is positive, ordered, and normalized.

    These invariants are required before integrating particle properties.
    """
    distribution = LogNormalDistribution(n=5, r=0.4, s=1.6)
    radius, weights, total = optical.get_radii(distribution, eta=1e-6, radii=80)
    assert radius.shape == weights.shape == (80,)
    assert np.all(np.diff(radius) > 0)
    assert np.all(radius > 0)
    assert total == pytest.approx(5, rel=2e-4)


def test_user_supplied_refractive_index_scalar_and_array():
    """Verify scalar and array refractive indices are accepted.

    Scalar values should broadcast consistently over the spectral grid.
    """
    scalar = optical.get_ri(
        "ri", refractive_index=1.5 - 0.01j, wave=[1, 2], wave_size=2
    )
    array = optical.get_ri(
        "ri", refractive_index=[1.5 - 0.01j, 1.6 - 0.02j], wave=[1, 2], wave_size=2
    )
    np.testing.assert_array_equal(scalar, [1.5 - 0.01j, 1.5 - 0.01j])
    np.testing.assert_array_equal(array, [1.5 - 0.01j, 1.6 - 0.02j])


def test_user_refractive_index_validation():
    """Verify malformed user refractive indices are rejected.

    Shape and type validation prevents mismatched optical calculations.
    """
    with pytest.raises(RuntimeError, match="not defined"):
        optical.get_ri("ri", wave=[1], wave_size=1)
    with pytest.raises(ValueError, match="match wave_size"):
        optical.get_ri("ri", refractive_index=[1.5 - 0.01j], wave=[1, 2], wave_size=2)
    with pytest.raises(ValueError, match="invalid"):
        optical.get_ri("ri", refractive_index=1.5 + 0.01j, wave=[1], wave_size=1)


def test_python_mie_solver_physical_properties_and_shape():
    """Verify the Python Mie solver returns physical shaped outputs.

    Extinction, scattering, and phase arrays are checked for core invariants.
    """
    qext, qsca, phase = optical.mie_ewp(2.0, 1.5 - 0.01j, np.array([-1.0, 0.0, 1.0]))
    assert qext > 0
    assert 0 <= qsca <= qext
    assert phase.shape == (3,)
    assert np.all(np.isfinite(phase))
    assert np.all(phase >= 0)


def test_python_mie_solver_validates_angle_dimensions_and_size():
    """Verify the Python Mie solver validates angular inputs.

    Unsupported dimensionality and undersized angle grids must fail clearly.
    """
    with pytest.raises(ValueError, match="1D"):
        optical.mie_ewp(2, 1.5 - 0.01j, np.ones((2, 2)))
    with pytest.raises(ValueError, match="105000"):
        optical.mie_ewp(105001, 1.5 - 0.01j, [0])


def test_regrid_interpolates_all_optical_properties_and_preserves_shapes():
    """Verify regridding interpolates every optical property consistently.

    Scalar spectra and phase-dependent arrays must retain compatible shapes.
    """
    source = {
        "wavelengths": np.array([1.0, 2.0, 3.0]),
        "beta_ext": np.array([1.0, 2.0, 3.0]),
        "ssalb": np.array([0.2, 0.4, 0.6]),
        "phase_function": np.array([[1, 2], [2, 4], [3, 6]], dtype=float),
        "legendre_coefficient": np.array([[1, 0], [1, 0.1], [1, 0.2]], dtype=float),
    }
    result, diff = optical.regrid(source, np.array([1.5, 2.5]))
    assert diff == {}
    np.testing.assert_allclose(result["beta_ext"], [1.5, 2.5])
    np.testing.assert_allclose(result["phase_function"], [[1.5, 3], [2.5, 5]])
    assert result["legendre_coefficient"].shape == (2, 2)


def test_regrid_rejects_unknown_fields_and_nonmonotonic_grid():
    """Verify regridding rejects unknown fields and unordered grids.

    Both errors would otherwise produce ambiguous interpolation results.
    """
    with pytest.raises(RuntimeError, match="Invalid key"):
        optical.regrid({"wavelengths": [1, 2], "unknown": [1, 2]}, [1, 2])
