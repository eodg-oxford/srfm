import numpy as np
import pytest
from scipy.special import eval_legendre

from srfm import quadrature as q


pytestmark = pytest.mark.unit


def test_bessel_zero_preserves_embedded_reference_value():
    assert q.bessel_zero(3) == pytest.approx(10.173467949597212, rel=2e-15)


@pytest.mark.parametrize(
    ("kind", "expected"),
    [("G", 0.6515842860871356), ("R", 0.6497710863213735), ("L", 0.6413165972728538)],
)
def test_first_guess_preserves_embedded_reference_values(kind, expected):
    assert q.first_guess(kind, 181, 50) == pytest.approx(expected, rel=2e-15)


@pytest.mark.parametrize(
    ("kind", "expected"),
    [("G", -0.004336062494010145), ("R", 0.0003714827832841103), ("L", -0.005288777969312314)],
)
def test_newton_correction_preserves_embedded_reference_values(kind, expected):
    assert q.newton_g(kind, 181, 0.1) == pytest.approx(expected, rel=3e-14)


def test_legendre_base_cases_and_reference_recurrence():
    assert q.legendre(0, 0.2) == (None, None, 1.0)
    assert q.legendre(1, 0.2) == (None, 1.0, 0.2)
    result = q.legendre(181, 0.1)
    expected = (0.0456111402767717, 0.0427731037482122, -0.03682815582601351)
    np.testing.assert_allclose(result, expected, rtol=3e-14, atol=1e-16)
    assert result[-1] == pytest.approx(eval_legendre(181, 0.1), rel=2e-13)


@pytest.mark.parametrize("kind", ["G", "R", "L", "T"])
def test_quadrature_nodes_weights_and_polynomial_integrals(kind):
    nodes, weights = q.quadrature101(kind, 8)
    assert nodes.shape == weights.shape == (8,)
    assert np.all(np.diff(nodes) > 0)
    assert np.all(weights > 0)
    assert np.sum(weights) == pytest.approx(2.0, rel=2e-13)
    np.testing.assert_allclose(np.sum(weights * nodes), 0.0, atol=3e-15)
    if kind in {"G", "L", "T"}:
        assert nodes[0] == pytest.approx(-1.0) if kind in {"L", "T"} else nodes[0] > -1
        assert nodes[-1] == pytest.approx(1.0) if kind in {"L", "T"} else nodes[-1] < 1


@pytest.mark.parametrize(
    ("kind", "node", "weight"),
    [
        ("G", -0.6383546791355884, 0.013323424039823995),
        ("R", -0.6431669568316702, 0.013290800444751914),
        ("L", -0.641312349855755, 0.01335472617564372),
    ],
)
def test_quadrature101_preserves_embedded_181_point_references(kind, node, weight):
    nodes, weights = q.quadrature101(kind, 181)
    assert nodes[50] == pytest.approx(node, rel=2e-14)
    assert weights[50] == pytest.approx(weight, rel=3e-14)


def test_shift_quadrature_preserves_integrals_and_reference():
    nodes, weights = q.quadrature101("L", 181)
    shifted_nodes, shifted_weights = q.shift_quadrature(nodes, weights, 0, 180)
    assert shifted_nodes[50] == pytest.approx(32.28188851298205, rel=2e-14)
    assert shifted_weights[50] == pytest.approx(1.2019253558079348, rel=3e-14)
    assert shifted_nodes[[0, -1]] == pytest.approx([0, 180])
    assert shifted_weights.sum() == pytest.approx(180.0, rel=2e-13)


def test_quadrature_integrates_on_transformed_interval():
    nodes, weights = q.quadrature("G", 6, 2.0, 5.0)
    assert np.sum(weights * nodes**3) == pytest.approx((5**4 - 2**4) / 4, rel=2e-14)


@pytest.mark.parametrize("function", [q.first_guess, q.newton_g])
def test_invalid_quadrature_type_rejected_by_helpers(function):
    args = ("bad", 8, 1 if function is q.first_guess else 0.1)
    with pytest.raises(ValueError, match="Invalid quadrature type"):
        function(*args)


@pytest.mark.parametrize("kind", ["S", "bad"])
def test_unsupported_quadrature_type_rejected(kind):
    with pytest.raises(ValueError):
        q.quadrature101(kind, 4)


@pytest.mark.parametrize("kind,npts", [("T", 1), ("L", 1), ("G", 0)])
def test_invalid_quadrature_sizes_rejected(kind, npts):
    with pytest.raises(ValueError):
        q.quadrature101(kind, npts)


def test_non_integer_quadrature_size_rejected():
    with pytest.raises(TypeError):
        q.quadrature101("G", 4.5)
