import numpy as np
import pytest

from srfm.forward_model import SRFM

pytestmark = pytest.mark.unit


def test_srfm_iasi_convolution_preserves_shape_and_spreads_impulse(tiny_iasi_ils):
    """Verify IASI convolution preserves shape and spreads an impulse.

    The symmetric synthetic ILS should conserve the array dimensions while
    distributing radiance into neighboring samples.

    Args:
        tiny_iasi_ils: Synthetic instrument-line-shape fixture.
    """
    model = SRFM()
    model.wvnm = np.array([999.5, 999.75, 1000.0, 1000.25, 1000.5])
    model.uu = np.array([0.0, 0.0, 1.0, 0.0, 0.0]).reshape(5, 1, 1, 1)

    model.convolve_with_iasi(str(tiny_iasi_ils))

    assert model.uu.shape == (5, 1, 1, 1)
    assert np.all(model.uu >= 0)
    assert model.uu[2, 0, 0, 0] < 1.0
    assert model.uu[1, 0, 0, 0] > 0
    assert model.uu[3, 0, 0, 0] > 0
    assert model.uu[:, 0, 0, 0].sum() == pytest.approx(1.0)


def test_srfm_iasi_convolution_rejects_irregular_grid(tiny_iasi_ils):
    """Verify IASI convolution rejects an irregular spectral grid.

    The convolution implementation derives one spacing and therefore requires
    regularly sampled wavenumbers.

    Args:
        tiny_iasi_ils: Synthetic instrument-line-shape fixture.
    """
    model = SRFM()
    model.wvnm = np.array([999.5, 999.75, 1000.1, 1000.5])
    model.uu = np.ones((4, 1, 1, 1))

    with pytest.raises(ValueError, match="not regular"):
        model.convolve_with_iasi(str(tiny_iasi_ils))
