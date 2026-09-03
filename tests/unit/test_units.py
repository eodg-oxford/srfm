import numpy as np
import pytest

from srfm import units

pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    ("decimal", "expected"),
    [
        (0, (0, 0, 0)),
        (12.5, (12, 30, 0)),
        (359.999, (359, 59, 56)),
        (-12.5, (-12, 30, 0)),
        (-0.5, (0, -30, 0)),
    ],
)
def test_decimal_degree_to_dms(decimal, expected):
    """Verify decimal degrees convert to the expected DMS tuple.

    Positive and negative inputs preserve their sign and fractional components.

    Args:
        decimal: Parameterized decimal-degree value.
        expected: Expected degree, minute, and second components.
    """
    assert units.decimal_degree_to_DMS(decimal) == expected


@pytest.mark.parametrize("wavenumber", [1.0, 500.0, 1000.0, 50000.0])
def test_micron_wavenumber_round_trip(wavenumber):
    """Verify micron conversion round-trips representative wavenumbers.

    Reciprocal conversion should recover the original spectral coordinate.

    Args:
        wavenumber: Parameterized wavenumber in inverse centimetres.
    """
    wavelength = units.inv_cm_to_micron(wavenumber)
    assert units.micron_to_inv_cm(wavelength) == pytest.approx(wavenumber, rel=1e-14)


@pytest.mark.parametrize("wavenumber", [1.0, 500.0, 1000.0, 50000.0])
def test_nanometre_wavenumber_round_trip(wavenumber):
    """Verify nanometre conversion round-trips representative wavenumbers.

    Reciprocal conversion should recover the original spectral coordinate.

    Args:
        wavenumber: Parameterized wavenumber in inverse centimetres.
    """
    wavelength = units.inv_cm_to_nm(wavenumber)
    assert units.nm_to_inv_cm(wavelength) == pytest.approx(wavenumber, rel=1e-14)


def test_micron_and_nanometre_conversions_are_consistent():
    """Verify micron and nanometre conversions share one physical scale.

    Equivalent wavelengths should differ by exactly the expected factor.
    """
    np.testing.assert_allclose(
        units.inv_cm_to_nm(np.array([500.0, 1000.0])),
        1000 * units.inv_cm_to_micron(np.array([500.0, 1000.0])),
        rtol=1e-15,
    )


@pytest.mark.parametrize("du", [0, 1, 300.0])
def test_dobson_units_to_column_density(du):
    """Verify Dobson units convert to the expected column density.

    Representative magnitudes anchor the physical conversion constant.

    Args:
        du: Parameterized Dobson-unit quantity.
    """
    assert units.DU_to_col_den(du) == pytest.approx(du * 2.69e20)


def test_dobson_units_reject_non_numeric_input():
    """Verify Dobson conversion rejects non-numeric values.

    A type error avoids silently producing meaningless column densities.
    """
    with pytest.raises(TypeError, match="int or float"):
        units.DU_to_col_den("300")
