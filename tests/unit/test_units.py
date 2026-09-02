import numpy as np
import pytest

from srfm import units


pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    ("decimal", "expected"),
    [(0, (0, 0, 0)), (12.5, (12, 30, 0)), (359.999, (359, 59, 56)), (-12.5, (-12, 30, 0)), (-0.5, (0, -30, 0))],
)
def test_decimal_degree_to_dms(decimal, expected):
    assert units.decimal_degree_to_DMS(decimal) == expected


@pytest.mark.parametrize("wavenumber", [1.0, 500.0, 1000.0, 50000.0])
def test_micron_wavenumber_round_trip(wavenumber):
    wavelength = units.inv_cm_to_micron(wavenumber)
    assert units.micron_to_inv_cm(wavelength) == pytest.approx(wavenumber, rel=1e-14)


@pytest.mark.parametrize("wavenumber", [1.0, 500.0, 1000.0, 50000.0])
def test_nanometre_wavenumber_round_trip(wavenumber):
    wavelength = units.inv_cm_to_nm(wavenumber)
    assert units.nm_to_inv_cm(wavelength) == pytest.approx(wavenumber, rel=1e-14)


def test_micron_and_nanometre_conversions_are_consistent():
    np.testing.assert_allclose(
        units.inv_cm_to_nm(np.array([500.0, 1000.0])),
        1000 * units.inv_cm_to_micron(np.array([500.0, 1000.0])),
        rtol=1e-15,
    )


@pytest.mark.parametrize("du", [0, 1, 300.0])
def test_dobson_units_to_column_density(du):
    assert units.DU_to_col_den(du) == pytest.approx(du * 2.69e20)


def test_dobson_units_reject_non_numeric_input():
    with pytest.raises(TypeError, match="int or float"):
        units.DU_to_col_den("300")
