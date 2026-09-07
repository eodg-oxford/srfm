import numpy as np
import pytest

from srfm import utilities
from srfm.forward_model import DISORT, DisortResult, SRFM

pytestmark = pytest.mark.unit


def test_disort_instances_do_not_share_mutable_dictionaries():
    """Verify default DISORT state is independent between instances."""
    first = DISORT()
    second = DISORT()
    first.disort_input["value"] = 1
    first.disort_out[1000.0] = {}
    assert second.disort_input == {}
    assert second.disort_out == {}


def test_disort_current_result_can_disable_history_retention():
    """Verify low-memory runs return one result without growing history."""
    model = DISORT(retain_history=False)
    model.disort_input = {"wvnmlo": 999.5, "wvnmhi": 1000.5}
    model.wvnm = 1000.0
    model.wvl = 10.0
    native = tuple(np.full(2, index, dtype=float) for index in range(8))
    result = model._record_result(native)
    assert isinstance(result, DisortResult)
    assert model.current_output is result
    assert model.disort_out == {}
    np.testing.assert_allclose(result.uu, native[5])


def test_srfm_selective_output_allocation_and_storage():
    """Verify only requested full-grid output arrays are allocated."""
    disort = DISORT()
    disort.disort_input = {"maxulv": 2, "maxumu": 1, "maxphi": 1}
    model = SRFM()
    model.wvnm = np.array([1000.0, 1001.0])
    model.initialize_srfm_output_arrays_from_disort(
        disort, retain_outputs={"radiance", "flup"}
    )
    assert hasattr(model, "uu")
    assert hasattr(model, "flup")
    assert not hasattr(model, "rfldir")
    result = DisortResult(
        1000.0,
        10.0,
        *(np.ones(2) for _ in range(5)),
        np.ones((1, 2, 1)),
        np.ones(1),
        np.ones(1),
    )
    model.store_disort_result(result, 0)
    np.testing.assert_allclose(model.flup[0], 1)
    np.testing.assert_allclose(model.uu[0], 1)


def test_mixed_pmom_uses_rayleigh_structure_and_reusable_workspace():
    """Verify direct mixed moments match the analytical weighted result."""
    model = DISORT()
    model.disort_input = {"maxmom": 4, "maxcly": 3}
    tau_rayleigh = np.array([0.2, 0.0, 0.3])
    tau_particle = np.array([0.1, 0.2, 0.0])
    particle_albedo = np.array([0.5, 0.8, 0.0])
    particle = np.array([1.0, 0.4, 0.3, 0.2, 0.1])
    moments = model.set_mixed_pmom(
        tau_rayleigh,
        particle_albedo,
        tau_particle,
        {0: particle, 1: particle},
        prec="double",
    )
    rayleigh = np.zeros((5, 3))
    rayleigh[0] = 1.0
    rayleigh[2] = 0.1
    denominator = tau_rayleigh + particle_albedo * tau_particle
    expected = np.zeros((5, 3))
    expected[:, 0] = (
        tau_rayleigh[0] * rayleigh[:, 0]
        + particle_albedo[0] * tau_particle[0] * particle
    ) / denominator[0]
    expected[:, 1] = particle
    expected[:, 2] = rayleigh[:, 2]
    np.testing.assert_allclose(moments, expected)
    assert moments.flags.c_contiguous


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


def test_srfm_calc_bbt_broadcasts_over_all_output_dimensions():
    """Verify brightness temperature retains every multidimensional spectrum."""
    model = SRFM()
    model.wvnm = np.array([900.0, 1000.0, 1100.0])
    model.uu = np.arange(1.0, 25.0).reshape(3, 2, 2, 2) * 1e-3

    model.calc_bbt()

    expected_wvnm = model.wvnm.reshape(3, 1, 1, 1)
    expected = utilities.convert_spectral_radiance_to_bbt(model.uu, expected_wvnm)
    assert model.bbt.shape == model.uu.shape
    np.testing.assert_allclose(model.bbt, expected)


def test_srfm_iasi_convolution_processes_every_output_spectrum(tiny_iasi_ils):
    """Verify multidimensional convolution does not omit output-level spectra."""
    model = SRFM()
    model.wvnm = np.array([999.5, 999.75, 1000.0, 1000.25, 1000.5])
    impulse = np.array([0.0, 0.0, 1.0, 0.0, 0.0]).reshape(5, 1, 1, 1)
    scales = np.arange(1.0, 13.0).reshape(1, 2, 3, 2)
    model.uu = impulse * scales

    model.convolve_with_iasi(str(tiny_iasi_ils))

    assert model.uu.shape == (5, 2, 3, 2)
    for polar in range(2):
        for level in range(3):
            for azimuthal in range(2):
                spectrum = model.uu[:, polar, level, azimuthal]
                assert spectrum.sum() == pytest.approx(scales[0, polar, level, azimuthal])
                assert spectrum[1] > 0
                assert spectrum[3] > 0


def test_srfm_interp_handles_multidimensional_radiance_and_refreshes_bbt():
    """Verify linear interpolation preserves dimensions and synchronized BBT."""
    model = SRFM()
    model.wvnm = np.array([900.0, 1000.0, 1100.0])
    slopes = np.arange(1.0, 9.0).reshape(1, 2, 2, 2)
    model.uu = model.wvnm.reshape(3, 1, 1, 1) * slopes * 1e-5
    model.calc_bbt()

    model.interp(np.array([950.0, 1050.0]))

    expected_uu = np.array([950.0, 1050.0]).reshape(2, 1, 1, 1) * slopes * 1e-5
    np.testing.assert_allclose(model.uu, expected_uu)
    expected_bbt = utilities.convert_spectral_radiance_to_bbt(
        expected_uu, model.wvnm.reshape(2, 1, 1, 1)
    )
    np.testing.assert_allclose(model.bbt, expected_bbt)
    np.testing.assert_allclose(model.wvls, 1e4 / model.wvnm)


def test_srfm_interp_rejects_extrapolation():
    """Verify interpolation cannot silently manufacture spectra outside the model grid."""
    model = SRFM()
    model.wvnm = np.array([900.0, 1000.0, 1100.0])
    model.uu = np.ones((3, 1, 1, 1))

    with pytest.raises(ValueError, match="within the source grid"):
        model.interp(np.array([850.0, 950.0]))
