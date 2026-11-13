"""Tests for independent, in-place satellite-instrument processing."""

from inspect import signature

import numpy as np
import pytest

from srfm.forward_model import SRFM
from srfm.iasi_main import run_srfm as run_iasi_srfm
from srfm.instruments import (
    IASI,
    apply_banded_factor,
    create_instrument,
    create_instrument_from_inputs,
)
from srfm.instruments.base import ConvolvingInstrument, NoiseInstrument
from srfm.main import (
    _resolve_convolving_instrument,
    _resolve_final_grid,
    _resolve_noise_instrument,
    run_srfm,
)
from srfm.oxharp_main import run_srfm as run_oxharp_srfm

pytestmark = pytest.mark.unit


def test_runner_api_accepts_instruments_only_through_inputs():
    """None of the three runners exposes a separate instrument argument."""
    for runner in (run_srfm, run_oxharp_srfm, run_iasi_srfm):
        assert tuple(signature(runner).parameters) == ("inp",)


def test_iasi_convolution_mutates_and_returns_the_same_result(tiny_iasi_ils):
    """Convolution uses bounded workspace rather than copying the complete result."""
    result = SRFM()
    result.wvnm = np.array([999.5, 999.75, 1000.0, 1000.25, 1000.5])
    impulse = np.array([0.0, 0.0, 1.0, 0.0, 0.0]).reshape(5, 1, 1, 1)
    result.uu = impulse * np.arange(1.0, 13.0).reshape(1, 2, 3, 2)
    radiance_identity = id(result.uu)

    returned = IASI(tiny_iasi_ils).convolve(result)

    assert returned is result
    assert id(result.uu) == radiance_identity
    assert result.uu.shape == (5, 2, 3, 2)
    assert result.instrument_name == "iasi"
    assert result.noise_applied is False
    assert result.processing_history[-1]["operation"] == "convolution"
    for polar in range(2):
        for level in range(3):
            for azimuthal in range(2):
                scale = 1 + azimuthal + 2 * level + 6 * polar
                assert result.uu[:, polar, level, azimuthal].sum() == pytest.approx(
                    scale
                )


def test_iasi_convolution_refreshes_existing_brightness_temperature(tiny_iasi_ils):
    """Derived BBT remains synchronized without retaining its old full array."""
    result = SRFM()
    result.wvnm = np.array([999.5, 999.75, 1000.0, 1000.25, 1000.5])
    result.uu = np.array([1.0, 1.0, 2.0, 1.0, 1.0]).reshape(5, 1, 1, 1)
    result.calc_bbt()
    old_bbt = result.bbt

    IASI(tiny_iasi_ils).convolve(result)

    assert result.bbt is not old_bbt
    assert result.bbt.shape == result.uu.shape
    assert np.all(np.isfinite(result.bbt))


def test_apply_banded_factor_multiplies_compact_lower_diagonals():
    """Every compact row contributes to its corresponding lower diagonal."""
    factor_banded = np.array(
        [
            [1.0, 2.0, 3.0, 4.0],
            [10.0, 20.0, 30.0, 0.0],
            [100.0, 200.0, 0.0, 0.0],
        ]
    )

    correlated_noise = apply_banded_factor(factor_banded, np.ones(4))

    np.testing.assert_allclose(correlated_noise, [1.0, 12.0, 123.0, 234.0])


def test_apply_banded_factor_rejects_incompatible_shapes():
    """Malformed factors and independent samples fail with useful diagnostics."""
    with pytest.raises(ValueError, match="two-dimensional"):
        apply_banded_factor(np.ones(4), np.ones(4))
    with pytest.raises(ValueError, match=r"shape \(4,\)"):
        apply_banded_factor(np.ones((2, 4)), np.ones(3))


def _iasi_result(spectrum_count=1):
    """Create a small-geometry result on the complete IASI channel grid."""
    result = SRFM()
    result.wvnm = 645.0 + np.arange(8461) * 0.25
    result.uu = np.full((8461, spectrum_count, 1, 1), 0.1)
    result.retained_outputs = frozenset({"uu"})
    result.spectral_representation = "iasi_ils_convolved_spectral_radiance"
    return result


def test_iasi_noise_is_reproducible_and_updates_every_spectrum():
    """A supplied generator reproduces independent multidimensional spectra."""
    first = _iasi_result(spectrum_count=2)
    second = _iasi_result(spectrum_count=2)
    first.calc_bbt()
    original_brightness_temperature = first.bbt.copy()
    instrument = IASI(iasi_instrument="IASI-C", iasi_ifov=4)

    returned = instrument.add_noise(first, rng=np.random.default_rng(731))
    instrument.add_noise(second, rng=np.random.default_rng(731))

    assert returned is first
    assert first.noise_applied is True
    assert first.iasi_instrument == "M03"
    assert first.iasi_ifov == 4
    assert first.processing_history[-1]["noise_factor"] == "M03_4"
    assert not np.array_equal(first.uu[:, 0, 0, 0], first.uu[:, 1, 0, 0])
    np.testing.assert_array_equal(first.uu, second.uu)
    assert not np.array_equal(first.bbt, original_brightness_temperature)


def test_iasi_noise_accepts_partial_grid_and_preserves_full_factor_statistics():
    """A partial result receives the same selected realization as a full result."""
    instrument = IASI(iasi_instrument="M01", iasi_ifov=1)
    complete = _iasi_result()
    partial = _iasi_result()
    partial.wvnm = partial.wvnm[1228:1293]
    partial.uu = partial.uu[1228:1293].copy()

    instrument.add_noise(complete, rng=np.random.default_rng(1))
    instrument.add_noise(partial, rng=np.random.default_rng(1))

    np.testing.assert_array_equal(partial.uu, complete.uu[1228:1293])


def test_iasi_noise_rejects_unconvolved_or_misaligned_spectrum():
    """Noise requires prior IASI convolution and exact contiguous channels."""
    instrument = IASI(iasi_instrument="M01", iasi_ifov=1)
    unconvolved = _iasi_result()
    unconvolved.spectral_representation = "monochromatic_spectral_radiance"
    with pytest.raises(ValueError, match="set convolve='iasi'"):
        instrument.add_noise(unconvolved, rng=np.random.default_rng(1))

    misaligned = _iasi_result()
    misaligned.wvnm = misaligned.wvnm + 0.25
    with pytest.raises(ValueError, match="contiguous subset"):
        instrument.add_noise(misaligned, rng=np.random.default_rng(1))


def test_iasi_output_grid_uses_channels_inside_two_sided_margin():
    """Partial IASI convolution derives its channels from the calculation grid."""
    values = {
        "spc_wvnmlo": 950.0,
        "spc_wvnmhi": 970.0,
        "spc_units": "cm-1",
        "fin_wvnmlo": 950.0,
        "fin_wvnmhi": 970.0,
        "fin_res": 0.1,
    }

    with pytest.warns(UserWarning, match="952 to 968 cm-1"):
        output_grid = IASI().output_grid(values)

    assert output_grid[0] == 952.0
    assert output_grid[-1] == 968.0
    assert output_grid.size == 65
    np.testing.assert_array_equal(np.diff(output_grid), np.full(64, 0.25))


def test_iasi_output_grid_can_cover_every_instrument_channel():
    """The complete calculation range returns all 8461 IASI channels."""
    output_grid = IASI().output_grid(
        {
            "spc_wvnmlo": 643.0,
            "spc_wvnmhi": 2762.0,
            "spc_units": "cm-1",
        }
    )

    assert output_grid.size == 8461
    assert output_grid[0] == 645.0
    assert output_grid[-1] == 2760.0


def test_ordinary_final_grid_is_unchanged_without_instrument_convolution():
    """Non-IASI calculations continue to use every configured final-grid field."""
    output_grid = _resolve_final_grid(
        {
            "fin_wvnmlo": 800.0,
            "fin_wvnmhi": 801.0,
            "fin_res": 0.25,
        }
    )

    np.testing.assert_array_equal(
        output_grid, np.array([800.0, 800.25, 800.5, 800.75, 801.0])
    )


@pytest.mark.parametrize(
    ("lower", "upper", "message"),
    [
        (950.0, 953.75, "at least 4 cm-1"),
        (642.75, 970.0, "within 643 to 2762 cm-1"),
        (950.0, 2762.25, "within 643 to 2762 cm-1"),
    ],
)
def test_iasi_output_grid_rejects_unsupported_calculation_ranges(
    lower, upper, message
):
    """IASI output rejects missing margins and excess spectral coverage."""
    with pytest.raises(ValueError, match=message):
        IASI().output_grid(
            {
                "spc_wvnmlo": lower,
                "spc_wvnmhi": upper,
                "spc_units": "cm-1",
            }
        )


def test_iasi_noise_constructor_defaults_select_m03_ifov_one():
    """Standalone noise application defaults to the M03 IFOV-one factor."""
    result = _iasi_result()

    IASI().add_noise(result, rng=np.random.default_rng(1))

    assert result.processing_history[-1]["noise_factor"] == "M03_1"


def test_registry_constructs_iasi_and_rejects_unknown_names(tiny_iasi_ils):
    """Configuration names resolve through one explicit registry."""
    instrument = create_instrument("IASI", ils_file=tiny_iasi_ils)

    assert isinstance(instrument, IASI)
    assert isinstance(instrument, ConvolvingInstrument)
    with pytest.raises(ValueError, match="Available instruments: iasi"):
        create_instrument("unknown")


def test_configured_convolution_strategy_resolves_from_generic_inputs(tiny_iasi_ils):
    """The runner selector delegates IASI parameter extraction to the instrument."""
    configured = _resolve_convolving_instrument(
        {
            "convolve": "iasi",
            "iasi_ils": tiny_iasi_ils,
            "add_noise": None,
            "iasi_instrument": "unused-value",
            "iasi_ifov": 0,
        }
    )

    assert isinstance(configured, IASI)
    assert configured.ils_file == tiny_iasi_ils
    assert configured.iasi_instrument == "M03"
    assert configured.iasi_ifov == 1


def test_instrument_resolution_defaults_to_none():
    """An omitted configured instrument leaves the radiative-transfer result raw."""
    assert _resolve_convolving_instrument({"convolve": None}) is None


def test_iasi_noise_resolution_is_independent_from_convolution():
    """The noise switch configures IASI without enabling ILS convolution."""
    configured = _resolve_noise_instrument(
        {
            "convolve": None,
            "add_noise": "iasi",
            "iasi_instrument": "M02",
            "iasi_ifov": 3,
        }
    )

    assert isinstance(configured, IASI)
    assert configured.iasi_instrument == "M02"
    assert configured.iasi_ifov == 3


def test_iasi_extracts_tailored_runner_inputs_inside_instrument_package(
    tiny_iasi_ils,
):
    """The generic selector does not need to know the tailored runner's ILS key."""
    instrument = create_instrument_from_inputs(
        {
            "convolve": "iasi",
            "ils": tiny_iasi_ils,
        },
        selector="convolve",
    )

    assert isinstance(instrument, IASI)
    assert instrument.ils_file == tiny_iasi_ils


def test_iasi_uses_packaged_ils_when_no_override_is_configured():
    """Configured IASI processing works without an external response file."""
    instrument = create_instrument_from_inputs(
        {"convolve": "iasi"}, selector="convolve"
    )
    result = SRFM()
    result.wvnm = np.arange(998.0, 1002.25, 0.25)
    result.uu = np.zeros((result.wvnm.size, 1, 1, 1))
    result.uu[result.wvnm.size // 2, 0, 0, 0] = 1.0

    returned = instrument.convolve(result)

    assert returned is result
    assert instrument.ils_file is None
    assert instrument.ils_source == "package:srfm/RFM/rfm_files/iasi.ils"
    assert result.processing_history[-1]["ils_file"] == instrument.ils_source
    assert result.uu[:, 0, 0, 0].sum() == pytest.approx(1.0)


def test_registry_resolves_noise_independently_from_convolution(tiny_iasi_ils):
    """Each operation selects its instrument through its own input field."""
    noise_instrument = create_instrument_from_inputs(
        {
            "convolve": None,
            "add_noise": "iasi",
            "iasi_ils": tiny_iasi_ils,
        },
        selector="add_noise",
    )

    assert isinstance(noise_instrument, IASI)
    assert (
        create_instrument_from_inputs(
            {"convolve": None, "add_noise": "iasi"},
            selector="convolve",
        )
        is None
    )
