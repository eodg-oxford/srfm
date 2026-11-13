"""Scientific and interface checks for CrIS and the HIRAS generations."""

import numpy as np
import pytest
from scipy.integrate import trapezoid

from srfm.forward_model import SRFM
from srfm.instruments.cris import CrIS
from srfm.instruments.hiras import HIRAS
from srfm.instruments.hiras_2 import HIRAS2
from srfm.instruments.interferometer import _integrate_channels


def inputs(**changes):
    values = {
        "spc_units": "cm-1",
        "spc_res": 0.005,
        "spc_wvnmlo": 960.0,
        "spc_wvnmhi": 1040.0,
        "fin_wvnmlo": 999.9,
        "fin_wvnmhi": 1000.1,
        "fin_res": None,
    }
    values.update(changes)
    return values


@pytest.mark.parametrize("instrument", (CrIS, HIRAS, HIRAS2))
@pytest.mark.parametrize("apodization", (None, "Hamming"))
def test_flat_radiance_is_preserved(instrument, apodization):
    values = inputs(apodization=apodization)
    processor = instrument.from_inputs(values)
    grid = processor.output_grid(values)
    source_grid = np.linspace(960.0, 1040.0, 16001)
    result = SRFM(
        wvnm=source_grid,
        uu=np.full((len(source_grid), 2), 0.1),
        retained_outputs=frozenset({"uu"}),
    )
    processor.process_to_grid(result, grid)
    np.testing.assert_allclose(grid, [1000.0])
    np.testing.assert_allclose(result.uu, 0.1, atol=2e-8)


def test_hiras2_overlap_and_hamming_grid():
    values = inputs(spc_wvnmlo=600.0, spc_wvnmhi=2600.0,
                    fin_wvnmlo=None, fin_wvnmhi=None)
    for platform in ("FY3E", "FY3F", "FY3H"):
        native = HIRAS2.from_inputs(values | {"hiras2_platform": platform})
        grid = native.output_grid(values)
        assert len(grid) == 3045
        assert np.all(np.diff(grid) > 0)
        hamming_values = values | {"hiras2_platform": platform,
                                   "apodization": "Hamming"}
        hamming = HIRAS2.from_inputs(hamming_values)
        assert len(hamming.output_grid(hamming_values)) == 3041


def test_grid_error_explains_native_spacing():
    values = inputs(fin_res=0.25)
    with pytest.raises(ValueError, match="fin_res=0.25.*0.625.*fin_res=None"):
        CrIS.from_inputs(values).output_grid(values)


def test_narrow_line_matches_direct_quadrature():
    grid = np.linspace(970.0, 1030.0, 120001)
    line = np.exp(-0.5 * ((grid - 1000.123) / 0.003) ** 2)
    actual = _integrate_channels(grid, line, np.array([1000.0]), None, None, 30.0)[0]
    response = np.sinc((grid - 1000.0) / 0.625)
    reference = trapezoid(line * response, grid) / trapezoid(response, grid)
    assert actual == pytest.approx(reference, rel=3e-4)


def test_cris_noise_is_seeded_and_fov_specific():
    values = inputs(apodization="Hamming", noise_seed=12)
    processor = CrIS.from_inputs(values)
    source_grid = np.linspace(960.0, 1040.0, 16001)
    outputs = []
    for _ in range(2):
        result = SRFM(
            wvnm=source_grid,
            uu=np.full(len(source_grid), 0.1),
            retained_outputs=frozenset({"uu"}),
        )
        processor.process_to_grid(result, processor.output_grid(values))
        processor.add_noise(result)
        outputs.append(result.uu.copy())
    np.testing.assert_array_equal(*outputs)
    assert outputs[0][0] != 0.1


def test_fy3d_noise_units_require_explicit_confirmation():
    values = inputs()
    processor = HIRAS.from_inputs(values)
    source_grid = np.linspace(960.0, 1040.0, 16001)
    result = SRFM(
        wvnm=source_grid,
        uu=np.full(len(source_grid), 0.1),
        retained_outputs=frozenset({"uu"}),
    )
    processor.process_to_grid(result, processor.output_grid(values))
    with pytest.raises(ValueError, match="units.*unverified"):
        processor.add_noise(result)


@pytest.mark.parametrize(
    "instrument,settings",
    ((HIRAS, {"hiras_noise_units": "mW"}), (HIRAS2, {})),
)
def test_fengyun_noise_profiles_can_be_used(instrument, settings):
    values = inputs(**settings, noise_seed=42)
    processor = instrument.from_inputs(values)
    source_grid = np.linspace(960.0, 1040.0, 16001)
    result = SRFM(
        wvnm=source_grid,
        uu=np.full(len(source_grid), 0.1),
        retained_outputs=frozenset({"uu"}),
    )
    processor.process_to_grid(result, processor.output_grid(values))
    processor.add_noise(result)
    assert result.noise_applied
    assert np.all(np.isfinite(result.uu))


def test_fy3f_requires_its_own_noise_source():
    values = inputs(hiras2_platform="FY3F")
    processor = HIRAS2.from_inputs(values)
    source_grid = np.linspace(960.0, 1040.0, 16001)
    result = SRFM(
        wvnm=source_grid,
        uu=np.full(len(source_grid), 0.1),
        retained_outputs=frozenset({"uu"}),
    )
    processor.process_to_grid(result, processor.output_grid(values))
    with pytest.raises(ValueError, match="matching noise_file"):
        processor.add_noise(result)
