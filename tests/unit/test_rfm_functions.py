from pathlib import Path

import numpy as np
import pytest

from srfm import rfm_functions


pytestmark = pytest.mark.unit


def test_read_regular_rfm_spectrum(tmp_path):
    path = tmp_path / "regular.asc"
    path.write_text(
        "!header one\n!header two\n!header three\n"
        "3 1000.0 0.5 1001.0 'OPTICAL DEPTH'\n"
        "0.1 0.2\n0.3\n",
        encoding="utf-8",
    )
    output = rfm_functions.read_output(path)
    assert [output[f"header{i}"] for i in range(1, 4)] == ["header one", "header two", "header three"]
    np.testing.assert_allclose(output["WNO"], [1000, 1000.5, 1001])
    np.testing.assert_allclose(output["SPC"], [0.1, 0.2, 0.3])
    assert output["LABSPC"] == "OPTICAL DEPTH"


def test_read_irregular_rfm_spectrum(tmp_path):
    path = tmp_path / "irregular.asc"
    path.write_text(
        "!one\n!two\n!three\n3 1000 0 1003 'IRREGULAR'\n"
        "1000 0.1\n1001.5 0.2\n1003 0.4\n",
        encoding="utf-8",
    )
    output = rfm_functions.read_output(path)
    np.testing.assert_allclose(output["WNO"], [1000, 1001.5, 1003])
    np.testing.assert_allclose(output["SPC"], [0.1, 0.2, 0.4])


def test_read_regular_spectrum_validates_point_count(tmp_path):
    path = tmp_path / "truncated.asc"
    path.write_text("!1\n!2\n!3\n3 1 1 3 'X'\n0.1 0.2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="NPNT"):
        rfm_functions.read_output(path)


def test_read_profile_sections_and_comma_separators(tmp_path):
    path = tmp_path / "prf.asc"
    path.write_text("*HGT [km]\n0, 1, 2\n*PRE [mb]\n1000 900 800\n*END\n", encoding="utf-8")
    assert rfm_functions.read_output_prf(path) == {
        "HGT [km]": [0.0, 1.0, 2.0],
        "PRE [mb]": [1000.0, 900.0, 800.0],
    }


def _driver_input():
    return {"HDR": "pytest", "FLG": "OPT NAD", "SPC": "1000 1001 1", "GAS": "F11", "ATM": "tiny.atm", "TAN": "1.0"}


def test_construct_driver_table_contents_and_does_not_mutate_input(tmp_path):
    source = _driver_input()
    original = source.copy()
    rfm_functions.construct_rfm_driver_table(source, tmp_path)
    text = (tmp_path / "rfm.drv").read_text(encoding="utf-8")
    assert text.startswith("*HDR\n  pytest\n*FLG")
    assert "*TAN\n  1.0" in text
    assert text.endswith("*END")
    assert source == original


def test_construct_driver_table_requires_mandatory_and_one_geometry_section(tmp_path):
    for invalid in ({}, _driver_input() | {"DIM": "1"}):
        with pytest.raises(ValueError):
            rfm_functions.construct_rfm_driver_table(invalid, tmp_path)


def test_construct_driver_force_false_preserves_previous_file(tmp_path):
    (tmp_path / "rfm.drv").write_text("old", encoding="utf-8")
    rfm_functions.construct_rfm_driver_table(_driver_input(), tmp_path, force=False)
    assert (tmp_path / "old_rfm.drv").read_text(encoding="utf-8") == "old"
    assert (tmp_path / "rfm.drv").read_text(encoding="utf-8").startswith("*HDR")


def test_construct_levels_file_formats_and_preserves_old_copy(tmp_path):
    rfm_functions.construct_rfm_output_levels_file([2, 0, 1], tmp_path)
    assert (tmp_path / "alts.lev").read_text(encoding="utf-8") == "!Output levels\n2.0000    0.0000    1.0000    \n*END"
    rfm_functions.construct_rfm_output_levels_file([3, 4], tmp_path, force=False)
    assert (tmp_path / "old_alts.lev").exists()


def test_construct_irregular_grid_file(tmp_path):
    relative = rfm_functions.construct_rfm_grid_file([1000, 1000.5, 1002], rfm_fldr=tmp_path)
    assert relative == "./rfm_files/grid.spc"
    lines = (tmp_path / "rfm_files" / "grid.spc").read_text(encoding="utf-8").splitlines()
    assert lines[3] == "3 1000.0000 0 1002.0000"
    assert lines[-1] == "1002.0000 0"


def test_atmosphere_write_read_round_trip(tiny_atmosphere):
    output = rfm_functions.read_atm_file(tiny_atmosphere)
    np.testing.assert_allclose(output["HGT [km]"], [0, 1, 2])
    np.testing.assert_allclose(output["PRE [mb]"], [1013, 900, 800])
    np.testing.assert_allclose(output["F11 [ppmv]"], [1, 1, 1])


def test_write_atmosphere_validates_height_and_lengths(tmp_path):
    with pytest.raises(ValueError, match="one key"):
        rfm_functions.write_atm_file({"PRE [mb]": [1]}, tmp_path / "x.atm")
    with pytest.raises(ValueError, match="expected 2"):
        rfm_functions.write_atm_file({"HGT [km]": [0, 1], "TEM [K]": [280]}, tmp_path / "x.atm")


def test_write_xsc_uses_hitran_fixed_width_format(tiny_xsc_file):
    lines = tiny_xsc_file.read_text(encoding="utf-8").splitlines()
    assert lines[0][:20].strip() == "F11"
    assert float(lines[0][20:30]) == pytest.approx(999)
    assert float(lines[0][30:40]) == pytest.approx(1001)
    assert int(lines[0][40:47]) == 3
    np.testing.assert_allclose([float(value) for value in lines[1].split()], [1e-20, 2e-20, 1e-20])


def test_write_xsc_validates_required_fields_and_lengths(tmp_path):
    with pytest.raises(ValueError, match="missing required"):
        rfm_functions.write_xsc_file({"x": {"molec": "F11"}}, tmp_path / "x.xsc")
    with pytest.raises(ValueError, match="length"):
        rfm_functions.write_xsc_file(
            {"x": {"molec": "F11", "low_spc": 1, "upp_spc": 2, "npts": 2, "temperature": 280, "pressure": 760, "beta_ext": [1]}},
            tmp_path / "x.xsc",
        )
