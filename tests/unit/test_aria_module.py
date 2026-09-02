from pathlib import Path

import numpy as np
import pytest

from srfm.ARIA_module import (
    RI,
    ReadError,
    find_ri_files,
    get_ri_filepathname,
    read_ri_file,
)


pytestmark = pytest.mark.unit


def test_find_ri_files_recurses_and_filters_extensions(tmp_path):
    nested = tmp_path / "nested"
    nested.mkdir()
    (tmp_path / "a.ri").write_text("", encoding="utf-8")
    (nested / "b.ri").write_text("", encoding="utf-8")
    (nested / "ignored.txt").write_text("", encoding="utf-8")
    assert {Path(path).name for path in find_ri_files(tmp_path)} == {"a.ri", "b.ri"}


def test_bundled_generic_composition_resolves():
    assert Path(get_ri_filepathname("ash")).name == "eyjafjallajokull-ash_Reed.ri"


def test_unknown_bundled_composition_fails_clearly():
    with pytest.raises(FileNotFoundError, match="not found"):
        get_ri_filepathname("not-a-real-composition")


def test_read_wavelength_file_and_derive_wavenumber(tiny_ri_file):
    ri = RI()
    ri.read(tiny_ri_file)
    assert ri.header["FORMAT"] == "wavl n k"
    np.testing.assert_allclose(ri.data["wavn"], [10000, 5000, 2500])
    wave, n, k = ri.select()
    np.testing.assert_array_equal(wave, [1, 2, 4])
    np.testing.assert_allclose(n, [1.4, 1.5, 1.7])
    np.testing.assert_allclose(k, [0.01, 0.02, 0.04])


def test_read_wavenumber_file_and_derive_wavelength(descending_ri_file):
    wave, n, k = read_ri_file(descending_ri_file, mode="wavenumber")
    np.testing.assert_array_equal(wave, [10000, 5000, 2500])
    np.testing.assert_allclose(n, [1.4, 1.5, 1.7])
    np.testing.assert_allclose(k, [0.01, 0.02, 0.04])


@pytest.mark.parametrize("fixture_name", ["tiny_ri_file", "descending_ri_file"])
def test_interpolation_handles_ascending_and_descending_data(request, fixture_name):
    path = request.getfixturevalue(fixture_name)
    mode = "wavelength" if fixture_name == "tiny_ri_file" else "wavenumber"
    target = [1.5, 3.0] if mode == "wavelength" else [3750, 7500]
    n, k = read_ri_file(path, wave=target, mode=mode)
    expected_n = [1.45, 1.6] if mode == "wavelength" else [1.6, 1.45]
    expected_k = [0.015, 0.03] if mode == "wavelength" else [0.03, 0.015]
    np.testing.assert_allclose(n, expected_n)
    np.testing.assert_allclose(k, expected_k)


def test_out_of_range_error_clip_and_nan(tiny_ri_file):
    ri = RI()
    ri.read(tiny_ri_file)
    with pytest.raises(ValueError, match="outside"):
        ri.select([0.5, 2], out_of_range="error")
    n_clip, k_clip = ri.select([0.5, 5], out_of_range="clip")
    np.testing.assert_allclose(n_clip, [1.4, 1.7])
    np.testing.assert_allclose(k_clip, [0.01, 0.04])
    n_nan, k_nan = ri.select([0.5, 2, 5], out_of_range="nan")
    np.testing.assert_allclose(n_nan, [np.nan, 1.5, np.nan], equal_nan=True)
    np.testing.assert_allclose(k_nan, [np.nan, 0.02, np.nan], equal_nan=True)


@pytest.mark.parametrize("mode", ["frequency", "", None])
def test_invalid_selection_mode_rejected(tiny_ri_file, mode):
    ri = RI()
    ri.read(tiny_ri_file)
    with pytest.raises(ValueError, match="Invalid mode"):
        ri.select([2], mode=mode)


def test_invalid_out_of_range_policy_rejected_even_for_in_range_data(tiny_ri_file):
    ri = RI()
    ri.read(tiny_ri_file)
    with pytest.raises(ValueError, match="out_of_range"):
        ri.select([2], out_of_range="extrapolate")


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("1 1.4 .01\n", "No header"),
        ("# DESCRIPTION = empty\n", "No data"),
        ("# DESCRIPTION = missing format\n1 1.4 .01\n", "No FORMAT"),
        ("# FORMAT = wavl n\n1 1.4\n", "requires n and k"),
        ("# FORMAT = n k\n1.4 .01\n", "requires wavl or wavn"),
        ("# FORMAT = wavl n mystery k\n1 1.4 9 .01\n", "Unknown FORMAT"),
        ("# FORMAT = wavl n k\n1 1.4\n", "expected 3"),
        ("# FORMAT = wavl n k\n1 nope .01\n", "Non-numeric"),
        ("# FORMAT = wavl n k\n1 1.4 .01\n# COMMENT = late\n", "not contiguous"),
    ],
)
def test_malformed_ri_files_raise_read_error(tmp_path, content, message):
    path = tmp_path / "bad.ri"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(ReadError, match=message):
        RI().read(path)


def test_ri_instance_can_be_reused_without_mixing_data(tmp_path, tiny_ri_file):
    second = tmp_path / "second.ri"
    second.write_text("# FORMAT = wavl n k\n10 2.0 0.2\n", encoding="utf-8")
    ri = RI()
    ri.read(tiny_ri_file)
    ri.read(second)
    assert ri.data["n"] == [2.0]


def test_load_refractive_indices_uses_composition_lookup(monkeypatch, tiny_ri_file):
    monkeypatch.setattr("srfm.ARIA_module.get_ri_filepathname", lambda _: str(tiny_ri_file))
    n, k = RI().load_refractive_indices("synthetic", wave=[2.0])
    np.testing.assert_allclose(n, [1.5])
    np.testing.assert_allclose(k, [0.02])
