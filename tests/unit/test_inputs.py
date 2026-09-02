from pathlib import Path

import pytest

from srfm.inputs import Inputs


pytestmark = pytest.mark.unit


def test_inputs_constructor_collects_keyword_values():
    inp = Inputs(alpha=1, nested={"beta": 2})
    assert inp.values == {"alpha": 1, "nested": {"beta": 2}}


@pytest.mark.parametrize("path_type", [str, Path])
def test_read_srfm_driver_table(tmp_path, path_type):
    driver = tmp_path / "driver.py"
    driver.write_text("inputs = {'answer': 42, 'flags': ('OPT',)}\n", encoding="utf-8")
    inp = Inputs()
    inp.read_srfm_drv(path_type(driver), validate=False)
    assert inp.values["answer"] == 42


def test_driver_modules_do_not_collide_in_import_cache(tmp_path):
    first = tmp_path / "first.py"
    second = tmp_path / "second.py"
    first.write_text("inputs = {'value': 1}\n", encoding="utf-8")
    second.write_text("inputs = {'value': 2}\n", encoding="utf-8")
    inp = Inputs()
    inp.read_srfm_drv(first, validate=False)
    inp.read_srfm_drv(second, validate=False)
    assert inp.values == {"value": 2}


@pytest.mark.parametrize("value", [123, None, object()])
def test_driver_path_type_validation(value):
    with pytest.raises(TypeError, match="str or pathlib.Path"):
        Inputs().read_srfm_drv(value)


def test_missing_driver_file_rejected(tmp_path):
    with pytest.raises(FileNotFoundError, match="not found"):
        Inputs().read_srfm_drv(tmp_path / "missing.py")


@pytest.mark.parametrize(
    ("source", "error"),
    [("value = 1\n", ValueError), ("inputs = []\n", TypeError)],
)
def test_malformed_srfm_driver_tables(tmp_path, source, error):
    path = tmp_path / "bad_driver.py"
    path.write_text(source, encoding="utf-8")
    with pytest.raises(error):
        Inputs().read_srfm_drv(path)


def test_read_oxharp_driver_merges_state_and_ancillary(tmp_path):
    path = tmp_path / "oxharp.py"
    path.write_text(
        "STATE = {'shared': {'state': 1}, 'a': 2}\n"
        "ANCILLARY = {'shared': {'ancillary': 3}, 'b': 4}\n",
        encoding="utf-8",
    )
    inp = Inputs()
    inp.read_oxharp_drv(path)
    assert inp.values == {"shared": {"state": 1, "ancillary": 3}, "a": 2, "b": 4}
