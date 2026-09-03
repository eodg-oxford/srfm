from pathlib import Path

import pytest

from srfm.inputs import Inputs

pytestmark = pytest.mark.unit


def test_inputs_constructor_collects_keyword_values():
    """Verify that the constructor stores arbitrary keyword inputs.

    Nested mappings must be retained unchanged for later schema validation.
    """
    inp = Inputs(alpha=1, nested={"beta": 2})
    assert inp.values == {"alpha": 1, "nested": {"beta": 2}}


@pytest.mark.parametrize("path_type", [str, Path])
def test_read_srfm_driver_table(tmp_path, path_type):
    """Verify SRFM driver paths accept strings and ``Path`` objects.

    Validation is disabled because this test isolates module loading.

    Args:
        tmp_path: Pytest temporary-path fixture.
        path_type: Parameterized path representation.
    """
    driver = tmp_path / "driver.py"
    driver.write_text("inputs = {'answer': 42, 'flags': ('OPT',)}\n", encoding="utf-8")
    inp = Inputs()
    inp.read_srfm_drv(path_type(driver), validate=False)
    assert inp.values["answer"] == 42


def test_driver_modules_do_not_collide_in_import_cache(tmp_path):
    """Verify successive driver modules use unique import identities.

    Loading a second file must replace values rather than reuse the first
    module from Python's import cache.

    Args:
        tmp_path: Pytest temporary-path fixture.
    """
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
    """Verify unsupported driver path types raise a clear error.

    Only strings and ``Path`` instances have defined loading semantics.

    Args:
        value: Parameterized non-path object.
    """
    with pytest.raises(TypeError, match="str or pathlib.Path"):
        Inputs().read_srfm_drv(value)


def test_missing_driver_file_rejected(tmp_path):
    """Verify a nonexistent driver path is rejected.

    The failure should occur before import machinery is invoked.

    Args:
        tmp_path: Pytest temporary-path fixture.
    """
    with pytest.raises(FileNotFoundError, match="not found"):
        Inputs().read_srfm_drv(tmp_path / "missing.py")


@pytest.mark.parametrize(
    ("source", "error"),
    [("value = 1\n", ValueError), ("inputs = []\n", TypeError)],
)
def test_malformed_srfm_driver_tables(tmp_path, source, error):
    """Verify malformed driver modules raise the documented exception.

    Missing mappings and mappings of the wrong type are covered separately.

    Args:
        tmp_path: Pytest temporary-path fixture.
        source: Parameterized invalid driver source code.
        error: Expected exception type.
    """
    path = tmp_path / "bad_driver.py"
    path.write_text(source, encoding="utf-8")
    with pytest.raises(error):
        Inputs().read_srfm_drv(path)


def test_read_oxharp_driver_merges_state_and_ancillary(tmp_path):
    """Verify OXHARP state and ancillary mappings merge recursively.

    Values unique to either section and nested shared values must survive.

    Args:
        tmp_path: Pytest temporary-path fixture.
    """
    path = tmp_path / "oxharp.py"
    path.write_text(
        "STATE = {'shared': {'state': 1}, 'a': 2}\n"
        "ANCILLARY = {'shared': {'ancillary': 3}, 'b': 4}\n",
        encoding="utf-8",
    )
    inp = Inputs()
    inp.read_oxharp_drv(path)
    assert inp.values == {"shared": {"state": 1, "ancillary": 3}, "a": 2, "b": 4}
