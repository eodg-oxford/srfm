from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from srfm import rfm_helper


pytestmark = pytest.mark.unit


def test_run_result_ok_property_and_defaults():
    success = rfm_helper.RunResult(0, [])
    failure = rfm_helper.RunResult(2, [Path("rfm.log")])
    assert success.ok is True
    assert success.output is None
    assert failure.ok is False


def test_clean_outputs_removes_only_matching_regular_files(tmp_path):
    for name in ("a.asc", "rfm.log", "keep.txt"):
        (tmp_path / name).write_text(name, encoding="utf-8")
    (tmp_path / "directory.asc").mkdir()
    removed = rfm_helper.clean_outputs(tmp_path)
    assert {path.name for path in removed} == {"a.asc", "rfm.log"}
    assert (tmp_path / "keep.txt").exists()
    assert (tmp_path / "directory.asc").is_dir()


@pytest.mark.parametrize(
    ("success", "stream_name"),
    [(True, "out"), (False, "err")],
)
def test_rfm_log_is_emitted_to_scheduler_stream(tmp_path, capsys, success, stream_name):
    log = tmp_path / "rfm.log_test"
    log.write_text("R-RFM: scientific diagnostic\n", encoding="utf-8")

    rfm_helper._emit_rfm_log(tmp_path, "_test", success=success)

    captured = capsys.readouterr()
    selected = getattr(captured, stream_name)
    other = captured.err if stream_name == "out" else captured.out
    assert "--- rfm.log_test ---" in selected
    assert "R-RFM: scientific diagnostic" in selected
    assert other == ""
    assert not log.exists()


def test_failed_rfm_run_emits_and_removes_native_log(monkeypatch, tmp_path, capsys):
    class FailingRfm:
        @staticmethod
        def rfm_run(**kwargs):
            Path("rfm.log_failed").write_text(
                "F-RFM: synthetic native failure\n", encoding="utf-8"
            )
            return 7

        @staticmethod
        def _rfm_py_error():
            return "synthetic failure"

    monkeypatch.setattr(rfm_helper, "rfm_py", FailingRfm())

    with pytest.raises(RuntimeError, match="synthetic failure"):
        rfm_helper._run_rfm_impl(
            run_id="_failed",
            directory=tmp_path,
            driver_lines=("*HDR", "  pytest", "*END"),
            output_mode="capture",
        )

    streams = capsys.readouterr()
    assert streams.out == ""
    assert "F-RFM: synthetic native failure" in streams.err
    assert not (tmp_path / "rfm.log_failed").exists()


def test_driver_dataclasses_format_records(tmp_path):
    assert rfm_helper.SpectralRange(1000, 1002, 0.25).as_record() == "1000 1002 0.25"
    assert rfm_helper.SpectralRange(1000, 1002, 1, "band").as_record() == "band 1000 1002 1"
    spectral_file = rfm_helper.SpectralFile(tmp_path / "grid.spc", "irregular")
    assert isinstance(spectral_file.path, Path)
    assert spectral_file.as_record().endswith("grid.spc")
    assert rfm_helper.SectionLine("TEM", "=", 280.0).as_record() == "TEM = 280"


@pytest.mark.parametrize("args", [(2, 1, 0.1), (1, 2, 0), (1, float("nan"), 1)])
def test_spectral_range_validation(args):
    with pytest.raises(ValueError):
        rfm_helper.SpectralRange(*args)


def test_input_catalog_exposes_mandatory_sections_flags_and_aliases():
    catalog = rfm_helper.get_rfm_input_catalog()
    assert catalog["mandatory_sections"][:5] == ("*HDR", "*FLG", "*SPC", "*GAS", "*ATM")
    assert "OPT" in catalog["flag_codes"]
    assert "CIA" in catalog["flag_codes"]
    assert catalog["aliases"]["*SEC"] == "*TAN"
    assert any(section.key == "*XSC" for section in catalog["sections"])


@pytest.mark.parametrize(
    ("payload", "expected"),
    [
        (None, []),
        (" value ", ["value"]),
        (3.5, ["3.5"]),
        ({"TEM": 280.0}, ["TEM = 280"]),
        (("N2", "O2"), ["N2", "O2"]),
        (((1000, 1002, 1),), ["1000 1002 1"]),
    ],
)
def test_normalise_section_data(payload, expected):
    assert rfm_helper._normalize_section_data(payload) == expected


def test_compose_driver_sections_has_deterministic_order_and_end_marker():
    lines = rfm_helper._compose_driver_sections(
        header="pytest", flags=["OPT", "NAD"], spectral=[rfm_helper.SpectralRange(1000, 1001, 1)], gases=["F11"], atmosphere=["tiny.atm"], tangent=["1"], tab_dimensions=None, cia=None, fin=None, fov=None, grd=None, hit=None, ils=None, jac=None, lev=[0, 1], lut=None, nte=None, obs=None, out=None, phy=None, rej=None, sfc=None, shp=None, svd=None, xsc=["f11.xsc"], extra_sections=None, uses_tab=False
    )
    assert lines[:4] == ["*HDR", "  pytest", "*FLG", "  OPT NAD"]
    assert lines.index("*TAN") < lines.index("*LEV") < lines.index("*XSC")
    assert lines[-1] == "*END"


def test_compose_driver_rejects_duplicate_extra_section():
    kwargs = dict(header="x", flags=["OPT"], spectral=["1 2 1"], gases=["F11"], atmosphere=["a.atm"], tangent=["1"], tab_dimensions=None, cia=None, fin=None, fov=None, grd=None, hit=None, ils=None, jac=None, lev=None, lut=None, nte=None, obs=None, out=None, phy=None, rej=None, sfc=None, shp=None, svd=None, xsc=["x.xsc"], uses_tab=False)
    with pytest.raises(ValueError, match="more than once"):
        rfm_helper._compose_driver_sections(**kwargs, extra_sections={"XSC": "other.xsc"})


def test_run_with_parameters_generates_driver_without_running_native(monkeypatch, tmp_path):
    captured = {}

    def fake_run_rfm(**kwargs):
        captured.update(kwargs)
        return rfm_helper.RunResult(0, [])

    monkeypatch.setattr(rfm_helper, "run_rfm", fake_run_rfm)
    result = rfm_helper.run_rfm_with_parameters(
        header="pytest", flags=("opt", "nad"), spectral=(rfm_helper.SpectralRange(1000, 1001, 1),), gases=("F11",), atmosphere=(tmp_path / "a.atm",), tangent=("1",), lev=(0, 1), xsc=(tmp_path / "f11.xsc",), directory=tmp_path, output_mode="capture", optical_levels=(0, 1)
    )
    assert result.ok
    assert captured["driver_lines"][0:4] == ["*HDR", "  pytest", "*FLG", "  OPT NAD"]
    assert captured["output_mode"] == "capture"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"flags": ()},
        {"flags": ("NOTREAL",)},
        {"flags": ("OPT",), "tangent": None},
        {"flags": ("TAB",), "tangent": ("1",), "tab_dimensions": ("P",)},
        {"flags": ("TAB",), "tangent": None, "tab_dimensions": None},
    ],
)
def test_run_with_parameters_validates_configuration(monkeypatch, tmp_path, kwargs):
    values = dict(header="x", flags=("OPT",), spectral=("1 2 1",), gases=("F11",), atmosphere=("a.atm",), tangent=("1",), directory=tmp_path)
    values.update(kwargs)
    with pytest.raises(ValueError):
        rfm_helper.run_rfm_with_parameters(**values)


class _FakeRfm:
    @staticmethod
    def rfm_get_optical_grid_size(ispc):
        return 0, 3, 2

    @staticmethod
    def rfm_get_optical_grid(n_levels, n_points, ispc):
        return (
            0,
            np.array([0.0, 1.0, 2.0]),
            np.array([1000.0, 900.0, 800.0]),
            np.array([290.0, 280.0, 270.0]),
            np.array([1000.0, 1001.0]),
            np.array([[0.3, 0.6], [0.2, 0.4], [0.0, 0.0]]),
        )


def test_captured_optical_depths_transform_cumulative_levels_to_layers(monkeypatch):
    monkeypatch.setattr(rfm_helper, "rfm_py", _FakeRfm())
    with pytest.warns(UserWarning, match="Temperature step"):
        result = rfm_helper.get_captured_optical_depths([0, 1, 2])
    assert list(result.columns[-4:]) == ["dOD_1000.0000", "dOD_1001.0000", "iOD_1000.0000", "iOD_1001.0000"]
    np.testing.assert_allclose(result[["dOD_1000.0000", "dOD_1001.0000"]], [[0.2, 0.4], [0.1, 0.2]])
    np.testing.assert_allclose(result[["iOD_1000.0000", "iOD_1001.0000"]], [[0.2, 0.4], [0.3, 0.6]])
    np.testing.assert_allclose(result["h_avg (km)"], [1.5, 0.5])


def test_captured_optical_depths_report_missing_level(monkeypatch):
    monkeypatch.setattr(rfm_helper, "rfm_py", _FakeRfm())
    with pytest.raises(ValueError, match="3"):
        rfm_helper.get_captured_optical_depths([0, 3])


def test_load_driver_lines_preserves_spaces_but_removes_newlines(tmp_path):
    path = tmp_path / "rfm.drv"
    path.write_text("*HDR  \r\n  content  \n*END\n", encoding="utf-8")
    assert rfm_helper._load_driver_lines(path) == ["*HDR  ", "  content  ", "*END"]
