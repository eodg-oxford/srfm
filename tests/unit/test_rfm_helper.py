from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from srfm import rfm_helper

pytestmark = pytest.mark.unit


def test_run_result_ok_property_and_defaults():
    """Verify run-result defaults and success classification.

    Status zero should be successful while preserving optional output defaults.
    """
    success = rfm_helper.RunResult(0, [])
    failure = rfm_helper.RunResult(2, [Path("rfm.log")])
    assert success.ok is True
    assert success.output is None
    assert failure.ok is False


def test_clean_outputs_removes_only_matching_regular_files(tmp_path):
    """Verify cleanup removes only matching regular files.

    Unrelated files and similarly named directories must survive.

    Args:
        tmp_path: Pytest temporary-path fixture.
    """
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
    """Verify RFM logs are emitted to the appropriate scheduler stream.

    Successful logs use stdout while failure diagnostics use stderr.

    Args:
        tmp_path: Pytest temporary-path fixture.
        capsys: Pytest fixture capturing standard streams.
        success: Whether the simulated RFM run succeeded.
        stream_name: Captured stream expected to receive the log.
    """
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
    """Verify failed native runs emit and remove their diagnostic log.

    The raised exception should retain native detail without stale artifacts.

    Args:
        monkeypatch: Pytest fixture replacing the native RFM module.
        tmp_path: Pytest temporary-path fixture.
        capsys: Pytest fixture capturing standard streams.
    """

    class FailingRfm:
        @staticmethod
        def rfm_run(**kwargs):
            """Simulate a failing native call and create its log.

            The stub mirrors the side effects and status of a native failure.

            Args:
                **kwargs: Native-call arguments accepted for API compatibility.

            Returns:
                Nonzero synthetic RFM status.
            """
            Path("rfm.log_failed").write_text(
                "F-RFM: synthetic native failure\n", encoding="utf-8"
            )
            return 7

        @staticmethod
        def _rfm_py_error():
            """Return the synthetic native failure detail.

            The helper mirrors the compiled wrapper's error accessor.

            Returns:
                Human-readable failure message.
            """
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
    """Verify structured driver dataclasses format RFM records.

    Paths, labels, ranges, and assignment lines cover the supported forms.

    Args:
        tmp_path: Pytest temporary-path fixture.
    """
    assert rfm_helper.SpectralRange(1000, 1002, 0.25).as_record() == "1000 1002 0.25"
    assert (
        rfm_helper.SpectralRange(1000, 1002, 1, "band").as_record()
        == "band 1000 1002 1"
    )
    spectral_file = rfm_helper.SpectralFile(tmp_path / "grid.spc", "irregular")
    assert isinstance(spectral_file.path, Path)
    assert spectral_file.as_record().endswith("grid.spc")
    assert rfm_helper.SectionLine("TEM", "=", 280.0).as_record() == "TEM = 280"


@pytest.mark.parametrize("args", [(2, 1, 0.1), (1, 2, 0), (1, float("nan"), 1)])
def test_spectral_range_validation(args):
    """Verify invalid spectral ranges are rejected.

    Reversed bounds, zero resolution, and nonfinite values are covered.

    Args:
        args: Parameterized spectral-range constructor arguments.
    """
    with pytest.raises(ValueError):
        rfm_helper.SpectralRange(*args)


def test_input_catalog_exposes_mandatory_sections_flags_and_aliases():
    """Verify the input catalogue exposes sections, flags, and aliases.

    The catalogue serves documentation and structured driver validation.
    """
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
    """Verify heterogeneous section data normalize to record strings.

    Scalars, mappings, tuples, and structured records share one output form.

    Args:
        payload: Parameterized section value.
        expected: Expected normalized list of records.
    """
    assert rfm_helper._normalize_section_data(payload) == expected


def test_compose_driver_sections_has_deterministic_order_and_end_marker():
    """Verify composed driver sections have stable order and termination.

    Deterministic ordering keeps generated tables readable and reproducible.
    """
    lines = rfm_helper._compose_driver_sections(
        header="pytest",
        flags=["OPT", "NAD"],
        spectral=[rfm_helper.SpectralRange(1000, 1001, 1)],
        gases=["F11"],
        atmosphere=["tiny.atm"],
        tangent=["1"],
        tab_dimensions=None,
        cia=None,
        fin=None,
        fov=None,
        grd=None,
        hit=None,
        ils=None,
        jac=None,
        lev=[0, 1],
        lut=None,
        nte=None,
        obs=None,
        out=None,
        phy=None,
        rej=None,
        sfc=None,
        shp=None,
        svd=None,
        xsc=["f11.xsc"],
        extra_sections=None,
        uses_tab=False,
    )
    assert lines[:4] == ["*HDR", "  pytest", "*FLG", "  OPT NAD"]
    assert lines.index("*TAN") < lines.index("*LEV") < lines.index("*XSC")
    assert lines[-1] == "*END"


def test_compose_driver_rejects_duplicate_extra_section():
    """Verify extra sections cannot duplicate standard driver sections.

    Duplicate RFM sections would make the generated configuration ambiguous.
    """
    kwargs = dict(
        header="x",
        flags=["OPT"],
        spectral=["1 2 1"],
        gases=["F11"],
        atmosphere=["a.atm"],
        tangent=["1"],
        tab_dimensions=None,
        cia=None,
        fin=None,
        fov=None,
        grd=None,
        hit=None,
        ils=None,
        jac=None,
        lev=None,
        lut=None,
        nte=None,
        obs=None,
        out=None,
        phy=None,
        rej=None,
        sfc=None,
        shp=None,
        svd=None,
        xsc=["x.xsc"],
        uses_tab=False,
    )
    with pytest.raises(ValueError, match="more than once"):
        rfm_helper._compose_driver_sections(
            **kwargs, extra_sections={"XSC": "other.xsc"}
        )


def test_run_with_parameters_generates_driver_without_running_native(
    monkeypatch, tmp_path
):
    """Verify structured parameters generate a driver for native execution.

    A stub captures the composed lines without invoking compiled code.

    Args:
        monkeypatch: Pytest fixture replacing native execution.
        tmp_path: Pytest temporary-path fixture.
    """
    captured = {}

    def fake_run_rfm(**kwargs):
        """Capture generated native arguments without running RFM.

        The successful result lets the surrounding orchestration complete.

        Args:
            **kwargs: Generated arguments intended for ``run_rfm``.

        Returns:
            Successful synthetic run result.
        """
        captured.update(kwargs)
        return rfm_helper.RunResult(0, [])

    monkeypatch.setattr(rfm_helper, "run_rfm", fake_run_rfm)
    result = rfm_helper.run_rfm_with_parameters(
        header="pytest",
        flags=("opt", "nad"),
        spectral=(rfm_helper.SpectralRange(1000, 1001, 1),),
        gases=("F11",),
        atmosphere=(tmp_path / "a.atm",),
        tangent=("1",),
        lev=(0, 1),
        xsc=(tmp_path / "f11.xsc",),
        directory=tmp_path,
        output_mode="capture",
        optical_levels=(0, 1),
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
    """Verify structured execution rejects invalid configurations.

    Parameterized cases cover flags and mutually dependent geometry sections.

    Args:
        monkeypatch: Pytest fixture available for execution isolation.
        tmp_path: Pytest temporary-path fixture.
        kwargs: Parameterized invalid configuration updates.
    """
    values = dict(
        header="x",
        flags=("OPT",),
        spectral=("1 2 1",),
        gases=("F11",),
        atmosphere=("a.atm",),
        tangent=("1",),
        directory=tmp_path,
    )
    values.update(kwargs)
    with pytest.raises(ValueError):
        rfm_helper.run_rfm_with_parameters(**values)


class _FakeRfm:
    @staticmethod
    def rfm_get_optical_grid_size(ispc):
        """Return dimensions for the synthetic captured grid.

        The fixed dimensions match arrays returned by the companion stub.

        Args:
            ispc: Requested native spectrum index.

        Returns:
            Synthetic status, level count, and point count.
        """
        return 0, 3, 2

    @staticmethod
    def rfm_get_optical_grid(n_levels, n_points, ispc):
        """Return a synthetic cumulative optical-depth grid.

        Values are chosen so expected differential depths remain simple.

        Args:
            n_levels: Number of requested atmospheric levels.
            n_points: Number of requested spectral points.
            ispc: Requested native spectrum index.

        Returns:
            Native-style status and captured profile arrays.
        """
        return (
            0,
            np.array([0.0, 1.0, 2.0]),
            np.array([1000.0, 900.0, 800.0]),
            np.array([290.0, 280.0, 270.0]),
            np.array([1000.0, 1001.0]),
            np.array([[0.3, 0.6], [0.2, 0.4], [0.0, 0.0]]),
        )


def test_captured_optical_depths_transform_cumulative_levels_to_layers(monkeypatch):
    """Verify captured cumulative depths transform to layer depths.

    Output columns and both cumulative/differential values are asserted.

    Args:
        monkeypatch: Pytest fixture replacing captured native arrays.
    """
    monkeypatch.setattr(rfm_helper, "rfm_py", _FakeRfm())
    with pytest.warns(UserWarning, match="Temperature step"):
        result = rfm_helper.get_captured_optical_depths([0, 1, 2])
    assert list(result.columns[-4:]) == [
        "dOD_1000.0000",
        "dOD_1001.0000",
        "iOD_1000.0000",
        "iOD_1001.0000",
    ]
    np.testing.assert_allclose(
        result[["dOD_1000.0000", "dOD_1001.0000"]], [[0.2, 0.4], [0.1, 0.2]]
    )
    np.testing.assert_allclose(
        result[["iOD_1000.0000", "iOD_1001.0000"]], [[0.2, 0.4], [0.3, 0.6]]
    )
    np.testing.assert_allclose(result["h_avg (km)"], [1.5, 0.5])


def test_captured_optical_depths_report_missing_level(monkeypatch):
    """Verify captured depths report requested levels that are absent.

    Missing native levels should be named in the validation error.

    Args:
        monkeypatch: Pytest fixture replacing captured native arrays.
    """
    monkeypatch.setattr(rfm_helper, "rfm_py", _FakeRfm())
    with pytest.raises(ValueError, match="3"):
        rfm_helper.get_captured_optical_depths([0, 3])


def test_load_driver_lines_preserves_spaces_but_removes_newlines(tmp_path):
    """Verify driver loading preserves spaces while stripping newlines.

    RFM record indentation is significant but platform line endings are not.

    Args:
        tmp_path: Pytest temporary-path fixture.
    """
    path = tmp_path / "rfm.drv"
    path.write_text("*HDR  \r\n  content  \n*END\n", encoding="utf-8")
    assert rfm_helper._load_driver_lines(path) == ["*HDR  ", "  content  ", "*END"]
