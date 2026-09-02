"""Tests for the authoritative schema around the existing driver-table format."""

from __future__ import annotations

from copy import deepcopy
import inspect
from pathlib import Path
import re

import numpy as np
import pytest

from srfm import rfm_helper
from srfm.input_schema import (
    InputValidationError,
    RFM_CONFIG_SCHEMA,
    RFM_DRIVER_SCHEMA,
    RFM_FLAG_CODES,
    SRFM_INPUT_SCHEMA,
    get_srfm_input_schema,
    validate_srfm_inputs,
)
from srfm.inputs import Inputs
from srfm.main import run_srfm


pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
BASIC_DRIVER = REPO_ROOT / "examples" / "basic_example" / "driver_table.py"


@pytest.fixture
def basic_values():
    """Load the existing public example without validation so tests may mutate it."""
    inputs = Inputs()
    inputs.read_srfm_drv(BASIC_DRIVER, validate=False)
    return deepcopy(inputs.values)


def test_existing_basic_driver_validates_without_changing_its_structure(basic_values):
    original_top_level = set(basic_values)
    original_rfm_fields = set(basic_values["driver_inputs"])
    original_layer_fields = {
        name: set(layer) for name, layer in basic_values["scat_lyrs_inputs"].items()
    }

    validated = validate_srfm_inputs(basic_values)

    assert set(validated) == original_top_level
    assert set(validated["driver_inputs"]) == original_rfm_fields
    assert {
        name: set(layer) for name, layer in validated["scat_lyrs_inputs"].items()
    } == original_layer_fields


def test_complete_driver_is_validated_automatically_when_read():
    inputs = Inputs()
    inputs.read_srfm_drv(BASIC_DRIVER)
    assert inputs.values["driver_inputs"]["header"] == "SRFM_standard_run"
    assert "Sulphuric_acid_1" in inputs.values["scat_lyrs_inputs"]


def test_schema_covers_every_rfm_driver_field():
    runtime_parameters = {
        "run_id",
        "clean_before",
        "directory",
        "patterns",
        "output_mode",
        "enable_capture",
        "optical_levels",
        "optical_spectrum_index",
        "optical_match_tol",
    }
    helper_parameters = set(inspect.signature(rfm_helper.run_rfm_with_parameters).parameters)
    assert set(RFM_DRIVER_SCHEMA) == helper_parameters - runtime_parameters


def test_rfm_flag_schema_matches_fortran_driver_parser():
    """Protect every code accepted by DRVFLG's authoritative SELECT CASE."""
    source = (
        REPO_ROOT / "src" / "srfm" / "RFM" / "source" / "drvflg_sub.f90"
    ).read_text(encoding="utf-8")
    fortran_flags = set(
        re.findall(r"^\s*CASE \( '([A-Z0-9]{3})' \)", source, flags=re.MULTILINE)
    )

    assert fortran_flags
    assert set(RFM_FLAG_CODES) == fortran_flags
    assert RFM_DRIVER_SCHEMA["flags"].item_choices == fortran_flags
    assert rfm_helper.FLAG_CODES is RFM_FLAG_CODES


def test_deprecated_cia_flag_remains_accepted_by_schema(basic_values):
    flags = basic_values["driver_inputs"]["flags"]
    basic_values["driver_inputs"]["flags"] = (*flags, "CIA")
    validate_srfm_inputs(basic_values)


def test_public_schema_copy_cannot_change_authoritative_mapping():
    schema_copy = get_srfm_input_schema()
    schema_copy.pop("results_fldr")
    assert "results_fldr" in SRFM_INPUT_SCHEMA
    assert set(RFM_CONFIG_SCHEMA) >= {
        "output_mode",
        "driver_path",
        "generate_driver",
        "capture_files_content",
    }


def test_validation_reports_all_top_level_problems_together():
    with pytest.raises(InputValidationError) as caught:
        validate_srfm_inputs({"results_fldr": 12, "result_folder": "typo"})

    message = str(caught.value)
    assert "result_folder: unknown field" in message
    assert "results_fldr: expected str, PathLike, got int" in message
    assert "base_plots: required field is missing" in message
    assert len(caught.value.issues) > 3


def test_nested_typos_and_invalid_rfm_flags_have_dotted_paths(basic_values):
    basic_values["rfm_config"]["capture_file_content"] = True
    basic_values["driver_inputs"]["flags"] = ("OPT", "LEV", "NOT_A_FLAG")

    with pytest.raises(InputValidationError) as caught:
        validate_srfm_inputs(basic_values)

    assert "rfm_config.capture_file_content: unknown field" in caught.value.issues
    assert any(
        issue.startswith("driver_inputs.flags: unknown RFM flag code(s)")
        for issue in caught.value.issues
    )


@pytest.mark.parametrize("missing_flag", ["OPT", "LEV"])
def test_run_srfm_requires_rfm_optical_capture_flags(basic_values, missing_flag):
    basic_values["driver_inputs"]["flags"] = tuple(
        flag
        for flag in basic_values["driver_inputs"]["flags"]
        if flag != missing_flag
    )
    with pytest.raises(InputValidationError, match=f"requires {missing_flag}"):
        validate_srfm_inputs(basic_values)


def test_run_srfm_rejects_file_only_rfm_mode_before_execution(basic_values):
    basic_values["rfm_config"]["output_mode"] = "files"
    with pytest.raises(InputValidationError, match="requires 'capture'"):
        validate_srfm_inputs(basic_values)


def test_numpy_levels_and_utau_are_valid_existing_input_forms(basic_values):
    basic_values["levels"] = np.array([0.0, 1.0, 2.0])
    basic_values["utau"] = np.array([0.0, 0.5])
    validated = validate_srfm_inputs(basic_values)
    np.testing.assert_array_equal(validated["levels"], [0.0, 1.0, 2.0])
    np.testing.assert_array_equal(validated["utau"], [0.0, 0.5])


def test_flux_only_disort_allows_zero_user_angle_dimensions(basic_values):
    basic_values.update(onlyfl=True, maxumu=0, maxphi=0)
    validate_srfm_inputs(basic_values)

    basic_values["onlyfl"] = False
    with pytest.raises(InputValidationError) as caught:
        validate_srfm_inputs(basic_values)
    assert "maxumu: may be zero only when onlyfl is True" in caught.value.issues
    assert "maxphi: may be zero only when onlyfl is True" in caught.value.issues


@pytest.mark.parametrize(
    ("updates", "problem"),
    [
        ({"fin_res": 0}, "fin_res: must be greater than zero"),
        ({"fin_res": np.nan}, "fin_res: must be finite"),
        ({"spc_wvnmlo": 900, "spc_wvnmhi": 800}, "spc_wvnmlo"),
        ({"convolve_iasi": True, "iasi_ils": None}, "iasi_ils"),
        ({"out_mode": "txt", "rad": False, "bbt": False}, "out_mode"),
        ({"date": (2025, 2, 29)}, "date: day is out of range"),
        ({"header": "x" * 127}, "header: must contain fewer"),
    ],
)
def test_cross_field_validation(basic_values, updates, problem):
    basic_values.update(updates)
    with pytest.raises(InputValidationError, match=problem):
        validate_srfm_inputs(basic_values)


def test_scattering_layer_errors_identify_layer_and_field(basic_values):
    layer = basic_values["scat_lyrs_inputs"]["Ash_1"]
    layer["thick"] = -1
    layer["mass_loading"] = -0.2
    layer["unexpected"] = 1

    with pytest.raises(InputValidationError) as caught:
        validate_srfm_inputs(basic_values)

    assert "scat_lyrs_inputs.Ash_1.unexpected: unknown field" in caught.value.issues
    assert "scat_lyrs_inputs.Ash_1.thick: must be greater than zero" in caught.value.issues
    assert "scat_lyrs_inputs.Ash_1.mass_loading: must be non-negative" in caught.value.issues


@pytest.mark.parametrize(
    ("updates", "problem"),
    [
        ({"dist_type": "gaussian"}, "gaussian is recognized but not implemented"),
        ({"rho": "water"}, "named density must be pumice, glass, mineral, or rock"),
        ({"center_alt": "high", "thick": []}, "center_alt: expected Real"),
    ],
)
def test_layer_validation_rejects_unsupported_or_malformed_values(
    basic_values, updates, problem
):
    basic_values["scat_lyrs_inputs"]["Ash_1"].update(updates)
    with pytest.raises(InputValidationError, match=problem):
        validate_srfm_inputs(basic_values)


def test_run_srfm_revalidates_programmatic_inputs_before_side_effects(tmp_path):
    results = tmp_path / "must-not-be-created"
    inputs = Inputs(results_fldr=str(results), misspelled_field=True)

    with pytest.raises(InputValidationError, match="misspelled_field: unknown field"):
        run_srfm(inputs)

    assert not results.exists()


def test_partial_driver_loading_requires_explicit_validation_opt_out(tmp_path):
    driver = tmp_path / "partial.py"
    driver.write_text("inputs = {'tooling_value': 7}\n", encoding="utf-8")

    with pytest.raises(InputValidationError):
        Inputs().read_srfm_drv(driver)

    inputs = Inputs()
    inputs.read_srfm_drv(driver, validate=False)
    assert inputs.values == {"tooling_value": 7}
