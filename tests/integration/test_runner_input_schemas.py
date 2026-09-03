"""Integration tests for top-level runners and their validation contracts."""

from __future__ import annotations

import pytest

from srfm.iasi_main import run_srfm as run_iasi_srfm
from srfm.input_schema import InputValidationError
from srfm.inputs import Inputs
from srfm.oxharp_main import run_srfm as run_oxharp_srfm

pytestmark = pytest.mark.integration


def test_oxharp_runner_uses_derived_geometry_schema_before_side_effects(tmp_path):
    """Verify that ``oxharp_main`` applies its schema before creating output.

    A deliberately incomplete mapping distinguishes the missing derived angle
    from the generic runner's direct-angle requirement.

    Args:
        tmp_path: Pytest temporary-path fixture.
    """
    results = tmp_path / "oxharp-must-not-exist"
    inputs = Inputs(results_fldr=str(results), sza_cos=1.0)

    with pytest.raises(InputValidationError) as caught:
        run_oxharp_srfm(inputs)

    assert "zen_cos: required field is missing" in caught.value.issues
    assert "zen: required field is missing" not in caught.value.issues
    assert not results.exists()


def test_iasi_runner_uses_observation_schema_before_side_effects(tmp_path):
    """Verify that ``iasi_main`` applies its schema before creating output.

    The missing processed-observation fields should be reported without
    requiring generic top-level solar or viewing geometry.

    Args:
        tmp_path: Pytest temporary-path fixture.
    """
    results = tmp_path / "iasi-must-not-exist"
    inputs = Inputs(results_fldr=str(results))

    with pytest.raises(InputValidationError) as caught:
        run_iasi_srfm(inputs)

    assert "iasi_fl: required field is missing" in caught.value.issues
    assert "ils: required field is missing" in caught.value.issues
    assert "sun: required field is missing" not in caught.value.issues
    assert not results.exists()
