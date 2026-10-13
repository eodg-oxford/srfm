import numpy as np
import pytest

from srfm.DISORT import disort_module_s
from srfm.DISORT_dbl import disort_module_d

pytestmark = pytest.mark.integration


@pytest.mark.parametrize(
    ("precision", "module", "name"),
    [
        ("single", disort_module_s, "single-precision DISORT"),
        ("double", disort_module_d, "double-precision DISORT"),
    ],
)
def test_python_disort_wrapper_runs_compiled_solver_in_isolated_process(
    precision, module, name, require_native, run_native_case
):
    """Verify that the Python DISORT wrapper runs each native precision.

    Isolation ensures a native solver failure becomes a normal pytest failure
    instead of terminating the main test process.

    Args:
        precision: Precision name passed to the wrapper.
        module: Compiled DISORT module for the selected precision.
        name: Human-readable native module name.
        require_native: Fixture helper for extension availability.
        run_native_case: Fixture helper that isolates native execution.
    """
    require_native(module, name)

    result, _ = run_native_case("disort", {"precision": precision})

    assert result["status"] == "DISORT run completed."
    assert result["uu"].shape == (1, 1, 1)
    assert np.all(np.isfinite(result["uu"]))
    assert result["uu"].item() > 0


def test_native_henyey_greenstein_moments_confirm_analytic_convention(
    require_native, run_native_case
):
    """Check all eight native moment configurations in one isolated process."""
    require_native(disort_module_d, "double-precision DISORT")
    cases = {
        f"nmom={nmom}, asymmetry={asymmetry}": {
            "case": "getmom",
            "payload": {"precision": "double", "asymmetry": asymmetry, "nmom": nmom},
        }
        for nmom in (2, 7)
        for asymmetry in (-0.4, 0.0, 0.35, 0.8)
    }
    results, _ = run_native_case("batch", {"cases": cases})
    assert results.keys() == cases.keys()
    for identifier, configuration in cases.items():
        parameters = configuration["payload"]
        np.testing.assert_allclose(
            results[identifier],
            parameters["asymmetry"] ** np.arange(parameters["nmom"] + 1),
            rtol=1e-13,
            atol=1e-15,
            err_msg=identifier,
        )
