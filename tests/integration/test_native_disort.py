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
