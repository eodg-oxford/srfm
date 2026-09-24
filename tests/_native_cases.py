"""Child-process entry points for tests that call f2py extensions."""

from __future__ import annotations

from pathlib import Path
import pickle
import sys
import tempfile
import traceback

import numpy as np


def _disort_case(payload):
    """Run a minimal compiled DISORT calculation.

    The child-process boundary converts native crashes into test failures.

    Args:
        payload: Mapping containing the requested floating-point precision.

    Returns:
        Serializable DISORT status and radiance output.
    """
    from srfm.forward_model import DISORT

    model = DISORT(disort_input={}, disort_out={})
    model.add_disort_empty_input()
    model.set_maxcmu(4)
    model.set_maxmom(4)
    model.set_rhoq(np.zeros((2, 3, 4)))
    model.set_rhou(np.zeros((1, 3, 4)))
    model.set_bemst(np.zeros(2))
    model.set_prnt([False] * 5)
    model.set_plank(True)
    model.set_dtauc_manually(np.array([0.1]))
    model.set_ssalb_manually(np.array([0.0]))
    model.set_pmom_manually(np.zeros((5, 1)))
    model.set_temper(np.array([270.0, 290.0]))
    model.set_wvnm_range(999.5, 1000.5)
    model.set_utau([0.0])
    model.set_umu0(0.5)
    model.set_umu([0.8])
    model.set_btemp(290.0)
    model.set_ttemp(270.0)
    model.set_temis(1.0)
    model.set_wvnm(1000.0)
    model.set_wvl(10.0)
    model.initialize_disort_output_arrays()
    model.run_disort(prec=payload["precision"], adjust_maxcmu=False)
    output = model.disort_out[1000.0]
    return {"status": model.status, "uu": output["uu"]}


def _getmom_case(payload):
    """Evaluate native DISORT HG moments as a convention oracle."""
    if payload.get("precision", "double") == "single":
        from srfm.DISORT import disort_module_s as module
    else:
        from srfm.DISORT_dbl import disort_module_d as module
    nmom = int(payload["nmom"])
    workspace = np.zeros(nmom + 1)
    return module.getmom(
        iphas=3,
        gg=float(payload["asymmetry"]),
        nmom=nmom,
        pmom=workspace,
    )


def _mie_case(payload):
    """Run a minimal compiled Mie optical-properties calculation.

    Synthetic inputs keep native scattering coverage deterministic.

    Args:
        payload: Mapping selecting optional refractive-index inputs.

    Returns:
        Serializable optical-properties output from the native solver.
    """
    from srfm import ARIA_module, optical_properties
    from srfm.size_distribution import LogNormalDistribution

    if payload.get("ri_file"):
        n, k = ARIA_module.read_ri_file(
            payload["ri_file"], wave=[2.0], mode="wavelength"
        )
        wavelength = 2.0
        refractive_index = n - 1j * k
        concentration = 2
    else:
        wavelength = 10.0
        refractive_index = 1.5 - 0.01j
        concentration = 10

    distribution = LogNormalDistribution(
        n=concentration, r=payload.get("radius", 0.3), s=payload.get("spread", 1.5)
    )
    return optical_properties.ewp_hs(
        wavelengths=np.array([wavelength]),
        composition="ri",
        distribution=distribution,
        refractive_index=refractive_index,
        angle=np.array([0.0, 90.0, 180.0]),
        legendre_coefficients_flag=False,
        radii=12,
        eta=1e-5,
    )


def _rfm_flag_matrix_case(payload):
    """Exercise many RFM ``*FLG`` configurations in one native process.

    Every deliberately incomplete driver stops immediately after flag parsing.
    This isolates the option parser and its compatibility rules from external
    HITRAN, CIA, LUT, SVD, FOV, and instrument files.  Reusing the process also
    verifies that the f2py wrapper resets Fortran units and flag state after a
    rejected configuration.

    Args:
        payload: Mapping whose ``cases`` value contains ``id`` and ``flags``.

    Returns:
        Mapping from case identifier to native status and complete RFM log.
    """
    from srfm.RFM import rfm_py

    if rfm_py is None:
        raise RuntimeError("RFM extension is unavailable")

    results = {}
    original = Path.cwd()
    with tempfile.TemporaryDirectory(prefix="srfm-rfm-flags-") as directory:
        root = Path(directory)
        try:
            # RFM writes its diagnostic log in the process working directory.
            # The child-process boundary makes this process-global operation safe.
            import os

            os.chdir(root)
            for case in payload["cases"]:
                flags = tuple(case["flags"])
                driver_lines = (
                    "*HDR",
                    f"  pytest RFM flag configuration {case['id']}",
                    "*FLG",
                    "  " + " ".join(flags),
                    "*END",
                )
                status = int(rfm_py.rfm_run(driver_lines=driver_lines))
                log_path = root / "rfm.log"
                log = log_path.read_text(encoding="utf-8", errors="replace")
                log_path.unlink()
                results[case["id"]] = {"status": status, "log": log}
        finally:
            os.chdir(original)
    return results


def _e2e_case(payload):
    """Run one complete SRFM pathway with the selected top-level runner.

    All runners return the same serializable result contract to pytest.

    Args:
        payload: Mapping containing inputs and an optional runner identifier.

    Returns:
        Serializable spectral results and effective input keys.
    """
    from srfm.inputs import Inputs

    runner_name = payload.get("runner", "main")
    if runner_name == "main":
        from srfm.main import run_srfm
    elif runner_name == "oxharp":
        from srfm.oxharp_main import run_srfm
    elif runner_name == "iasi":
        from srfm.iasi_main import run_srfm
    else:
        raise ValueError(f"Unknown E2E runner: {runner_name!r}")

    if payload.get("driver_path"):
        inputs = Inputs()
        inputs.read_srfm_drv(payload["driver_path"])
    else:
        inputs = Inputs(**payload["values"])

    model = run_srfm(inputs)
    return {
        "wvnm": model.wvnm,
        "uu": getattr(model, "uu", None),
        "bbt": getattr(model, "bbt", None),
        "flup": getattr(model, "flup", None),
        "rfldir": getattr(model, "rfldir", None),
        "rfldn": getattr(model, "rfldn", None),
        "output_values": getattr(model, "output_values", None),
        "input_keys": tuple(inputs.values),
    }


CASES = {
    "disort": _disort_case,
    "getmom": _getmom_case,
    "mie": _mie_case,
    "rfm_flag_matrix": _rfm_flag_matrix_case,
    "e2e": _e2e_case,
}


def main() -> int:
    """Dispatch a native test case and pickle its result for pytest.

    Exceptions are printed for the parent process and converted to failure.

    Returns:
        Process exit status, with nonzero indicating a Python exception.
    """
    case, args_name, result_name = sys.argv[1:4]
    with Path(args_name).open("rb") as handle:
        payload = pickle.load(handle)
    try:
        result = CASES[case](payload)
    except Exception:
        traceback.print_exc(file=sys.stderr)
        return 1
    with Path(result_name).open("wb") as handle:
        pickle.dump(result, handle)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
