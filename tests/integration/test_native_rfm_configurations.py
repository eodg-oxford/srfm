"""Exhaustive tests for the finite RFM ``*FLG`` configuration rules."""

from __future__ import annotations

from pathlib import Path
import re

import pytest

from srfm.input_schema import RFM_FLAG_CODES
from srfm.RFM import rfm_py

pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).resolve().parents[2]
DRVFLG_SOURCE = REPO_ROOT / "src" / "srfm" / "RFM" / "source" / "drvflg_sub.f90"


# Minimal valid contexts for flags that have dependencies in DRVFLG.  All
# other flags are valid by themselves at the *FLG parsing stage.
VALID_FLAG_CONTEXTS = {
    "BFX": ("BFX", "RAD", "SFC"),
    "C32": ("C32", "CTM"),
    "C41": ("C41", "CTM"),
    "C4C": ("C4C", "CTM"),
    "COO": ("COO", "FLX", "ZEN"),
    "FIN": ("FIN", "AVG"),
    "FLX": ("FLX", "ZEN"),
    "FVZ": ("FVZ", "FOV"),
    "HYD": ("HYD", "ZEN"),
    "JTP": ("JTP", "JAC"),
    "MTX": ("MTX", "FLX", "SFC"),
    "NAD": ("NAD", "SFC"),
    "VRT": ("VRT", "FLX", "ZEN"),
}


# One configuration reaches each non-pairwise dependency branch.  The exact
# diagnostic is asserted so a reordered or weakened native rule cannot pass.
DEPENDENCY_VIOLATIONS = (
    (("FLX", "OPT", "SFC"), "F-DRVFLG: FLX and OPT flags also require VRT flag"),
    (("BFX",), "F-DRVFLG: BFX flag also requires BBT, RAD or RJT flags"),
    (("C32",), "F-DRVFLG: C32 flag also requires CTM flag"),
    (("C41",), "F-DRVFLG: C41 flag also requires CTM flag"),
    (("C4C",), "F-DRVFLG: C4C flag also requires CTM flag"),
    (("COO",), "F-DRVFLG: COO flag also requires FLX flag"),
    (("FIN",), "F-DRVFLG: FIN flag also requires ILS or AVG flag"),
    (("ABS", "FLX", "SFC"), "F-DRVFLG: ABS+FLX-MTX flags also requires ZEN or NAD"),
    (("TRA", "FLX", "SFC"), "F-DRVFLG: TRA+FLX-MTX flags also requires ZEN or NAD"),
    (("MTX", "FLX", "RAD"), "F-DRVFLG: MTX+RAD flags also requires SFC flag"),
    (("FLX",), "F-DRVFLG: FLX-ZEN flags also require SFC flag"),
    (("FVZ",), "F-DRVFLG: FVZ flag also requires FOV flag"),
    (("HYD",), "F-DRVFLG: HYD flag also requires NAD or ZEN flag"),
    (("JTP",), "F-DRVFLG: JTP flag also requires JAC flag"),
    (("MTX",), "F-DRVFLG: MTX flag also requires FLX flag"),
    (("NAD",), "F-DRVFLG: NAD flag also requires SFC flag"),
    (("VRT",), "F-DRVFLG: VRT flag also requires FLX flag"),
)


def _pairwise_incompatibilities():
    """Read every explicit incompatible pair from the authoritative parser."""
    source = DRVFLG_SOURCE.read_text(encoding="utf-8")
    matches = re.findall(
        r"ERRMSG = '(F-DRVFLG: ([A-Z0-9]{3}) and ([A-Z0-9]{3}) flags are incompatible)'",
        source,
    )
    return tuple(((left, right), diagnostic) for diagnostic, left, right in matches)


def _run_matrix(require_native, run_native_case, cases):
    """Run flag configurations through the compiled parser in one child."""
    require_native(rfm_py, "RFM")
    payload = {
        "cases": [
            {"id": f"case-{index:03d}", "flags": flags}
            for index, flags in enumerate(cases)
        ]
    }
    results, _ = run_native_case("rfm_flag_matrix", payload)
    return [results[item["id"]] for item in payload["cases"]]


def test_every_rfm_flag_has_a_native_accepted_configuration(
    require_native, run_native_case
):
    """Pass all 56 flags through DRVFLG in a valid dependency context."""
    configurations = [VALID_FLAG_CONTEXTS.get(flag, (flag,)) for flag in RFM_FLAG_CODES]

    # Exercise every alternative in the parser's OR-dependencies as well.
    configurations.extend(
        (
            ("BFX", "BBT"),
            ("BFX", "RAD", "SFC"),
            ("BFX", "RJT"),
            ("FIN", "ILS"),
            ("FIN", "AVG"),
            ("HYD", "NAD", "SFC"),
            ("HYD", "ZEN"),
            ("ABS", "FLX", "NAD", "SFC"),
            ("ABS", "FLX", "ZEN"),
            ("TRA", "FLX", "NAD", "SFC"),
            ("TRA", "FLX", "ZEN"),
        )
    )

    results = _run_matrix(require_native, run_native_case, configurations)

    for flags, result in zip(configurations, results, strict=True):
        # Drivers deliberately end before *SPC. Reaching DRVKEY proves that
        # DRVFLG accepted the configuration and all dependencies were reset.
        assert result["status"] == 1
        assert "F-DRVFLG:" not in result["log"], (flags, result["log"])
        assert "F-DRVKEY: Expected *SPC" in result["log"], (flags, result["log"])


def test_every_rfm_incompatibility_and_dependency_rule_rejects_natively(
    require_native, run_native_case
):
    """Reach every finite incompatibility/dependency branch in DRVFLG."""
    configurations = _pairwise_incompatibilities() + DEPENDENCY_VIOLATIONS
    results = _run_matrix(
        require_native,
        run_native_case,
        [flags for flags, _ in configurations],
    )

    for (flags, diagnostic), result in zip(configurations, results, strict=True):
        assert result["status"] == 1
        assert diagnostic in result["log"], (flags, result["log"])

    # Keep the test matrix exhaustive when native rules are added or removed.
    source = DRVFLG_SOURCE.read_text(encoding="utf-8")
    compatibility_source = source.split(
        "! Check for inconsistent flag combinations", maxsplit=1
    )[1]
    native_rule_diagnostics = set(
        re.findall(r"ERRMSG = '(F-DRVFLG:[^']+)'", compatibility_source)
    )
    tested_diagnostics = {diagnostic for _, diagnostic in configurations}
    assert tested_diagnostics == native_rule_diagnostics


def test_rfm_flag_parser_rejects_malformed_unknown_and_repeated_flags(
    require_native, run_native_case
):
    """Exercise the three syntax/error paths preceding compatibility checks."""
    configurations = (("AB",), ("XYZ",), ("OPT", "OPT"))
    expected = (
        "F-DRVFLG: Flag not a C*3 string",
        "F-DRVFLG: *FLG section contains unrecognised flag=XYZ",
        "F-DRVFLG: *FLG section contains repeated flag: OPT",
    )
    results = _run_matrix(require_native, run_native_case, configurations)
    for diagnostic, result in zip(expected, results, strict=True):
        assert result["status"] == 1
        assert diagnostic in result["log"]
