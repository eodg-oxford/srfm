"""Authoritative validation schema for the existing SRFM input dictionary.

The schema deliberately validates the dictionary already used by SRFM.  It
does not introduce a replacement file format or rename driver-table fields.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import datetime as dt
from numbers import Integral, Real
import os
from typing import Any

import numpy as np
from .spectral_fields import GRID_UNITS, SOLAR_VALUE_UNITS

_MISSING = object()


class InputValidationError(ValueError):
    """Raised with all problems found in an SRFM input mapping."""

    def __init__(self, issues: Sequence[str]):
        """Initialize an error containing every detected input problem.

        Keeping all issues on the exception lets callers report the complete
        validation result instead of fixing one field at a time.

        Args:
            issues: Human-readable validation problems in discovery order.
        """
        self.issues = tuple(issues)
        detail = "\n".join(f"- {issue}" for issue in self.issues)
        super().__init__(f"Invalid SRFM inputs:\n{detail}")


@dataclass(frozen=True)
class FieldSpec:
    """Describe one field while retaining the driver's flat dictionary layout."""

    expected: tuple[type, ...] | None = None
    required: bool = False
    nullable: bool = False
    default: Any = _MISSING
    choices: frozenset[Any] | None = None
    item_choices: frozenset[Any] | None = None
    minimum: Real | None = None
    maximum: Real | None = None
    minimum_inclusive: bool = True
    maximum_inclusive: bool = True
    permitted: str | None = None


PATH_TYPES = (str, os.PathLike)
NUMBER_TYPES = (Real,)
INTEGER_TYPES = (int,)
MAPPING_TYPES = (Mapping,)


# Authoritative list from the SELECT CASE in RFM/source/drvflg_sub.f90.  CIA is
# a deprecated no-op, but the bundled RFM still accepts it and emits a warning,
# so Python validation must accept it as well.
RFM_FLAG_CODES: tuple[str, ...] = (
    "ABS",
    "AVG",
    "BBT",
    "BFX",
    "BIN",
    "C32",
    "C41",
    "C4C",
    "CHI",
    "CIA",
    "CLC",
    "COO",
    "CTM",
    "DBL",
    "FIN",
    "FLX",
    "FOV",
    "FVZ",
    "GEO",
    "GHZ",
    "GRA",
    "GRD",
    "HOM",
    "HYD",
    "ILS",
    "JAC",
    "JTP",
    "LAY",
    "LEV",
    "LIN",
    "LOS",
    "LUT",
    "MIX",
    "MTX",
    "NAD",
    "NEW",
    "NTE",
    "OBS",
    "OPT",
    "PRF",
    "PTH",
    "QAD",
    "RAD",
    "REJ",
    "REX",
    "RJT",
    "SFC",
    "SHH",
    "SHP",
    "SVD",
    "TAB",
    "TRA",
    "VRT",
    "VVW",
    "WID",
    "ZEN",
)
RFM_FLAG_CODES_SET = frozenset(RFM_FLAG_CODES)


SRFM_INPUT_SCHEMA: dict[str, FieldSpec] = {
    # General/provenance fields retained from existing driver tables.
    "fwd_model": FieldSpec((str,)),
    "instrument": FieldSpec((str,)),
    "plot_profiles": FieldSpec((bool,)),
    "iasi_spc_fldr": FieldSpec(PATH_TYPES, nullable=True),
    "iasi_fl": FieldSpec((str,), nullable=True),
    "px": FieldSpec(INTEGER_TYPES, nullable=True),
    "g_rtv": FieldSpec(nullable=True),
    "ils": FieldSpec(nullable=True),
    "nedt": FieldSpec(nullable=True),
    "sza_cos": FieldSpec(NUMBER_TYPES, nullable=True),
    "zen_cos": FieldSpec(NUMBER_TYPES, nullable=True),
    "zen_sec": FieldSpec(NUMBER_TYPES, nullable=True),
    # Output and run locations.
    "results_fldr": FieldSpec(PATH_TYPES, required=True),
    "base_plots": FieldSpec((bool,), required=True),
    "out_mode": FieldSpec(
        (str,), required=True, nullable=True, choices=frozenset({"txt", "netcdf"})
    ),
    "show_plots": FieldSpec((bool,), required=True),
    "out_fname": FieldSpec((str,), nullable=True),
    "convolve_iasi": FieldSpec((bool,), required=True),
    "iasi_ils": FieldSpec(PATH_TYPES, required=True, nullable=True),
    # Memory controls.
    "scattering_block_size": FieldSpec(INTEGER_TYPES, default=10000),
    "retain_phase_functions": FieldSpec((bool,), default=False),
    "retain_outputs": FieldSpec((list, tuple, set, frozenset), required=True),
    # Spectral grids.
    "fin_wvnmlo": FieldSpec(NUMBER_TYPES, required=True),
    "fin_wvnmhi": FieldSpec(NUMBER_TYPES, required=True),
    "fin_res": FieldSpec(NUMBER_TYPES, required=True),
    "spc_wvnmlo": FieldSpec(NUMBER_TYPES, required=True),
    "spc_wvnmhi": FieldSpec(NUMBER_TYPES, required=True),
    "spc_res": FieldSpec(NUMBER_TYPES, required=True),
    "spc_units": FieldSpec(
        (str,), required=True, choices=frozenset({"cm-1", "um", "nm"})
    ),
    # Output geometry
    "out_fmt": FieldSpec((str,), required=True, choices=frozenset({"altitude", "tau"})),
    "out": FieldSpec(required=True, nullable=True),
    "out_toa": FieldSpec((bool,), required=True),
    # RFM structures.
    "rfm_config": FieldSpec(MAPPING_TYPES, required=True),
    "driver_inputs": FieldSpec(MAPPING_TYPES, required=True),
    "levels": FieldSpec(nullable=True),
    # DISORT configuration.
    "fisot": FieldSpec(
        NUMBER_TYPES,
        required=True,
        minimum=0,
        permitted="a finite real number greater than or equal to 0",
    ),
    "albedo": FieldSpec(
        (Real, Mapping),
        required=True,
        permitted=(
            "a finite real number in [0, 1] or an in-memory/file-backed "
            "spectral-field mapping"
        ),
    ),
    "temis": FieldSpec(
        NUMBER_TYPES,
        required=True,
        minimum=0,
        maximum=1,
        permitted="a finite real number in the inclusive range [0, 1]",
    ),
    "earth_radius": FieldSpec(
        NUMBER_TYPES,
        minimum=0,
        minimum_inclusive=False,
        permitted="a finite real number greater than 0 km",
    ),
    "nmom": FieldSpec(
        INTEGER_TYPES,
        required=True,
        minimum=0,
        permitted="an integer greater than or equal to 0",
    ),
    "maxcmu": FieldSpec(
        INTEGER_TYPES,
        required=True,
        minimum=4,
        permitted="an even integer greater than or equal to 4",
    ),
    "maxumu": FieldSpec(
        INTEGER_TYPES,
        required=True,
        minimum=0,
        permitted=(
            "a non-negative integer that is positive unless onlyfl is True; "
            "when usrang is False and onlyfl is False it must be at least maxcmu"
        ),
    ),
    "maxphi": FieldSpec(
        INTEGER_TYPES,
        required=True,
        minimum=0,
        permitted="a non-negative integer that is positive unless onlyfl is True",
    ),
    "usrang": FieldSpec((bool,), required=True, permitted="True or False"),
    "usrtau": FieldSpec(
        (bool,), required=True, permitted="True for SRFM's user-output-level workflow"
    ),
    "ibcnd": FieldSpec(
        INTEGER_TYPES,
        required=True,
        choices=frozenset({0, 1}),
        permitted="0 (general case) or 1 (albedo/transmissivity case)",
    ),
    "onlyfl": FieldSpec(
        (bool,),
        required=True,
        permitted="True or False; it must be False when ibcnd is 1",
    ),
    "prnt": FieldSpec(
        (list,), required=True, permitted="a list containing exactly five booleans"
    ),
    "planck": FieldSpec((bool,), required=True, permitted="True or False"),
    "lamber": FieldSpec((bool,), required=True, permitted="True or False"),
    "deltamplus": FieldSpec((bool,), required=True, permitted="True or False"),
    "do_pseudo_sphere": FieldSpec(
        (bool,), required=True, permitted="True or False"
    ),
    "disort_precision": FieldSpec(
        (str,),
        required=True,
        choices=frozenset({"single", "double"}),
        permitted="'single' or 'double'",
    ),
    "header": FieldSpec(
        (str,),
        permitted="a string containing at most 127 characters, including an empty string",
    ),
    "adjust_maxcmu": FieldSpec((bool,), required=True, permitted="True or False"),
    "btemp": FieldSpec(
        NUMBER_TYPES,
        minimum=0,
        permitted="a finite real number greater than or equal to 0 K",
    ),
    "ttemp": FieldSpec(
        NUMBER_TYPES,
        minimum=0,
        permitted="a finite real number greater than or equal to 0 K",
    ),
    # Scattering and geometry.
    "scat_lyrs_inputs": FieldSpec(MAPPING_TYPES, nullable=True),
    "prescribed_lyrs_inputs": FieldSpec(MAPPING_TYPES, nullable=True),
    "gbc_lyrs_inputs": FieldSpec(MAPPING_TYPES, nullable=True),
    "solar_spectrum": FieldSpec(MAPPING_TYPES, nullable=True),
    "date": FieldSpec((dt.datetime, tuple)),
    "sun": FieldSpec((bool,), required=True, permitted="True or False"),
    "sza": FieldSpec(
        NUMBER_TYPES,
        required=True,
        minimum=0,
        maximum=180,
        permitted=(
            "a finite angle in [0, 180] degrees and, when sun is True, "
            "in [0, 90) degrees"
        ),
    ),
    "saa": FieldSpec(
        NUMBER_TYPES,
        required=True,
        minimum=0,
        maximum=360,
        permitted="a finite angle in the inclusive range [0, 360] degrees",
    ),
    "zen": FieldSpec(
        NUMBER_TYPES,
        required=True,
        minimum=0,
        maximum=180,
        permitted=(
            "a finite angle in [0, 180] degrees other than 90 degrees when "
            "angular intensities are requested; with ibcnd=1 and usrang=True, "
            "an angle in [0, 90) degrees"
        ),
    ),
    "azi": FieldSpec(
        NUMBER_TYPES,
        required=True,
        minimum=0,
        maximum=360,
        permitted="a finite angle in the inclusive range [0, 360] degrees",
    ),
}


# OXHARP passes a merged state/ancillary mapping to ``oxharp_main``.  The
# generic fields remain accepted for compatibility with existing driver tables,
# but OXHARP's derived geometry is the contract actually consumed by the runner.
OXHARP_INPUT_SCHEMA: dict[str, FieldSpec] = {
    **SRFM_INPUT_SCHEMA,
    "plot_profiles": FieldSpec((bool,)),
    "base_plots": FieldSpec((bool,)),
    "show_plots": FieldSpec((bool,)),
    "sza": FieldSpec(NUMBER_TYPES, minimum=0, maximum=180),
    "saa": FieldSpec(NUMBER_TYPES, minimum=0, maximum=360),
    "zen": FieldSpec(NUMBER_TYPES, minimum=0, maximum=180),
    "sza_cos": FieldSpec(
        NUMBER_TYPES,
        minimum=-1,
        maximum=1,
        permitted="a finite cosine in [-1, 1] and, when sun is True, in (0, 1]",
    ),
    "zen_cos": FieldSpec(
        NUMBER_TYPES,
        required=True,
        minimum=-1,
        maximum=1,
        permitted=(
            "a finite cosine in [-1, 1] other than 0 when angular intensities "
            "are requested; with ibcnd=1 and usrang=True, a cosine in (0, 1]"
        ),
    ),
    "zen_sec": FieldSpec(NUMBER_TYPES, required=True),
    "sza_deg": FieldSpec(NUMBER_TYPES),
    "sza_rad": FieldSpec(NUMBER_TYPES),
    "sza_sec": FieldSpec(NUMBER_TYPES),
    "zen_deg": FieldSpec(NUMBER_TYPES),
    "zen_rad": FieldSpec(NUMBER_TYPES),
    "atmosphere": FieldSpec(PATH_TYPES),
    "g_rtv": FieldSpec(MAPPING_TYPES, nullable=True),
    "maxcount": FieldSpec(INTEGER_TYPES),
    "epsilon": FieldSpec(NUMBER_TYPES),
    "gamma": FieldSpec(NUMBER_TYPES),
    "parallel": FieldSpec((bool,)),
    "eps": FieldSpec(NUMBER_TYPES),
    "max_workers": FieldSpec(INTEGER_TYPES, nullable=True),
}


# ``iasi_main`` reads observation geometry and profiles from a processed IASI
# pickle.  Consequently its file/pixel inputs are mandatory, whereas the
# generic runner's explicit geometry and optional-convolution switch are not.
IASI_INPUT_SCHEMA: dict[str, FieldSpec] = {
    **SRFM_INPUT_SCHEMA,
    "plot_profiles": FieldSpec((bool,), required=True),
    "iasi_spc_fldr": FieldSpec(PATH_TYPES, required=True),
    "iasi_fl": FieldSpec((str,), required=True),
    "px": FieldSpec(INTEGER_TYPES, required=True),
    "g_rtv": FieldSpec(MAPPING_TYPES, nullable=True),
    "ils": FieldSpec(PATH_TYPES, required=True),
    "nedt": FieldSpec(PATH_TYPES, required=True),
    "convolve_iasi": FieldSpec((bool,)),
    "iasi_ils": FieldSpec(PATH_TYPES, nullable=True),
    "sun": FieldSpec((bool,)),
    "sza": FieldSpec(NUMBER_TYPES),
    "saa": FieldSpec(NUMBER_TYPES),
    "zen": FieldSpec(NUMBER_TYPES),
    "azi": FieldSpec(NUMBER_TYPES),
}


RFM_CONFIG_SCHEMA: dict[str, FieldSpec] = {
    "output_mode": FieldSpec((str,), choices=frozenset({"files", "capture"})),
    "driver_path": FieldSpec(PATH_TYPES, nullable=True),
    "generate_driver": FieldSpec((bool,)),
    "verbose": FieldSpec((bool,)),
    "capture_files_content": FieldSpec((bool,)),
    "clean_before": FieldSpec((bool,)),
    "run_id": FieldSpec((str,), nullable=True),
    "patterns": FieldSpec(nullable=True),
    "optical_spectrum_index": FieldSpec(INTEGER_TYPES),
    "optical_match_tol": FieldSpec(NUMBER_TYPES),
    "optical_depth_format": FieldSpec(
        (str,), choices=frozenset({"compact", "dataframe"})
    ),
}


RFM_DRIVER_SCHEMA: dict[str, FieldSpec] = {
    "header": FieldSpec(required=True),
    "flags": FieldSpec(required=True, item_choices=RFM_FLAG_CODES_SET),
    "spectral": FieldSpec(required=True),
    "gases": FieldSpec(required=True),
    "atmosphere": FieldSpec(required=True),
    "tangent": FieldSpec(nullable=True),
    "tab_dimensions": FieldSpec(nullable=True),
    "cia": FieldSpec(nullable=True),
    "fin": FieldSpec(nullable=True),
    "fov": FieldSpec(nullable=True),
    "grd": FieldSpec(nullable=True),
    "hit": FieldSpec(nullable=True),
    "ils": FieldSpec(nullable=True),
    "jac": FieldSpec(nullable=True),
    "lev": FieldSpec(nullable=True),
    "lut": FieldSpec(nullable=True),
    "nte": FieldSpec(nullable=True),
    "obs": FieldSpec(nullable=True),
    "out": FieldSpec(nullable=True),
    "phy": FieldSpec(nullable=True),
    "rej": FieldSpec(nullable=True),
    "sfc": FieldSpec(nullable=True),
    "shp": FieldSpec(nullable=True),
    "svd": FieldSpec(nullable=True),
    "xsc": FieldSpec(nullable=True),
    "extra_sections": FieldSpec(MAPPING_TYPES, nullable=True),
}


LAYER_SCHEMA: dict[str, FieldSpec] = {
    "name": FieldSpec((str,), required=True),
    "low_spc": FieldSpec(NUMBER_TYPES, required=True),
    "upp_spc": FieldSpec(NUMBER_TYPES, required=True),
    "res": FieldSpec(NUMBER_TYPES, required=True),
    "spec_units": FieldSpec(
        (str,), required=True, choices=frozenset({"cm-1", "um", "nm"})
    ),
    "mass_loading": FieldSpec(NUMBER_TYPES, required=True, nullable=True),
    "n": FieldSpec(NUMBER_TYPES, required=True, nullable=True),
    "r": FieldSpec(NUMBER_TYPES, required=True),
    "s": FieldSpec(NUMBER_TYPES, required=True),
    "rho": FieldSpec((Real, str), required=True),
    "s_a_den": FieldSpec(NUMBER_TYPES, required=True, nullable=True),
    "v_den": FieldSpec(NUMBER_TYPES, required=True, nullable=True),
    "dist_type": FieldSpec(
        (str,), required=True, choices=frozenset({"log_normal", "gaussian"})
    ),
    "comp": FieldSpec((str,), required=True),
    "refractive_index": FieldSpec(nullable=True),
    "center_alt": FieldSpec(NUMBER_TYPES, nullable=True),
    "thick": FieldSpec(NUMBER_TYPES, nullable=True),
    "alt_upp": FieldSpec(NUMBER_TYPES, nullable=True),
    "alt_low": FieldSpec(NUMBER_TYPES, nullable=True),
    "radii": FieldSpec(INTEGER_TYPES, required=True),
    "eta": FieldSpec(NUMBER_TYPES, required=True),
    "phase_quad_N": FieldSpec(INTEGER_TYPES, required=True),
    "phase_quad_type": FieldSpec(
        (str,), required=True, choices=frozenset({"G", "R", "L", "T"})
    ),
    "radii_quad_type": FieldSpec(
        (str,), required=True, choices=frozenset({"G", "R", "L", "T"})
    ),
    "leg_coeffs": FieldSpec((bool,), required=True),
    "leg_coeffs_type": FieldSpec(
        (str,), required=True, choices=frozenset({"regular", "normalised"})
    ),
    "multiprocess": FieldSpec((bool,), required=True),
    "angle": FieldSpec(nullable=True),
}


GREY_BODY_LAYER_SCHEMA: dict[str, FieldSpec] = {
    "name": FieldSpec((str,), required=True),
    "low_spc": FieldSpec(NUMBER_TYPES, required=True),
    "upp_spc": FieldSpec(NUMBER_TYPES, required=True),
    "res": FieldSpec(NUMBER_TYPES, required=True),
    "spec_units": FieldSpec(
        (str,), required=True, choices=frozenset({"cm-1", "um", "nm"})
    ),
    "center_alt": FieldSpec(NUMBER_TYPES, nullable=True),
    "thick": FieldSpec(NUMBER_TYPES, nullable=True),
    "alt_upp": FieldSpec(NUMBER_TYPES, nullable=True),
    "alt_low": FieldSpec(NUMBER_TYPES, nullable=True),
    "emis": FieldSpec(NUMBER_TYPES, required=True),
    "inp_tau": FieldSpec(NUMBER_TYPES, required=True),
}


PRESCRIBED_LAYER_SCHEMA: dict[str, FieldSpec] = {
    "name": FieldSpec((str,), required=True),
    "alt_low": FieldSpec(NUMBER_TYPES, required=True),
    "alt_upp": FieldSpec(NUMBER_TYPES, required=True),
    "optical_depth": FieldSpec(MAPPING_TYPES, required=True),
    "ssalb": FieldSpec((Real, Mapping), required=True),
    "phase_function": FieldSpec(MAPPING_TYPES, nullable=True),
}


_SPECTRAL_FIELD_KEYS = {
    "grid",
    "values",
    "file",
    "grid_column",
    "value_column",
    "value_columns",
    "skiprows",
    "grid_units",
    "value_units",
}


def _is_sequence(value: Any) -> bool:
    """Return whether a value is a non-string sequence.

    Strings and bytes are excluded because schema fields such as levels and
    RFM sections require containers of distinct values.

    Args:
        value: Object to classify.

    Returns:
        ``True`` when the object is a sequence other than text or bytes.
    """
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes))


def _is_array_like(value: Any) -> bool:
    """Return whether a value is a supported sequence or NumPy array.

    The validation layer accepts ordinary Python sequences and NumPy arrays
    for numerical vector fields.

    Args:
        value: Object to classify.

    Returns:
        ``True`` when the value can represent a schema vector.
    """
    return _is_sequence(value) or isinstance(value, np.ndarray)


def _append_value_issue(
    issues: list[str], field_path: str, value: Any, permitted: str
) -> None:
    """Append a diagnostic containing the field, supplied value, and constraint."""
    issue = (
        f"{field_path}: current value {value!r}; "
        f"permitted values: {permitted}"
    )
    if issue not in issues:
        issues.append(issue)


def _check_mapping(
    value: Mapping[str, Any],
    schema: Mapping[str, FieldSpec],
    path: str,
    issues: list[str],
) -> None:
    """Validate mapping keys and scalar field specifications.

    Problems are appended to a shared list so nested and top-level failures
    can be reported together after all validation passes have completed.

    Args:
        value: Mapping being validated.
        schema: Field specifications allowed at this mapping level.
        path: Dotted prefix used in validation messages.
        issues: Mutable collection receiving detected problems.
    """
    unknown = sorted(set(value) - set(schema))
    for key in unknown:
        issues.append(f"{path}{key}: unknown field")

    for key, spec in schema.items():
        field_path = f"{path}{key}"
        if key not in value:
            if spec.required and spec.default is _MISSING:
                issues.append(f"{field_path}: required field is missing")
            continue
        item = value[key]
        if item is None:
            if not spec.nullable:
                if spec.permitted is None:
                    issues.append(f"{field_path}: may not be None")
                else:
                    _append_value_issue(issues, field_path, item, spec.permitted)
            continue
        if spec.expected is not None:
            valid_type = isinstance(item, spec.expected)
            if spec.expected == (bool,):
                valid_type = type(item) is bool
            elif any(
                isinstance(expected_type, type) and issubclass(expected_type, Real)
                for expected_type in spec.expected
            ):
                valid_type = valid_type and not isinstance(item, (bool, np.bool_))
            if not valid_type:
                if spec.permitted is None:
                    expected = ", ".join(t.__name__ for t in spec.expected)
                    issues.append(
                        f"{field_path}: expected {expected}, got {type(item).__name__}"
                    )
                else:
                    _append_value_issue(issues, field_path, item, spec.permitted)
                continue
            if isinstance(item, Real) and not np.isfinite(item):
                if spec.permitted is None:
                    issues.append(f"{field_path}: must be finite")
                else:
                    _append_value_issue(issues, field_path, item, spec.permitted)
                continue
        if spec.choices is not None and item not in spec.choices:
            if spec.permitted is None:
                choices = ", ".join(
                    repr(choice) for choice in sorted(spec.choices, key=str)
                )
                issues.append(f"{field_path}: expected one of {choices}, got {item!r}")
            else:
                _append_value_issue(issues, field_path, item, spec.permitted)
            continue
        below_minimum = spec.minimum is not None and (
            item < spec.minimum
            if spec.minimum_inclusive
            else item <= spec.minimum
        )
        above_maximum = spec.maximum is not None and (
            item > spec.maximum
            if spec.maximum_inclusive
            else item >= spec.maximum
        )
        if below_minimum or above_maximum:
            permitted = spec.permitted
            if permitted is None:
                lower = "-infinity" if spec.minimum is None else repr(spec.minimum)
                upper = "infinity" if spec.maximum is None else repr(spec.maximum)
                permitted = f"a value between {lower} and {upper}"
            _append_value_issue(issues, field_path, item, permitted)


def _validate_positive(mapping, keys, path, issues, *, allow_zero=False):
    """Validate selected numeric fields as positive or non-negative.

    Missing and non-numeric values are left to the structural type checks.

    Args:
        mapping: Input mapping containing the selected fields.
        keys: Field names whose numerical signs should be checked.
        path: Dotted prefix used in validation messages.
        issues: Mutable collection receiving detected problems.
        allow_zero: Accept zero and reject only negative values when true.
    """
    for key in keys:
        value = mapping.get(key)
        if isinstance(value, Real) and not isinstance(value, (bool, np.bool_)):
            invalid = value < 0 if allow_zero else value <= 0
            if invalid:
                relation = "non-negative" if allow_zero else "greater than zero"
                issues.append(f"{path}{key}: must be {relation}")


def _validate_disort_inputs(
    mapping: Mapping[str, Any], runner: str, issues: list[str]
) -> None:
    """Validate user-controlled values that reach native DISORT arguments.

    The constraints mirror fatal checks in ``DISORT.f`` plus the stricter
    requirements introduced by SRFM's conversion from degrees to angle cosines.
    Every diagnostic uses the same actionable format: input name, supplied value,
    and permitted values.

    Args:
        mapping: Normalized top-level input mapping.
        runner: Runner identifier used to select its geometry representation.
        issues: Mutable collection receiving detected problems.
    """
    maxcmu = mapping.get("maxcmu")
    if type(maxcmu) is int and maxcmu % 2:
        _append_value_issue(
            issues, "maxcmu", maxcmu, SRFM_INPUT_SCHEMA["maxcmu"].permitted
        )

    onlyfl = mapping.get("onlyfl")
    for key in ("maxumu", "maxphi"):
        value = mapping.get(key)
        if onlyfl is False and type(value) is int and value == 0:
            _append_value_issue(
                issues, key, value, SRFM_INPUT_SCHEMA[key].permitted
            )

    maxumu = mapping.get("maxumu")
    if (
        mapping.get("usrang") is False
        and onlyfl is False
        and type(maxumu) is int
        and type(maxcmu) is int
        and maxumu < maxcmu
    ):
        _append_value_issue(
            issues, "maxumu", maxumu, SRFM_INPUT_SCHEMA["maxumu"].permitted
        )

    if mapping.get("usrtau") is False:
        _append_value_issue(
            issues,
            "usrtau",
            False,
            SRFM_INPUT_SCHEMA["usrtau"].permitted,
        )

    if mapping.get("ibcnd") == 1 and onlyfl is True:
        _append_value_issue(
            issues,
            "onlyfl",
            True,
            SRFM_INPUT_SCHEMA["onlyfl"].permitted,
        )

    output = mapping.get("out")
    if mapping.get("ibcnd") == 1 and isinstance(output, list) and all(
        isinstance(value, Real) for value in output
    ):
        output_count = len(output)
        if mapping.get("out_toa") is True:
            if mapping.get("out_fmt") == "tau":
                contains_toa = any(np.isclose(value, 0.0) for value in output)
            else:
                levels = mapping.get("levels")
                contains_toa = False
                if _is_array_like(levels) and len(levels):
                    top_altitude = levels[-1]
                    if isinstance(top_altitude, Real):
                        contains_toa = any(
                            np.isclose(value, top_altitude) for value in output
                        )
            if not contains_toa:
                output_count += 1
        if output_count < 2:
            _append_value_issue(
                issues,
                "out",
                output,
                (
                    "at least two resolved output levels when ibcnd is 1 "
                    "(native DISORT requires effective MAXULV >= 2)"
                ),
            )

    prnt = mapping.get("prnt")
    if "prnt" in mapping and (
        not isinstance(prnt, list)
        or len(prnt) != 5
        or not all(type(value) is bool for value in prnt)
    ):
        _append_value_issue(
            issues, "prnt", prnt, SRFM_INPUT_SCHEMA["prnt"].permitted
        )

    header = mapping.get("header")
    if isinstance(header, str) and len(header) > 127:
        _append_value_issue(
            issues, "header", header, SRFM_INPUT_SCHEMA["header"].permitted
        )

    angular_intensities = mapping.get("usrang") is True and onlyfl is False
    if runner == "srfm":
        solar_zenith = mapping.get("sza")
        if (
            mapping.get("sun") is True
            and isinstance(solar_zenith, Real)
            and not isinstance(solar_zenith, (bool, np.bool_))
            and solar_zenith >= 90
        ):
            _append_value_issue(
                issues,
                "sza",
                solar_zenith,
                SRFM_INPUT_SCHEMA["sza"].permitted,
            )

        viewing_zenith = mapping.get("zen")
        invalid_view = (
            angular_intensities
            and isinstance(viewing_zenith, Real)
            and np.isclose(viewing_zenith, 90)
        )
        invalid_ibcnd_view = (
            mapping.get("ibcnd") == 1
            and mapping.get("usrang") is True
            and isinstance(viewing_zenith, Real)
            and viewing_zenith >= 90
        )
        if invalid_view or invalid_ibcnd_view:
            _append_value_issue(
                issues,
                "zen",
                viewing_zenith,
                SRFM_INPUT_SCHEMA["zen"].permitted,
            )
    elif runner == "oxharp":
        solar_cosine = mapping.get("sza_cos")
        if (
            mapping.get("sun") is True
            and isinstance(solar_cosine, Real)
            and solar_cosine <= 0
        ):
            _append_value_issue(
                issues,
                "sza_cos",
                solar_cosine,
                OXHARP_INPUT_SCHEMA["sza_cos"].permitted,
            )

        viewing_cosine = mapping.get("zen_cos")
        invalid_view = angular_intensities and viewing_cosine == 0
        invalid_ibcnd_view = (
            mapping.get("ibcnd") == 1
            and mapping.get("usrang") is True
            and isinstance(viewing_cosine, Real)
            and viewing_cosine <= 0
        )
        if invalid_view or invalid_ibcnd_view:
            _append_value_issue(
                issues,
                "zen_cos",
                viewing_cosine,
                OXHARP_INPUT_SCHEMA["zen_cos"].permitted,
            )


def _validate_rfm_config(config: Any, issues: list[str]) -> None:
    """Validate nested RFM execution configuration.

    Structural type errors are recorded by the top-level schema, so this
    helper only descends into mappings.

    Args:
        config: Candidate RFM configuration mapping.
        issues: Mutable collection receiving detected problems.
    """
    if not isinstance(config, Mapping):
        return
    _check_mapping(config, RFM_CONFIG_SCHEMA, "rfm_config.", issues)
    patterns = config.get("patterns")
    if patterns is not None and not _is_sequence(patterns):
        issues.append("rfm_config.patterns: expected a non-string sequence or None")
    _validate_positive(
        config,
        ("optical_spectrum_index", "optical_match_tol"),
        "rfm_config.",
        issues,
    )


def _validate_driver(
    driver: Any, issues: list[str], *, minimum_atmospheres: int = 1
) -> None:
    """Validate nested RFM driver sections and required capture flags.

    The IASI runner needs two pre-existing atmosphere entries because it reads
    the profile at index one before appending its generated profile.

    Args:
        driver: Candidate RFM driver-section mapping.
        issues: Mutable collection receiving detected problems.
        minimum_atmospheres: Minimum number of atmosphere section entries.
    """
    if not isinstance(driver, Mapping):
        return
    _check_mapping(driver, RFM_DRIVER_SCHEMA, "driver_inputs.", issues)
    header = driver.get("header")
    if header is not None and not (
        isinstance(header, str)
        or (_is_sequence(header) and all(isinstance(line, str) for line in header))
    ):
        issues.append("driver_inputs.header: must be a string or sequence of strings")
    for key in ("flags", "spectral", "gases", "atmosphere"):
        value = driver.get(key)
        if key in driver and (not _is_sequence(value) or len(value) == 0):
            issues.append(f"driver_inputs.{key}: must be a non-empty sequence")
    atmosphere = driver.get("atmosphere")
    if _is_sequence(atmosphere) and len(atmosphere) < minimum_atmospheres:
        issues.append(
            "driver_inputs.atmosphere: must contain at least "
            f"{minimum_atmospheres} entries"
        )
    flags = driver.get("flags")
    if _is_sequence(flags):
        if not all(isinstance(flag, str) for flag in flags):
            issues.append("driver_inputs.flags: every flag must be a string")
        else:
            accepted_flags = RFM_DRIVER_SCHEMA["flags"].item_choices
            assert accepted_flags is not None
            unknown = sorted({flag.upper() for flag in flags} - accepted_flags)
            if unknown:
                issues.append(
                    "driver_inputs.flags: unknown RFM flag code(s): "
                    + ", ".join(unknown)
                )
            required = {"OPT", "LEV"} - {flag.upper() for flag in flags}
            if required:
                issues.append(
                    "driver_inputs.flags: run_srfm requires "
                    + " and ".join(sorted(required))
                    + " for optical-depth capture"
                )


def _validate_spectral_field_mapping(
    value: Any,
    path: str,
    issues: list[str],
    *,
    matrix_values: bool = False,
    minimum: float | None = None,
    maximum: float | None = None,
    require_value_units: bool = False,
    extra_keys: set[str] | None = None,
) -> None:
    """Validate one in-memory or file-backed spectral-field mapping.

    File contents are deliberately loaded by the runner preflight, where coverage
    can be checked against the computational grid before model side effects.

    Args:
        value: Candidate spectral-field mapping.
        path: Dotted path used in diagnostics.
        issues: Mutable collection receiving validation failures.
        matrix_values: Require a spectral-by-component value matrix.
        minimum: Optional inclusive lower value bound.
        maximum: Optional inclusive upper value bound.
        require_value_units: Require one supported solar-density unit spelling.
        extra_keys: Representation-specific keys accepted in addition to the
            common spectral-field keys.
    """
    if not isinstance(value, Mapping):
        issues.append(f"{path}: expected a mapping")
        return
    allowed = _SPECTRAL_FIELD_KEYS | (extra_keys or set())
    for unknown in sorted(set(value) - allowed):
        issues.append(f"{path}.{unknown}: unknown field")

    has_file = "file" in value
    has_grid = "grid" in value
    has_values = "values" in value
    if has_file and (has_grid or has_values):
        issues.append(f"{path}: file and grid/values representations are mutually exclusive")
    elif not has_file and not (has_grid and has_values):
        issues.append(f"{path}: provide either file or both grid and values")
    elif has_grid != has_values:
        issues.append(f"{path}: grid and values must be supplied together")

    grid_units = value.get("grid_units")
    if grid_units not in GRID_UNITS:
        issues.append(
            f"{path}.grid_units: expected one of "
            + ", ".join(repr(unit) for unit in sorted(GRID_UNITS))
        )
    if require_value_units:
        value_units = value.get("value_units")
        if value_units not in SOLAR_VALUE_UNITS:
            issues.append(
                f"{path}.value_units: expected one of "
                + ", ".join(repr(unit) for unit in sorted(SOLAR_VALUE_UNITS))
            )
    elif "value_units" in value:
        issues.append(f"{path}.value_units: only solar_spectrum accepts value units")

    if has_file:
        if not isinstance(value.get("file"), PATH_TYPES):
            issues.append(f"{path}.file: expected str or path-like value")
        grid_column = value.get("grid_column", 0)
        if type(grid_column) is not int or grid_column < 0:
            issues.append(f"{path}.grid_column: must be a non-negative integer")
        skiprows = value.get("skiprows", 0)
        if type(skiprows) is not int or skiprows < 0:
            issues.append(f"{path}.skiprows: must be a non-negative integer")
        if matrix_values:
            columns = value.get("value_columns")
            if isinstance(columns, np.ndarray):
                columns = columns.tolist()
            if (
                not _is_sequence(columns)
                or len(columns) == 0
                or any(type(column) is not int or column < 0 for column in columns)
            ):
                issues.append(
                    f"{path}.value_columns: must be a non-empty sequence of "
                    "non-negative integers"
                )
            elif len(set(columns)) != len(columns):
                issues.append(f"{path}.value_columns: values must be distinct")
            if "value_column" in value:
                issues.append(f"{path}.value_column: use value_columns for moment data")
        else:
            value_column = value.get("value_column", 1)
            if type(value_column) is not int or value_column < 0:
                issues.append(f"{path}.value_column: must be a non-negative integer")
            if "value_columns" in value:
                issues.append(f"{path}.value_columns: use value_column for scalar data")
        return

    if not (has_grid and has_values):
        return
    try:
        grid = np.asarray(value["grid"], dtype=float)
        values = np.asarray(value["values"], dtype=float)
    except (TypeError, ValueError):
        issues.append(f"{path}: grid and values must be numeric arrays")
        return
    if grid.ndim != 1 or grid.size < 2 or not np.all(np.isfinite(grid)):
        issues.append(f"{path}.grid: must be a finite 1D array with at least two points")
    elif np.any(grid <= 0) or not (
        np.all(np.diff(grid) > 0) or np.all(np.diff(grid) < 0)
    ):
        issues.append(f"{path}.grid: values must be positive and strictly monotonic")
    if matrix_values:
        valid_shape = values.ndim == 2 and values.shape[0] == grid.size and values.shape[1] > 0
        expected = "shape (spectral_points, moments)"
    else:
        valid_shape = values.ndim == 1 and values.shape == grid.shape
        expected = "the same one-dimensional shape as grid"
    if not valid_shape:
        issues.append(f"{path}.values: must have {expected}")
        return
    if not np.all(np.isfinite(values)):
        issues.append(f"{path}.values: all values must be finite")
    if minimum is not None and np.any(values < minimum):
        issues.append(f"{path}.values: values must be greater than or equal to {minimum:g}")
    if maximum is not None and np.any(values > maximum):
        issues.append(f"{path}.values: values must be less than or equal to {maximum:g}")


def _validate_boundary_spectral_inputs(mapping: Mapping[str, Any], issues: list[str]) -> None:
    """Validate scalar/spectral albedo and optional custom solar irradiance."""
    albedo = mapping.get("albedo")
    if isinstance(albedo, Real) and not isinstance(albedo, (bool, np.bool_)):
        if not np.isfinite(albedo) or not 0 <= albedo <= 1:
            _append_value_issue(
                issues,
                "albedo",
                albedo,
                "a finite scalar in [0, 1] or a spectral-field mapping",
            )
    elif isinstance(albedo, Mapping):
        _validate_spectral_field_mapping(
            albedo, "albedo", issues, minimum=0, maximum=1
        )

    solar_spectrum = mapping.get("solar_spectrum")
    if solar_spectrum is not None:
        _validate_spectral_field_mapping(
            solar_spectrum,
            "solar_spectrum",
            issues,
            minimum=0,
            require_value_units=True,
        )


def _validate_prescribed_layers(layers: Any, issues: list[str]) -> None:
    """Validate prescribed-layer structure and all in-memory numerical fields."""
    if layers is None or not isinstance(layers, Mapping):
        return
    for layer_name, prescribed in layers.items():
        path = f"prescribed_lyrs_inputs.{layer_name}."
        if not isinstance(layer_name, str):
            issues.append(
                f"prescribed_lyrs_inputs: layer name {layer_name!r} must be a string"
            )
            continue
        if not isinstance(prescribed, Mapping):
            issues.append(f"prescribed_lyrs_inputs.{layer_name}: expected a mapping")
            continue
        _check_mapping(prescribed, PRESCRIBED_LAYER_SCHEMA, path, issues)
        configured_name = prescribed.get("name")
        if isinstance(configured_name, str) and configured_name != layer_name:
            issues.append(f"{path}name: must match the containing layer name {layer_name!r}")
        alt_low = prescribed.get("alt_low")
        alt_upp = prescribed.get("alt_upp")
        if isinstance(alt_low, Real) and isinstance(alt_upp, Real) and alt_low >= alt_upp:
            issues.append(f"{path}alt_low: must be less than alt_upp")

        optical_depth = prescribed.get("optical_depth")
        optical_path = f"{path}optical_depth"
        if isinstance(optical_depth, Mapping):
            optical_type = optical_depth.get("type")
            if optical_type == "angstrom":
                allowed = {
                    "type",
                    "reference_value",
                    "reference_wavelength_um",
                    "angstrom_exponent",
                }
                for unknown in sorted(set(optical_depth) - allowed):
                    issues.append(f"{optical_path}.{unknown}: unknown field")
                for required in allowed - {"type"}:
                    if required not in optical_depth:
                        issues.append(f"{optical_path}.{required}: required field is missing")
                reference_value = optical_depth.get("reference_value")
                reference_wavelength = optical_depth.get("reference_wavelength_um")
                exponent = optical_depth.get("angstrom_exponent")
                for key, candidate in (
                    ("reference_value", reference_value),
                    ("reference_wavelength_um", reference_wavelength),
                    ("angstrom_exponent", exponent),
                ):
                    if (
                        not isinstance(candidate, Real)
                        or isinstance(candidate, (bool, np.bool_))
                        or not np.isfinite(candidate)
                    ):
                        issues.append(f"{optical_path}.{key}: must be a finite real number")
                if isinstance(reference_value, Real) and reference_value < 0:
                    issues.append(f"{optical_path}.reference_value: must be non-negative")
                if isinstance(reference_wavelength, Real) and reference_wavelength <= 0:
                    issues.append(
                        f"{optical_path}.reference_wavelength_um: must be greater than zero"
                    )
            elif optical_type == "tabulated":
                _validate_spectral_field_mapping(
                    optical_depth,
                    optical_path,
                    issues,
                    minimum=0,
                    extra_keys={"type"},
                )
            else:
                issues.append(f"{optical_path}.type: expected 'angstrom' or 'tabulated'")

        single_scattering_albedo = prescribed.get("ssalb")
        if isinstance(single_scattering_albedo, Real) and not isinstance(
            single_scattering_albedo, (bool, np.bool_)
        ):
            if not np.isfinite(single_scattering_albedo) or not 0 <= single_scattering_albedo <= 1:
                issues.append(f"{path}ssalb: scalar value must be finite and in [0, 1]")
        elif isinstance(single_scattering_albedo, Mapping):
            _validate_spectral_field_mapping(
                single_scattering_albedo,
                f"{path}ssalb",
                issues,
                minimum=0,
                maximum=1,
            )
        elif "ssalb" in prescribed:
            issues.append(
                f"{path}ssalb: expected a scalar or spectral-field mapping; "
                "bare arrays are not supported"
            )

        phase_function = prescribed.get("phase_function")
        if isinstance(phase_function, Mapping):
            phase_type = phase_function.get("type")
            phase_path = f"{path}phase_function"
            if phase_type == "henyey_greenstein":
                for unknown in sorted(set(phase_function) - {"type", "asymmetry"}):
                    issues.append(f"{phase_path}.{unknown}: unknown field")
                if "asymmetry" not in phase_function:
                    issues.append(f"{phase_path}.asymmetry: required field is missing")
                asymmetry = phase_function.get("asymmetry")
                if isinstance(asymmetry, Real) and not isinstance(asymmetry, (bool, np.bool_)):
                    if not np.isfinite(asymmetry) or not -1 <= asymmetry <= 1:
                        issues.append(
                            f"{phase_path}.asymmetry: scalar must be finite and in [-1, 1]"
                        )
                elif isinstance(asymmetry, Mapping):
                    _validate_spectral_field_mapping(
                        asymmetry,
                        f"{phase_path}.asymmetry",
                        issues,
                        minimum=-1,
                        maximum=1,
                    )
                elif "asymmetry" in phase_function:
                    issues.append(
                        f"{phase_path}.asymmetry: expected a scalar or spectral-field mapping"
                    )
            elif phase_type == "legendre_moments":
                if phase_function.get("convention") != "normalised":
                    issues.append(f"{phase_path}.convention: expected 'normalised'")
                _validate_spectral_field_mapping(
                    phase_function,
                    phase_path,
                    issues,
                    matrix_values=True,
                    extra_keys={"type", "convention"},
                )
            else:
                issues.append(
                    f"{phase_path}.type: expected 'henyey_greenstein' or 'legendre_moments'"
                )


def _validate_layers(layers: Any, issues: list[str]) -> None:
    """Validate every configured scattering layer and cross-field constraint.

    Layer names are included in dotted error paths so failures remain
    actionable when a run contains several aerosol or cloud layers.

    Args:
        layers: Candidate mapping of layer names to layer input mappings.
        issues: Mutable collection receiving detected problems.
    """
    if layers is None:
        return
    if not isinstance(layers, Mapping):
        return
    for layer_name, layer in layers.items():
        path = f"scat_lyrs_inputs.{layer_name}."
        if not isinstance(layer_name, str):
            issues.append(
                f"scat_lyrs_inputs: layer name {layer_name!r} must be a string"
            )
            continue
        if not isinstance(layer, Mapping):
            issues.append(f"scat_lyrs_inputs.{layer_name}: expected a mapping")
            continue
        _check_mapping(layer, LAYER_SCHEMA, path, issues)
        configured_name = layer.get("name")
        if isinstance(configured_name, str) and configured_name != layer_name:
            issues.append(
                f"{path}name: must match the containing layer name {layer_name!r}"
            )
        _validate_positive(
            layer,
            ("res", "r", "s", "rho", "radii", "phase_quad_N"),
            path,
            issues,
        )
        eta = layer.get("eta")
        if isinstance(eta, Real) and not 0 < eta < 1:
            issues.append(f"{path}eta: must satisfy 0 < eta < 1")
        low, high = layer.get("low_spc"), layer.get("upp_spc")
        if isinstance(low, Real) and isinstance(high, Real) and low >= high:
            issues.append(f"{path}low_spc: must be less than upp_spc")
        for key in ("low_spc", "upp_spc"):
            value = layer.get(key)
            if isinstance(value, Real) and value < 0:
                issues.append(f"{path}{key}: must be non-negative")
        spread = layer.get("s")
        if (
            layer.get("dist_type") == "log_normal"
            and isinstance(spread, Real)
            and spread <= 1
        ):
            issues.append(f"{path}s: must be greater than 1 for log_normal")
        if layer.get("dist_type") == "gaussian":
            issues.append(
                f"{path}dist_type: gaussian is recognized but not implemented "
                "for a complete MieLayer calculation"
            )
        density = layer.get("rho")
        if isinstance(density, str) and density not in {
            "pumice",
            "glass",
            "mineral",
            "rock",
        }:
            issues.append(
                f"{path}rho: named density must be pumice, glass, mineral, or rock"
            )
        for key in ("mass_loading", "n", "s_a_den", "v_den"):
            value = layer.get(key)
            if isinstance(value, Real) and value < 0:
                issues.append(f"{path}{key}: must be non-negative")
        if all(
            layer.get(key) is None for key in ("mass_loading", "n", "s_a_den", "v_den")
        ):
            issues.append(
                f"{path}mass_loading: one of mass_loading, n, s_a_den, or v_den is required"
            )
        centre_extent = (
            layer.get("center_alt") is not None and layer.get("thick") is not None
        )
        bound_extent = (
            layer.get("alt_low") is not None and layer.get("alt_upp") is not None
        )
        if not centre_extent and not bound_extent:
            issues.append(
                f"{path}center_alt: provide center_alt/thick or alt_low/alt_upp"
            )
        thickness = layer.get("thick")
        if thickness is not None and isinstance(thickness, Real) and thickness <= 0:
            issues.append(f"{path}thick: must be greater than zero")
        alt_low, alt_upp = layer.get("alt_low"), layer.get("alt_upp")
        if (
            isinstance(alt_low, Real)
            and isinstance(alt_upp, Real)
            and alt_low >= alt_upp
        ):
            issues.append(f"{path}alt_low: must be less than alt_upp")
        if (
            centre_extent
            and bound_extent
            and all(
                isinstance(layer.get(key), Real)
                for key in ("center_alt", "thick", "alt_low", "alt_upp")
            )
        ):
            expected_low = round(layer["center_alt"] - layer["thick"] / 2, 3)
            expected_upp = round(layer["center_alt"] + layer["thick"] / 2, 3)
            if not (
                np.isclose(layer["alt_low"], expected_low)
                and np.isclose(layer["alt_upp"], expected_upp)
            ):
                issues.append(
                    f"{path}alt_low: explicit bounds must match center_alt/thick"
                )
        if layer.get("comp") == "ri" and layer.get("refractive_index") is None:
            issues.append(f"{path}refractive_index: required when comp is 'ri'")


def _validate_grey_body_layers(layers: Any, issues: list[str]) -> None:
    """Validate every configured non-scattering grey-body cloud layer.

    Args:
        layers: Candidate mapping of layer names to grey-body input mappings.
        issues: Mutable collection receiving detected problems.
    """
    if layers is None or not isinstance(layers, Mapping):
        return

    for layer_name, layer in layers.items():
        path = f"gbc_lyrs_inputs.{layer_name}."
        if not isinstance(layer_name, str):
            issues.append(
                f"gbc_lyrs_inputs: layer name {layer_name!r} must be a string"
            )
            continue
        if not isinstance(layer, Mapping):
            issues.append(f"gbc_lyrs_inputs.{layer_name}: expected a mapping")
            continue

        _check_mapping(layer, GREY_BODY_LAYER_SCHEMA, path, issues)
        _validate_positive(layer, ("res",), path, issues)

        configured_name = layer.get("name")
        if isinstance(configured_name, str) and configured_name != layer_name:
            issues.append(
                f"{path}name: must match the containing layer name {layer_name!r}"
            )

        low, high = layer.get("low_spc"), layer.get("upp_spc")
        if isinstance(low, Real) and isinstance(high, Real) and low >= high:
            issues.append(f"{path}low_spc: must be less than upp_spc")
        for key in ("low_spc", "upp_spc", "inp_tau"):
            value = layer.get(key)
            if isinstance(value, Real) and value < 0:
                issues.append(f"{path}{key}: must be non-negative")

        emissivity = layer.get("emis")
        if isinstance(emissivity, Real) and not 0 <= emissivity <= 1:
            issues.append(f"{path}emis: must be between 0 and 1")

        centre_extent = (
            layer.get("center_alt") is not None and layer.get("thick") is not None
        )
        bound_extent = (
            layer.get("alt_low") is not None and layer.get("alt_upp") is not None
        )
        if not centre_extent and not bound_extent:
            issues.append(
                f"{path}center_alt: provide center_alt/thick or alt_low/alt_upp"
            )

        thickness = layer.get("thick")
        if isinstance(thickness, Real) and thickness < 0.002:
            issues.append(f"{path}thick: must be at least 0.002 km")

        alt_low, alt_upp = layer.get("alt_low"), layer.get("alt_upp")
        if (
            isinstance(alt_low, Real)
            and isinstance(alt_upp, Real)
            and alt_low >= alt_upp
        ):
            issues.append(f"{path}alt_low: must be less than alt_upp")

        if (
            centre_extent
            and bound_extent
            and all(
                isinstance(layer.get(key), Real)
                for key in ("center_alt", "thick", "alt_low", "alt_upp")
            )
        ):
            expected_low = round(layer["center_alt"] - layer["thick"] / 2, 3)
            expected_upp = round(layer["center_alt"] + layer["thick"] / 2, 3)
            if not (
                np.isclose(layer["alt_low"], expected_low)
                and np.isclose(layer["alt_upp"], expected_upp)
            ):
                issues.append(
                    f"{path}alt_low: explicit bounds must match center_alt/thick"
                )


def _configured_layer_bounds(layer: Mapping[str, Any]) -> tuple[float, float] | None:
    """Return effective layer bounds when they can be resolved safely."""
    alt_low, alt_upp = layer.get("alt_low"), layer.get("alt_upp")
    if (
        isinstance(alt_low, Real)
        and isinstance(alt_upp, Real)
        and np.isfinite(alt_low)
        and np.isfinite(alt_upp)
        and alt_low < alt_upp
    ):
        return float(alt_low), float(alt_upp)
    center, thickness = layer.get("center_alt"), layer.get("thick")
    if (
        isinstance(center, Real)
        and isinstance(thickness, Real)
        and np.isfinite(center)
        and np.isfinite(thickness)
        and thickness > 0
    ):
        return (
            float(round(center - thickness / 2, 3)),
            float(round(center + thickness / 2, 3)),
        )
    return None


def _validate_optical_layer_structure(mapping: Mapping[str, Any], issues: list[str]) -> None:
    """Reject duplicate names, overlaps, and shared optical-layer boundaries."""
    groups = (
        ("scat_lyrs_inputs", mapping.get("scat_lyrs_inputs")),
        ("prescribed_lyrs_inputs", mapping.get("prescribed_lyrs_inputs")),
        ("gbc_lyrs_inputs", mapping.get("gbc_lyrs_inputs")),
    )
    ownership: dict[str, str] = {}
    intervals: list[tuple[float, float, str, str]] = []
    for group_name, layers in groups:
        if not isinstance(layers, Mapping):
            continue
        for layer_name, configured in layers.items():
            if not isinstance(layer_name, str) or not isinstance(configured, Mapping):
                continue
            previous_group = ownership.get(layer_name)
            if previous_group is not None:
                issues.append(
                    f"{group_name}.{layer_name}: layer name is already used by {previous_group}"
                )
            else:
                ownership[layer_name] = group_name
            bounds = _configured_layer_bounds(configured)
            if bounds is not None:
                intervals.append((*bounds, group_name, layer_name))

    intervals.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
    for first_index, first in enumerate(intervals):
        first_low, first_upp, first_group, first_name = first
        for second in intervals[first_index + 1 :]:
            second_low, second_upp, second_group, second_name = second
            if second_low > first_upp:
                break
            issues.append(
                f"{second_group}.{second_name}: optical layer [{second_low:g}, "
                f"{second_upp:g}] km overlaps or shares a boundary with "
                f"{first_group}.{first_name} [{first_low:g}, {first_upp:g}] km"
            )


def _validate_inputs(
    values: Mapping[str, Any],
    schema: Mapping[str, FieldSpec],
    *,
    runner: str,
) -> dict[str, Any]:
    """Validate inputs against a runner-specific schema.

    Shared spectral, DISORT, output, RFM, and layer constraints are applied to
    all runners before their distinct geometry and observation rules.

    Args:
        values: Flat input mapping supplied to a runner.
        schema: Top-level field contract for that runner.
        runner: Runner identifier: ``srfm``, ``oxharp``, or ``iasi``.

    Returns:
        A shallow normalized copy of the supplied mapping.

    Raises:
        InputValidationError: If the mapping violates any schema constraint.
    """
    if not isinstance(values, Mapping):
        raise InputValidationError(
            [f"inputs: expected a mapping, got {type(values).__name__}"]
        )

    normalized = dict(values)
    for key, spec in schema.items():
        if key not in normalized and spec.default is not _MISSING:
            normalized[key] = spec.default

    issues: list[str] = []
    _check_mapping(normalized, schema, "", issues)

    for prefix in ("fin", "spc"):
        low = normalized.get(f"{prefix}_wvnmlo")
        high = normalized.get(f"{prefix}_wvnmhi")
        resolution = normalized.get(f"{prefix}_res")
        if isinstance(low, Real) and isinstance(high, Real) and low >= high:
            issues.append(f"{prefix}_wvnmlo: must be less than {prefix}_wvnmhi")
        if isinstance(resolution, Real) and resolution <= 0:
            issues.append(f"{prefix}_res: must be greater than zero")

    date = normalized.get("date")
    if isinstance(date, tuple) and (
        len(date) != 3 or not all(isinstance(item, Integral) for item in date)
    ):
        issues.append("date: tuple form must contain exactly (year, month, day)")
    elif isinstance(date, tuple):
        try:
            normalized_date = tuple(int(item) for item in date)
            dt.datetime(*normalized_date)
        except (TypeError, ValueError) as exc:
            issues.append(f"date: {exc}")
        else:
            normalized["date"] = normalized_date

    levels = normalized.get("levels")
    if levels is not None:
        if not _is_array_like(levels) or len(levels) < 2:
            issues.append("levels: must be a sequence containing at least two levels")
        else:
            try:
                level_array = np.asarray(levels, dtype=float)
            except (TypeError, ValueError):
                issues.append("levels: every level must be numeric")
            else:
                if level_array.ndim != 1 or not np.all(np.isfinite(level_array)):
                    issues.append("levels: must be a finite one-dimensional sequence")
                elif not np.all(np.diff(level_array) > 0):
                    issues.append("levels: values must be strictly increasing")

    output = normalized.get("out")
    if "out" in normalized and output is not None:
        if isinstance(output, Real) and not isinstance(output, (bool, np.bool_)):
            output_values = [output]
        elif _is_array_like(output):
            try:
                output_array = np.asarray(output)
            except (TypeError, ValueError):
                output_values = None
            else:
                output_values = output_array.tolist() if output_array.ndim == 1 else None
        else:
            output_values = None

        if not output_values and normalized.get("out_toa") is not True:
            issues.append(
                "out: must be a non-empty numeric scalar or one-dimensional sequence"
            )
        elif not all(
            isinstance(value, Real)
            and not isinstance(value, (bool, np.bool_))
            and np.isfinite(value)
            and value >= 0
            for value in output_values
        ):
            issues.append("out: values must be finite non-negative numbers")
        else:
            normalized["out"] = list(output_values)
    elif "out" in normalized and output is None:
        if normalized.get("out_toa") is True:
            normalized["out"] = []
        else:
            issues.append(
                "out: must be specified unless out_toa is True"
            )

    _validate_positive(normalized, ("scattering_block_size",), "", issues)
    _validate_disort_inputs(normalized, runner, issues)
    _validate_boundary_spectral_inputs(normalized, issues)

    if (
        runner != "iasi"
        and normalized.get("convolve_iasi") is True
        and not normalized.get("iasi_ils")
    ):
        issues.append("iasi_ils: required when convolve_iasi is True")
    retained_outputs = normalized.get("retain_outputs")
    if isinstance(retained_outputs, (list, tuple, set, frozenset)):
        allowed_outputs = {
            "rad",
            "radiance",
            "bbt",
            "rfldir",
            "rfldn",
            "flup",
            "dfdt",
            "uavg",
            "uu",
            "albmed",
            "trnmed",
        }
        invalid_outputs = sorted(
            {
                repr(name)
                for name in retained_outputs
                if not isinstance(name, str) or name not in allowed_outputs
            }
        )
        if invalid_outputs:
            issues.append(
                "retain_outputs: unknown output name(s): "
                + ", ".join(invalid_outputs)
            )
        selected_outputs = {
            "uu" if name in {"rad", "radiance"} else name
            for name in retained_outputs
            if isinstance(name, str) and name in allowed_outputs
        }
        if normalized.get("out_mode") == "txt" and not (
            {"bbt", "uu"} & selected_outputs
        ):
            issues.append(
                "out_mode: retain_outputs must include bbt or radiance "
                "when writing text output"
            )

    rfm_config = normalized.get("rfm_config")
    if isinstance(rfm_config, Mapping) and rfm_config.get("output_mode") == "files":
        issues.append(
            "rfm_config.output_mode: run_srfm requires 'capture' optical-depth output"
        )

    _validate_rfm_config(rfm_config, issues)
    _validate_driver(
        normalized.get("driver_inputs"),
        issues,
        minimum_atmospheres=2 if runner == "iasi" else 1,
    )
    _validate_layers(normalized.get("scat_lyrs_inputs"), issues)
    _validate_prescribed_layers(normalized.get("prescribed_lyrs_inputs"), issues)
    _validate_grey_body_layers(normalized.get("gbc_lyrs_inputs"), issues)
    _validate_optical_layer_structure(normalized, issues)

    if runner == "oxharp":
        if normalized.get("sun") is True:
            for key in ("sza_cos", "saa"):
                if normalized.get(key) is None:
                    issues.append(f"{key}: required when sun is True")
        zen_cos = normalized.get("zen_cos")
        zen_sec = normalized.get("zen_sec")
        if isinstance(zen_cos, Real) and isinstance(zen_sec, Real):
            if zen_cos == 0 or not np.isclose(zen_sec, 1 / zen_cos):
                issues.append("zen_sec: must be the reciprocal of zen_cos")
        _validate_positive(
            normalized, ("maxcount", "epsilon", "gamma", "eps"), "", issues
        )
        max_workers = normalized.get("max_workers")
        if isinstance(max_workers, Integral) and max_workers <= 0:
            issues.append("max_workers: must be greater than zero")

    if runner == "iasi":
        pixel = normalized.get("px")
        if (
            isinstance(pixel, Integral)
            and not isinstance(pixel, (bool, np.bool_))
            and pixel < 0
        ):
            issues.append("px: must be non-negative")
        filename = normalized.get("iasi_fl")
        if isinstance(filename, str):
            first_separator = filename.find("_")
            last_separator = filename.rfind("_")
            if first_separator < 0 or last_separator <= first_separator:
                issues.append(
                    "iasi_fl: must contain a YYYYMMDD date between underscores"
                )
            else:
                date_text = filename[first_separator + 1 : last_separator]
                try:
                    dt.datetime.strptime(date_text, "%Y%m%d")
                except ValueError:
                    issues.append(
                        "iasi_fl: must contain a valid YYYYMMDD date between underscores"
                    )

    if issues:
        raise InputValidationError(issues)
    return normalized


def validate_srfm_inputs(values: Mapping[str, Any]) -> dict[str, Any]:
    """Validate inputs for the generic SRFM runner.

    This public validator preserves the established flat driver-table format
    and applies the geometry consumed by :mod:`srfm.main`.

    Args:
        values: Flat SRFM input mapping.

    Returns:
        A shallow normalized copy of the supplied mapping.

    Raises:
        InputValidationError: If any generic-runner input is invalid.
    """
    return _validate_inputs(values, SRFM_INPUT_SCHEMA, runner="srfm")


def validate_oxharp_inputs(values: Mapping[str, Any]) -> dict[str, Any]:
    """Validate inputs for the OXHARP-tailored SRFM runner.

    The contract accepts retrieval metadata from the real OXHARP driver while
    requiring the precomputed cosine and secant geometry used by the runner.

    Args:
        values: Merged OXHARP state and ancillary mapping.

    Returns:
        A shallow normalized copy of the supplied mapping.

    Raises:
        InputValidationError: If any OXHARP-runner input is invalid.
    """
    return _validate_inputs(values, OXHARP_INPUT_SCHEMA, runner="oxharp")


def validate_iasi_inputs(values: Mapping[str, Any]) -> dict[str, Any]:
    """Validate inputs for the processed-IASI SRFM runner.

    This contract requires the observation file, pixel, noise table, and ILS
    that :mod:`srfm.iasi_main` reads directly.

    Args:
        values: Flat processed-IASI input mapping.

    Returns:
        A shallow normalized copy of the supplied mapping.

    Raises:
        InputValidationError: If any IASI-runner input is invalid.
    """
    return _validate_inputs(values, IASI_INPUT_SCHEMA, runner="iasi")


def get_srfm_input_schema() -> dict[str, FieldSpec]:
    """Return a copy of the generic runner's public schema.

    A copy prevents inspection and tooling code from mutating the authoritative
    mapping used during execution.

    Returns:
        A shallow copy of the generic SRFM top-level schema.
    """
    return dict(SRFM_INPUT_SCHEMA)


def get_oxharp_input_schema() -> dict[str, FieldSpec]:
    """Return a copy of the OXHARP runner's public schema.

    A copy keeps callers from changing validation globally while allowing
    documentation and tooling to inspect the contract.

    Returns:
        A shallow copy of the OXHARP top-level schema.
    """
    return dict(OXHARP_INPUT_SCHEMA)


def get_iasi_input_schema() -> dict[str, FieldSpec]:
    """Return a copy of the processed-IASI runner's public schema.

    A copy keeps callers from changing validation globally while allowing
    documentation and tooling to inspect the contract.

    Returns:
        A shallow copy of the IASI top-level schema.
    """
    return dict(IASI_INPUT_SCHEMA)
