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


_MISSING = object()


class InputValidationError(ValueError):
    """Raised with all problems found in an SRFM input mapping."""

    def __init__(self, issues: Sequence[str]):
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
    "LUN",
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
    "rad": FieldSpec((bool,), required=True),
    "bbt": FieldSpec((bool,), required=True),
    "rad_out_fname": FieldSpec((str,), required=True, nullable=True),
    "bbt_out_fname": FieldSpec((str,), required=True, nullable=True),
    "plot_type": FieldSpec(
        (str,), required=True, choices=frozenset({"rad", "bbt"})
    ),
    "convolve_iasi": FieldSpec((bool,), required=True),
    "iasi_ils": FieldSpec(PATH_TYPES, required=True, nullable=True),
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
    # RFM structures.
    "rfm_config": FieldSpec(MAPPING_TYPES, required=True),
    "driver_inputs": FieldSpec(MAPPING_TYPES, required=True),
    "levels": FieldSpec(nullable=True),
    # DISORT configuration.
    "fisot": FieldSpec(NUMBER_TYPES, required=True),
    "albedo": FieldSpec(NUMBER_TYPES, required=True),
    "temis": FieldSpec(NUMBER_TYPES, required=True),
    "earth_radius": FieldSpec(NUMBER_TYPES),
    "nmom": FieldSpec(INTEGER_TYPES, required=True),
    "maxcmu": FieldSpec(INTEGER_TYPES, required=True),
    "maxumu": FieldSpec(INTEGER_TYPES, required=True),
    "maxphi": FieldSpec(INTEGER_TYPES, required=True),
    "maxulv": FieldSpec(INTEGER_TYPES, required=True),
    "usrang": FieldSpec((bool,), required=True),
    "usrtau": FieldSpec((bool,), required=True),
    "ibcnd": FieldSpec(INTEGER_TYPES, required=True, choices=frozenset({0, 1})),
    "onlyfl": FieldSpec((bool,), required=True),
    "prnt": FieldSpec((list,), required=True),
    "planck": FieldSpec((bool,), required=True),
    "lamber": FieldSpec((bool,), required=True),
    "deltamplus": FieldSpec((bool,), required=True),
    "do_pseudo_sphere": FieldSpec((bool,), required=True),
    "utau": FieldSpec(required=True),
    "disort_precision": FieldSpec(
        (str,), required=True, choices=frozenset({"single", "double"})
    ),
    "header": FieldSpec((str,)),
    "adjust_maxcmu": FieldSpec((bool,), required=True),
    "btemp": FieldSpec(NUMBER_TYPES),
    "ttemp": FieldSpec(NUMBER_TYPES),
    # Scattering and geometry.
    "scat_lyrs_inputs": FieldSpec(MAPPING_TYPES),
    "date": FieldSpec((dt.datetime, tuple)),
    "sun": FieldSpec((bool,), required=True),
    "sza": FieldSpec(NUMBER_TYPES, required=True),
    "saa": FieldSpec(NUMBER_TYPES, required=True),
    "zen": FieldSpec(NUMBER_TYPES, required=True),
    "azi": FieldSpec(NUMBER_TYPES, required=True),
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
    "center_alt": FieldSpec(NUMBER_TYPES, required=True, nullable=True),
    "thick": FieldSpec(NUMBER_TYPES, required=True, nullable=True),
    "alt_upp": FieldSpec(NUMBER_TYPES, required=True, nullable=True),
    "alt_low": FieldSpec(NUMBER_TYPES, required=True, nullable=True),
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


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes))


def _is_array_like(value: Any) -> bool:
    return _is_sequence(value) or isinstance(value, np.ndarray)


def _check_mapping(
    value: Mapping[str, Any],
    schema: Mapping[str, FieldSpec],
    path: str,
    issues: list[str],
) -> None:
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
                issues.append(f"{field_path}: may not be None")
            continue
        if spec.expected is not None:
            valid_type = isinstance(item, spec.expected)
            if spec.expected == (bool,):
                valid_type = type(item) is bool
            elif any(
                isinstance(expected_type, type)
                and issubclass(expected_type, Real)
                for expected_type in spec.expected
            ):
                valid_type = valid_type and not isinstance(item, (bool, np.bool_))
            if not valid_type:
                expected = ", ".join(t.__name__ for t in spec.expected)
                issues.append(
                    f"{field_path}: expected {expected}, got {type(item).__name__}"
                )
                continue
            if isinstance(item, Real) and not np.isfinite(item):
                issues.append(f"{field_path}: must be finite")
                continue
        if spec.choices is not None and item not in spec.choices:
            choices = ", ".join(repr(choice) for choice in sorted(spec.choices, key=str))
            issues.append(f"{field_path}: expected one of {choices}, got {item!r}")


def _validate_positive(mapping, keys, path, issues, *, allow_zero=False):
    for key in keys:
        value = mapping.get(key)
        if isinstance(value, Real) and not isinstance(value, (bool, np.bool_)):
            invalid = value < 0 if allow_zero else value <= 0
            if invalid:
                relation = "non-negative" if allow_zero else "greater than zero"
                issues.append(f"{path}{key}: must be {relation}")


def _validate_rfm_config(config: Any, issues: list[str]) -> None:
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


def _validate_driver(driver: Any, issues: list[str]) -> None:
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


def _validate_layers(layers: Any, issues: list[str]) -> None:
    if layers is None:
        return
    if not isinstance(layers, Mapping):
        return
    for layer_name, layer in layers.items():
        path = f"scat_lyrs_inputs.{layer_name}."
        if not isinstance(layer_name, str):
            issues.append(f"scat_lyrs_inputs: layer name {layer_name!r} must be a string")
            continue
        if not isinstance(layer, Mapping):
            issues.append(f"scat_lyrs_inputs.{layer_name}: expected a mapping")
            continue
        _check_mapping(layer, LAYER_SCHEMA, path, issues)
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
        if all(layer.get(key) is None for key in ("mass_loading", "n", "s_a_den", "v_den")):
            issues.append(
                f"{path}mass_loading: one of mass_loading, n, s_a_den, or v_den is required"
            )
        centre_extent = layer.get("center_alt") is not None and layer.get("thick") is not None
        bound_extent = layer.get("alt_low") is not None and layer.get("alt_upp") is not None
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
        if centre_extent and bound_extent and all(
            isinstance(layer.get(key), Real)
            for key in ("center_alt", "thick", "alt_low", "alt_upp")
        ):
            expected_low = layer["center_alt"] - layer["thick"] / 2
            expected_upp = layer["center_alt"] + layer["thick"] / 2
            if not (
                np.isclose(layer["alt_low"], expected_low)
                and np.isclose(layer["alt_upp"], expected_upp)
            ):
                issues.append(
                    f"{path}alt_low: explicit bounds must match center_alt/thick"
                )
        if layer.get("comp") == "ri" and layer.get("refractive_index") is None:
            issues.append(
                f"{path}refractive_index: required when comp is 'ri'"
            )


def validate_srfm_inputs(values: Mapping[str, Any]) -> dict[str, Any]:
    """Validate an existing flat SRFM mapping and return a shallow normalized copy."""
    if not isinstance(values, Mapping):
        raise InputValidationError(
            [f"inputs: expected a mapping, got {type(values).__name__}"]
        )

    normalized = dict(values)
    for key, spec in SRFM_INPUT_SCHEMA.items():
        if key not in normalized and spec.default is not _MISSING:
            normalized[key] = spec.default

    issues: list[str] = []
    _check_mapping(normalized, SRFM_INPUT_SCHEMA, "", issues)

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

    prnt = normalized.get("prnt")
    if "prnt" in normalized and (
        not isinstance(prnt, list)
        or len(prnt) != 5
        or not all(type(v) is bool for v in prnt)
    ):
        issues.append("prnt: must be a list containing exactly five boolean values")
    utau = normalized.get("utau")
    if "utau" in normalized:
        if not _is_array_like(utau) or len(utau) == 0:
            issues.append("utau: must be a non-empty numeric sequence")
        elif not all(
            isinstance(v, Real)
            and not isinstance(v, (bool, np.bool_))
            and np.isfinite(v)
            and v >= 0
            for v in utau
        ):
            issues.append("utau: values must be finite non-negative numbers")

    _validate_positive(normalized, ("nmom", "maxcmu", "maxulv"), "", issues)
    _validate_positive(
        normalized, ("maxumu", "maxphi"), "", issues, allow_zero=True
    )
    if normalized.get("onlyfl") is False:
        for key in ("maxumu", "maxphi"):
            if normalized.get(key) == 0:
                issues.append(f"{key}: may be zero only when onlyfl is True")
    if isinstance(normalized.get("maxcmu"), Integral) and normalized["maxcmu"] % 2:
        issues.append("maxcmu: must be even")
    if isinstance(normalized.get("maxcmu"), Integral) and normalized["maxcmu"] < 2:
        issues.append("maxcmu: must be at least 2")
    for key, lower, upper in (
        ("albedo", 0, 1),
        ("temis", 0, 1),
        ("sza", 0, 180),
        ("zen", 0, 180),
        ("saa", 0, 360),
        ("azi", 0, 360),
    ):
        value = normalized.get(key)
        if isinstance(value, Real) and not lower <= value <= upper:
            issues.append(f"{key}: must be between {lower} and {upper}")

    if normalized.get("convolve_iasi") is True and not normalized.get("iasi_ils"):
        issues.append("iasi_ils: required when convolve_iasi is True")
    if normalized.get("out_mode") in {"txt", "netcdf"} and not (
        normalized.get("rad") or normalized.get("bbt")
    ):
        issues.append("out_mode: rad or bbt must be enabled when writing output")

    header = normalized.get("header")
    if isinstance(header, str) and len(header) >= 127:
        issues.append("header: must contain fewer than 127 characters")
    earth_radius = normalized.get("earth_radius")
    if isinstance(earth_radius, Real) and earth_radius <= 0:
        issues.append("earth_radius: must be greater than zero")

    rfm_config = normalized.get("rfm_config")
    if isinstance(rfm_config, Mapping) and rfm_config.get("output_mode") == "files":
        issues.append(
            "rfm_config.output_mode: run_srfm requires 'capture' optical-depth output"
        )

    _validate_rfm_config(rfm_config, issues)
    _validate_driver(normalized.get("driver_inputs"), issues)
    _validate_layers(normalized.get("scat_lyrs_inputs"), issues)

    if issues:
        raise InputValidationError(issues)
    return normalized


def get_srfm_input_schema() -> dict[str, FieldSpec]:
    """Return a copy of the public top-level schema for inspection/tooling."""
    return dict(SRFM_INPUT_SCHEMA)
