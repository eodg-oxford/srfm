"""Reconcile SRFM vertical geometry with RFM's single precision ATMLEV rules."""

from pathlib import Path
import re
import shlex

import numpy as np


def _driver_sections(driver_inputs, configuration):
    """Return ordered relevant driver records, including an existing driver file.

    Inputs are not modified. Generated sections follow ``_compose_driver_sections``;
    file sections retain their actual order so pre-*LEV insertions are replayed.
    """
    from .rfm_helper import _load_driver_lines, _normalize_section_data

    driver_path = (configuration or {}).get("driver_path")
    if (driver_path and Path(driver_path).expanduser().exists()
            and not (configuration or {}).get("generate_driver", False)):
        sections = {}
        current = None
        for line in _load_driver_lines(Path(driver_path).expanduser()):
            line = line.split("!", 1)[0].strip()
            if line.startswith("*"):
                current = line[:4].upper()
                sections[current] = []
            elif line and current is not None:
                sections[current].append(line)
        return sections
    return {
        section: _normalize_section_data(driver_inputs.get(key))
        for section, key in (
            ("*FLG", "flags"), ("*ATM", "atmosphere"), ("*TAN", "tangent"),
            ("*JAC", "jac"), ("*LEV", "lev"), ("*OBS", "obs"), ("*SFC", "sfc"),
        )
    }


def _record_fields(records):
    """Split RFM records into fields, retaining quoted filenames and assignments."""
    return shlex.split(" ".join(
        re.sub(r"\s*=\s*", "=", record.split("!", 1)[0]) for record in records
    ))


def _interpolate_r4(coordinates, values, targets):
    """Interpolate altitude profiles with RFM VAL1DI_RR arithmetic and extrapolation.

    All arrays use R4; coordinates may increase or decrease. No input is mutated.
    This is needed only for a pressure computational grid, whose HGT profile is
    interpolated in ln(pressure) by RFM before altitude matching.
    """
    # Preserve endpoint order even for descending pressure coordinates: reversing
    # the interpolation expression can change its R4 rounding by one ulp.
    if coordinates[0] > coordinates[-1]:
        indices = np.searchsorted(-coordinates, -targets, side="right") - 1
    else:
        indices = np.searchsorted(coordinates, targets, side="right") - 1
    indices = np.clip(indices, 0, len(coordinates) - 2)
    fraction = (targets - coordinates[indices]) / (coordinates[indices + 1] - coordinates[indices])
    return (np.float32(1) - fraction) * values[indices] + fraction * values[indices + 1]


def read_matching_grid(driver_inputs, configuration, directory):
    """Read the computational altitude grid present when RFM reaches *LEV.

    Args:
        driver_inputs (Mapping): Structured RFM sections, not modified.
        configuration (Mapping): RFM configuration, including driver_path and
            generate_driver. Existing drivers retain their actual section order.
        directory (path-like): RFM working directory for relative input filenames.

    Returns:
        numpy.ndarray: R4 atmospheric altitudes in km. The first atmospheric
        profile sets the grid; LAY subdivision and preceding height/Jacobian/
        observer/surface insertions are replayed. Later composition grids do not
        replace an established altitude grid.

    Raises:
        ValueError: If the grid is absent/invalid or a requested insertion is
            outside it.
        OSError: If a required existing input file cannot be read.
    """
    from .rfm_functions import read_atm_file

    directory = Path(directory).expanduser().resolve()
    sections = _driver_sections(driver_inputs, configuration)
    flags = {field.upper() for field in _record_fields(sections.get("*FLG", []))}
    atmospheric_levels = None
    for field in _record_fields(sections.get("*ATM", [])):
        if "=" in field:
            if field.split("=", 1)[0].upper().lstrip("*") == "HGT":
                raise ValueError("A constant HGT profile cannot define increasing SRFM altitudes.")
            continue
        filename, _, qualifier = field.partition("(")
        path = Path(filename).expanduser()
        profiles = {
            label.split()[0].upper(): np.asarray(values, dtype=np.float32)
            for label, values in read_atm_file(path if path.is_absolute() else directory / path).items()
        }
        if qualifier:
            qualifier = qualifier.rstrip(")")
            excluded = qualifier.startswith("-")
            names = {name.strip().lstrip("*").upper() for name in qualifier.lstrip("-").split(";")}
            profiles = {name: values for name, values in profiles.items()
                        if (name not in names if excluded else name in names)}
        if not profiles:
            continue
        grid_name = next(iter(profiles))
        if grid_name not in {"HGT", "PRE"}:
            raise ValueError("The first RFM atmospheric profile must define HGT or PRE.")
        grid = profiles[grid_name].copy()
        if "LAY" in flags:
            subdivided = np.empty(2 * len(grid) - 1, dtype=np.float32)
            subdivided[::2] = grid
            subdivided[1::2] = (np.sqrt(grid[:-1] * grid[1:]) if grid_name == "PRE"
                                 else np.float32(0.5) * (grid[:-1] + grid[1:]))
            grid = subdivided
        if grid_name == "HGT":
            atmospheric_levels = grid
        elif "HGT" in profiles:
            atmospheric_levels = _interpolate_r4(
                np.log(profiles["PRE"]), profiles["HGT"], np.log(grid)
            )
        else:
            raise ValueError("A pressure computational grid requires a HGT profile in its first file.")
        break
    if (atmospheric_levels is None or len(atmospheric_levels) < 2
            or not np.all(np.isfinite(atmospheric_levels))
            or np.any(np.diff(atmospheric_levels) <= 0)):
        raise ValueError("RFM atmospheric altitudes must be finite and strictly increasing in R4.")

    def altitude_fields(records):
        """Expand one numeric altitude file, as supported by RFM NXTFFL."""
        for field in _record_fields(records):
            if field.upper() == "BOA":
                yield atmospheric_levels[0]
            elif field.upper() == "TOA":
                yield atmospheric_levels[-1]
            else:
                try:
                    yield np.float32(field.replace("D", "E").replace("d", "e"))
                except ValueError:
                    path = Path(field).expanduser()
                    records = (path if path.is_absolute() else directory / path).read_text().splitlines()
                    for value in _record_fields(records):
                        yield np.float32(value.replace("D", "E").replace("d", "e"))

    for section, records in sections.items():
        if section == "*LEV":
            break
        insertions = []
        if section == "*HGT" or (section == "*TAN" and "FLX" in flags):
            insertions = list(altitude_fields(records))
        elif section == "*JAC" and "JAC" in flags:
            for record in records:
                fields = _record_fields([record])[1:]  # first field is the species
                for field in fields:
                    values = (altitude_fields(sections.get("*TAN", [])) if field.upper() == "TAN"
                              else altitude_fields([field]))
                    insertions.extend(value for value in values
                                      if atmospheric_levels[0] <= value <= atmospheric_levels[-1])
        elif section in {"*OBS", "*SFC"} and section[1:] in flags:
            positional_fields = 0
            for field in _record_fields(records):
                parameter, separator, value = field.partition("=")
                if not separator and section == "*OBS":
                    positional_fields += 1
                    parameter, value = ("HGTOBS" if positional_fields == 1 else "PSIOBS"), field
                if parameter.upper() in {"HGTOBS", "HGTSFC"}:
                    altitude = np.float32(value)
                    if parameter.upper() != "HGTOBS" or altitude <= atmospheric_levels[-1]:
                        insertions.append(altitude)
                elif parameter.upper() in {"PREOBS", "PRESFC"}:
                    # In the bundled ATMLEV the denominator retains its sign.
                    # On decreasing ln(p), ALPHA is nonpositive, so an in-range
                    # pressure always selects an existing level, never inserts.
                    # Leave pressure/range validation to RFM itself.
                    continue
        for altitude in insertions:
            _, atmospheric_levels = match_altitude(atmospheric_levels, altitude)
    return atmospheric_levels


def match_altitude(atmospheric_levels, requested_altitude):
    """Match or insert an altitude using RFM's relative tolerance.

    Args:
        atmospheric_levels (numpy.ndarray): Increasing R4 altitude grid in km,
            containing at least two levels. Not modified.
        requested_altitude (float): Altitude in km, rounded to R4 before matching.

    Returns:
        tuple: Effective R4 altitude and the working grid including any insertion.
            An exact match or a fractional distance <= R4(0.001) snaps to the
            nearest level, choosing the first index on a tie. The fraction uses
            adjacent spacing at endpoints and half the neighbour span inside.

    Raises:
        ValueError: If the requested altitude is nonfinite or outside the grid.
    """
    requested_altitude = np.float32(requested_altitude)
    if (
        not np.isfinite(requested_altitude)
        or requested_altitude < atmospheric_levels[0]
        or requested_altitude > atmospheric_levels[-1]
    ):
        raise ValueError(f"Altitude {requested_altitude} km is outside the atmospheric grid.")
    distances = np.abs(atmospheric_levels - requested_altitude)
    nearest_level_index = int(np.argmin(distances))
    distance = distances[nearest_level_index]
    if distance == 0:
        return atmospheric_levels[nearest_level_index], atmospheric_levels
    if nearest_level_index == 0:
        fraction = distance / (atmospheric_levels[1] - atmospheric_levels[0])
    elif nearest_level_index == len(atmospheric_levels) - 1:
        fraction = distance / (atmospheric_levels[-1] - atmospheric_levels[-2])
    else:
        fraction = (np.float32(2) * distance) / (
            atmospheric_levels[nearest_level_index + 1]
            - atmospheric_levels[nearest_level_index - 1]
        )
    if fraction <= np.float32(0.001):
        return atmospheric_levels[nearest_level_index], atmospheric_levels
    insertion_index = np.searchsorted(atmospheric_levels, requested_altitude)
    return requested_altitude, np.insert(atmospheric_levels, insertion_index, requested_altitude)


def reconcile_levels(atmospheric_levels, requested_levels, layer_boundaries):
    """Return canonical output levels and layer assignments without mutating inputs.

    Args:
        atmospheric_levels (array-like): Actual RFM computational altitudes in km,
            including preprocessing and levels inserted before *LEV.
        requested_levels (array-like): Requested output altitudes in km. As in
            the existing layer insertion helper, points within optical layers
            are replaced by the layer's boundaries.
        layer_boundaries (Mapping): Layer name to (lower, upper) altitude in km.
            Layers must be disjoint, as required by the public runner schemas.

    Returns:
        tuple: Increasing Python-float levels representing exact R4 altitudes,
            top-down interval assignments (names or None), and a mapping of
            effective layer boundaries in km. All tolerance arithmetic is R4;
            geometry calculations downstream retain their usual precision.

    Raises:
        ValueError: For invalid grids, out-of-range levels, collapsed layers,
            overlap or shared boundaries, or an unstable canonical sequence.
            Layer geometry errors include original and effective boundaries.
    """
    atmospheric_levels = np.asarray(atmospheric_levels, dtype=np.float32).copy()
    if (
        atmospheric_levels.ndim != 1 or atmospheric_levels.size < 2
        or not np.all(np.isfinite(atmospheric_levels))
        or np.any(np.diff(atmospheric_levels) <= 0)
    ):
        raise ValueError("RFM atmospheric altitudes must be finite and strictly increasing in R4.")
    candidates = [float(value) for value in requested_levels]
    for lower, upper in layer_boundaries.values():
        candidates = [value for value in candidates if not lower <= value <= upper]
    candidates.extend(value for bounds in layer_boundaries.values() for value in bounds)
    if not candidates or not np.all(np.isfinite(candidates)):
        raise ValueError("Output altitudes must be nonempty and finite.")
    working_levels = atmospheric_levels.copy()
    matched = {}
    for requested_altitude in sorted(set(candidates)):
        effective_altitude, working_levels = match_altitude(working_levels, requested_altitude)
        matched[requested_altitude] = float(effective_altitude)
    effective_layer_boundaries = {
        name: tuple(matched[float(value)] for value in bounds)
        for name, bounds in layer_boundaries.items()
    }

    def describe_layer(name):
        """Include both geometries in an actionable layer error."""
        return (f"layer {name!r}: original {tuple(layer_boundaries[name])} km, "
                f"effective {effective_layer_boundaries[name]} km")

    ordered_layers = sorted(effective_layer_boundaries, key=effective_layer_boundaries.get)
    for index, name in enumerate(ordered_layers):
        lower, upper = effective_layer_boundaries[name]
        if lower >= upper:
            raise ValueError(f"RFM altitude matching collapsed {describe_layer(name)}.")
        if index:
            previous = ordered_layers[index - 1]
            if lower <= effective_layer_boundaries[previous][1]:
                raise ValueError(
                    "RFM altitude matching created incompatible overlap or a shared boundary: "
                    f"{describe_layer(previous)}; {describe_layer(name)}."
                )
    levels = sorted(set(matched.values()))
    if len(levels) < 2:
        raise ValueError("RFM altitude matching leaves fewer than two distinct output levels.")
    # Replaying the canonical sequence from the original grid checks that removing
    # duplicate requests has not changed later matching decisions. This is bounded
    # to one pass, with no convergence loop or additional native model call.
    working_levels = atmospheric_levels.copy()
    for altitude in levels:
        effective_altitude, working_levels = match_altitude(working_levels, altitude)
        if effective_altitude != altitude:
            raise ValueError(f"Canonical RFM altitude {altitude} km would snap again.")
    tracking = [None] * (len(levels) - 1)
    for name, (lower, upper) in effective_layer_boundaries.items():
        lower_index = levels.index(lower)
        upper_index = levels.index(upper)
        if upper_index != lower_index + 1:
            raise ValueError(f"Incompatible interval assignment for {describe_layer(name)}.")
        tracking[lower_index] = name
    return levels, tracking[::-1], effective_layer_boundaries
