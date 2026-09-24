"""Load, validate, convert, and interpolate SRFM spectral input fields.

Spectral coordinates are converted to wavenumber in inverse centimetres before
interpolation.  Solar spectral densities are likewise converted to
``W m-2 (cm-1)-1`` at their source points, including the appropriate Jacobian,
so interpolation is always performed in SRFM's computational coordinate.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


GRID_UNITS = frozenset({"cm-1", "um", "nm"})
SOLAR_VALUE_UNITS = frozenset(
    {
        "W m-2 (cm-1)-1",
        "W m-2 um-1",
        "W m-2 nm-1",
        "mW cm-2 um-1",
    }
)


def _as_column_indices(value: Any, field_path: str) -> tuple[int, ...]:
    """Return validated zero-based file-column indices.

    Args:
        value: Candidate sequence of column indices.
        field_path: Dotted input path used in error messages.

    Returns:
        Tuple of non-negative integer column indices.

    Raises:
        ValueError: If the value is not a non-empty sequence of distinct,
            non-negative integers.
    """
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise ValueError(f"{field_path}.value_columns must be a non-empty sequence")
    if any(type(column) is not int or column < 0 for column in value):
        raise ValueError(
            f"{field_path}.value_columns must contain non-negative integers"
        )
    columns = tuple(value)
    if len(set(columns)) != len(columns):
        raise ValueError(f"{field_path}.value_columns must not contain duplicates")
    return columns


def _load_source_arrays(
    specification: Mapping[str, Any],
    field_path: str,
    *,
    matrix_values: bool,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Load a spectral coordinate and value array from memory or a text file."""
    has_file = "file" in specification
    has_memory = "grid" in specification or "values" in specification
    if has_file == has_memory:
        raise ValueError(
            f"{field_path} must provide exactly one of file or grid/values"
        )

    if has_file:
        filename = Path(specification["file"])
        grid_column = specification.get("grid_column", 0)
        skiprows = specification.get("skiprows", 0)
        if type(grid_column) is not int or grid_column < 0:
            raise ValueError(f"{field_path}.grid_column must be a non-negative integer")
        if type(skiprows) is not int or skiprows < 0:
            raise ValueError(f"{field_path}.skiprows must be a non-negative integer")
        if matrix_values:
            value_columns = _as_column_indices(
                specification.get("value_columns"), field_path
            )
        else:
            value_column = specification.get("value_column", 1)
            if type(value_column) is not int or value_column < 0:
                raise ValueError(
                    f"{field_path}.value_column must be a non-negative integer"
                )
            value_columns = (value_column,)
        selected_columns = (grid_column, *value_columns)
        if len(set(selected_columns)) != len(selected_columns):
            raise ValueError(
                f"{field_path}: grid and value columns must be distinct"
            )
        try:
            loaded = np.loadtxt(
                filename,
                skiprows=skiprows,
                usecols=selected_columns,
                ndmin=2,
            )
        except (OSError, ValueError) as exc:
            raise ValueError(f"{field_path}.file: unable to load {filename}: {exc}") from exc
        grid = loaded[:, 0]
        values = loaded[:, 1:] if matrix_values else loaded[:, 1]
        provenance = {
            "source": "file",
            "file": str(filename),
            "grid_column": grid_column,
            ("value_columns" if matrix_values else "value_column"): (
                list(value_columns) if matrix_values else value_columns[0]
            ),
            "skiprows": skiprows,
        }
    else:
        if "grid" not in specification or "values" not in specification:
            raise ValueError(f"{field_path} requires both grid and values")
        try:
            grid = np.asarray(specification["grid"], dtype=float)
            values = np.asarray(specification["values"], dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{field_path}: grid and values must be numeric") from exc
        provenance = {
            "source": "memory",
            "grid": {"shape": list(grid.shape), "stored_as": "NetCDF variable"},
            "values": {
                "shape": list(values.shape),
                "stored_as": "NetCDF variable",
            },
        }
    return np.asarray(grid, dtype=float), np.asarray(values, dtype=float), provenance


def wavenumber_from_grid(grid: np.ndarray, grid_units: str) -> np.ndarray:
    """Convert a supported spectral coordinate to inverse centimetres.

    Args:
        grid: Finite, positive spectral coordinate.
        grid_units: ``"cm-1"``, ``"um"``, or ``"nm"``.

    Returns:
        Spectral coordinate in ``cm-1``.

    Raises:
        ValueError: If units are unsupported or the coordinate is invalid.
    """
    coordinate = np.asarray(grid, dtype=float)
    if coordinate.ndim != 1 or coordinate.size == 0:
        raise ValueError("spectral grid must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(coordinate)) or np.any(coordinate <= 0):
        raise ValueError("spectral grid values must be finite and greater than zero")
    if grid_units == "cm-1":
        wavenumber = coordinate.copy()
    elif grid_units == "um":
        wavenumber = 1.0e4 / coordinate
    elif grid_units == "nm":
        wavenumber = 1.0e7 / coordinate
    else:
        raise ValueError(
            f"grid_units must be one of {', '.join(sorted(GRID_UNITS))}; got {grid_units!r}"
        )
    if not np.all(np.isfinite(wavenumber)) or np.any(wavenumber <= 0):
        raise ValueError("converted wavenumber grid must be finite and greater than zero")
    return wavenumber


def convert_solar_spectral_density(
    values: np.ndarray,
    value_units: str,
    wavenumber_cm_inverse: np.ndarray,
) -> np.ndarray:
    """Convert solar spectral density to ``W m-2 (cm-1)-1``.

    Args:
        values: Solar spectral irradiance at the source grid points.
        value_units: One of the exact spellings in :data:`SOLAR_VALUE_UNITS`.
        wavenumber_cm_inverse: Source coordinate in inverse centimetres.

    Returns:
        Beam-normal irradiance in ``W m-2 (cm-1)-1``.

    Raises:
        ValueError: If ``value_units`` is unsupported.
    """
    spectral_irradiance = np.asarray(values, dtype=float)
    wavelength_um = 1.0e4 / np.asarray(wavenumber_cm_inverse, dtype=float)
    if value_units == "W m-2 (cm-1)-1":
        return spectral_irradiance.copy()
    if value_units == "W m-2 um-1":
        return spectral_irradiance * wavelength_um**2 / 1.0e4
    if value_units == "W m-2 nm-1":
        return spectral_irradiance * wavelength_um**2 / 10.0
    if value_units == "mW cm-2 um-1":
        return spectral_irradiance * wavelength_um**2 / 1.0e3
    raise ValueError(
        "value_units must be one of "
        + ", ".join(repr(unit) for unit in sorted(SOLAR_VALUE_UNITS))
        + f"; got {value_units!r}"
    )


@dataclass(frozen=True)
class SpectralField:
    """Validated spectral data stored on an ascending wavenumber grid."""

    wavenumber_cm_inverse: np.ndarray
    values: np.ndarray
    field_path: str
    grid_units: str
    value_units: str | None
    source_metadata: Mapping[str, Any]

    @classmethod
    def from_specification(
        cls,
        specification: Mapping[str, Any],
        field_path: str,
        *,
        matrix_values: bool = False,
        minimum: float | None = None,
        maximum: float | None = None,
        value_units: str | None = None,
        solar_density: bool = False,
    ) -> "SpectralField":
        """Build a field from its in-memory or file-backed representation.

        Args:
            specification: Mapping containing a grid/value pair or text-file columns.
            field_path: Dotted configuration path used in diagnostics.
            matrix_values: Require values shaped ``(spectral, component)``.
            minimum: Optional inclusive lower bound for every value.
            maximum: Optional inclusive upper bound for every value.
            value_units: Physical units recorded for provenance.
            solar_density: Convert values to ``W m-2 (cm-1)-1``.

        Returns:
            Validated field ordered by increasing wavenumber.

        Raises:
            TypeError: If ``specification`` is not a mapping.
            ValueError: If data, units, dimensions, or bounds are invalid.
        """
        if not isinstance(specification, Mapping):
            raise TypeError(f"{field_path} must be a mapping")
        grid_units = specification.get("grid_units")
        if grid_units not in GRID_UNITS:
            raise ValueError(
                f"{field_path}.grid_units must be one of "
                + ", ".join(repr(unit) for unit in sorted(GRID_UNITS))
            )
        grid, values, source_metadata = _load_source_arrays(
            specification, field_path, matrix_values=matrix_values
        )
        if grid.ndim != 1 or grid.size < 2:
            raise ValueError(
                f"{field_path}.grid must be a one-dimensional array with at least two points"
            )
        if matrix_values:
            if values.ndim != 2 or values.shape[0] != grid.size or values.shape[1] == 0:
                raise ValueError(
                    f"{field_path}.values must have shape (spectral_points, moments)"
                )
        elif values.ndim != 1 or values.shape != grid.shape:
            raise ValueError(f"{field_path}.values must match the one-dimensional grid")
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{field_path}.values must contain only finite values")
        if minimum is not None and np.any(values < minimum):
            raise ValueError(f"{field_path}.values must be greater than or equal to {minimum:g}")
        if maximum is not None and np.any(values > maximum):
            raise ValueError(f"{field_path}.values must be less than or equal to {maximum:g}")

        wavenumber = wavenumber_from_grid(grid, grid_units)
        differences = np.diff(wavenumber)
        if not (np.all(differences > 0) or np.all(differences < 0)):
            raise ValueError(f"{field_path}.grid must be strictly monotonic")
        if solar_density:
            if value_units is None:
                raise ValueError(f"{field_path}.value_units is required")
            values = convert_solar_spectral_density(values, value_units, wavenumber)
            if not np.all(np.isfinite(values)) or np.any(values < 0):
                raise ValueError(
                    f"{field_path}.values are invalid after spectral-density conversion"
                )
        if wavenumber[0] > wavenumber[-1]:
            wavenumber = wavenumber[::-1].copy()
            values = values[::-1].copy()
        else:
            wavenumber = wavenumber.copy()
            values = values.copy()
        metadata = dict(source_metadata)
        metadata["grid_units"] = grid_units
        if value_units is not None:
            metadata["value_units"] = value_units
        return cls(
            wavenumber_cm_inverse=wavenumber,
            values=values,
            field_path=field_path,
            grid_units=grid_units,
            value_units=("W m-2 (cm-1)-1" if solar_density else value_units),
            source_metadata=metadata,
        )

    @property
    def component_count(self) -> int:
        """Return one for scalar fields or the matrix's second dimension."""
        return 1 if self.values.ndim == 1 else int(self.values.shape[1])

    def validate_coverage(self, target_wavenumber_cm_inverse: Any) -> None:
        """Require the field to bracket an entire target spectral grid.

        Args:
            target_wavenumber_cm_inverse: Computational grid in ``cm-1``.

        Raises:
            ValueError: If the target is invalid or requires extrapolation.
        """
        target = np.asarray(target_wavenumber_cm_inverse, dtype=float)
        if target.ndim != 1 or target.size == 0 or not np.all(np.isfinite(target)):
            raise ValueError("computational wavenumber grid must be a finite 1D array")
        target_low = float(np.min(target))
        target_high = float(np.max(target))
        scale = max(abs(target_low), abs(target_high), 1.0)
        tolerance = 1.0e-12 * scale
        source_low = float(self.wavenumber_cm_inverse[0])
        source_high = float(self.wavenumber_cm_inverse[-1])
        if source_low > target_low + tolerance or source_high < target_high - tolerance:
            raise ValueError(
                f"{self.field_path} covers {source_low:g} to {source_high:g} cm-1, "
                f"but the computational grid spans {target_low:g} to {target_high:g} cm-1"
            )

    def interpolate(self, target_wavenumber_cm_inverse: Any) -> np.ndarray:
        """Linearly interpolate the field in wavenumber without extrapolation."""
        target = np.asarray(target_wavenumber_cm_inverse, dtype=float)
        self.validate_coverage(target)
        if self.values.ndim == 1:
            return np.interp(target, self.wavenumber_cm_inverse, self.values)
        output = np.empty((target.size, self.values.shape[1]), dtype=float)
        for component_index in range(self.values.shape[1]):
            output[:, component_index] = np.interp(
                target,
                self.wavenumber_cm_inverse,
                self.values[:, component_index],
            )
        return output

    def provenance(self) -> dict[str, Any]:
        """Return compact JSON-safe source metadata without numerical arrays."""
        metadata = dict(self.source_metadata)
        metadata["canonical_grid_units"] = "cm-1"
        metadata["source_point_count"] = int(self.wavenumber_cm_inverse.size)
        metadata["component_count"] = self.component_count
        if self.value_units is not None:
            metadata["canonical_value_units"] = self.value_units
        return metadata
