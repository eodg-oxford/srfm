"""Helpers for running the RFM Fortran model from Python.

- Name: rfm_helper
- Parent package: srfm
- Author: Antonin Knizek
- Contributors:
- Date: 14 Nov 2025
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, List, Mapping, Sequence, Literal
import os
import pickle
import shutil
import subprocess
import sys
import tempfile
import textwrap
import traceback
import warnings

import numpy as np
import pandas as pd

from .RFM import rfm_py


@dataclass
class RunResult:
    """Encapsulates the outcome of a single RFM execution.

    Attributes:
        status (int): Exit status reported by the underlying executable or
            shared library (zero indicates success).
        removed_files (list[pathlib.Path]): Output artefacts deleted prior to
            the run when ``clean_before`` is enabled.
        output (Any | None): Optional payload gathered from the run, such as
            captured optical-depth data or a mapping of output filenames to
            their in-memory contents.
        output_df (Any | None): Optional dataframe payload reserved for
            optical-depth capture while still allowing ``output`` to represent
            file-based artefacts.
    """

    status: int
    removed_files: List[Path]
    output: Any | None = None
    output_df: Any | None = None

    @property
    def ok(self) -> bool:
        """True when the RFM run completed successfully."""
        return self.status == 0


def clean_outputs(
    directory: Path | str | None = None,
    patterns: Sequence[str] | None = None,
) -> List[Path]:
    """Delete known RFM artefacts prior to a new run.

    Args:
        directory (Path | str | None): Directory containing the outputs.
            Defaults to the current working directory.
        patterns (Sequence[str] | None): Filename glob patterns to remove.
            Defaults to ``("*.asc", "rfm.log*")``.

    Returns:
        list[pathlib.Path]: Paths that were removed.
    """

    root = Path(directory) if directory is not None else Path.cwd()
    globs: Iterable[str] = patterns or ("*.asc", "rfm.log*")

    removed: List[Path] = []
    for pattern in globs:
        for candidate in root.glob(pattern):
            if candidate.is_file():
                candidate.unlink()
                removed.append(candidate)

    return removed


def _run_rfm_impl(
    run_id: str | None = None,
    *,
    clean_before: bool = False,
    directory: Path | str | None = None,
    patterns: Sequence[str] | None = None,
    driver_lines: Sequence[str] | None = None,
    driver_path: Path | str | None = None,
    output_mode: Literal["files", "capture"] = "files",
    enable_capture: bool | None = None,
    optical_levels: Sequence[float] | None = None,
    optical_spectrum_index: int = 1,
    optical_match_tol: float = 1e-6,
) -> RunResult:
    """Execute the compiled RFM model via ``rfm_py.rfm_run``.

    Args:
        run_id (str | None): Identifier appended to output filenames.
        clean_before (bool): When ``True`` remove artefacts matching
            ``patterns`` before running.
        directory (Path | str | None): Working directory containing the driver
            and output files. Defaults to the current working directory or the
            parent directory of ``driver_path`` when supplied.
        patterns (Sequence[str] | None): Filename glob patterns passed to
            :func:`clean_outputs` when ``clean_before`` is set.
        driver_lines (Sequence[str] | None): Contents of ``rfm.drv`` provided
            in memory; when omitted the file on disk is used.
        driver_path (Path | str | None): Path to a driver table file to load.
            Mutually exclusive with ``driver_lines``.
        output_mode (Literal["files", "capture"]): When ``"capture"`` return
            in-memory payloads without leaving artefacts on disk; ``"files"``
            executes the legacy workflow that writes outputs alongside
            ``rfm.drv`` (default ``"files"``).
        enable_capture (bool | None): Explicitly enable or disable optical-depth
            capture in the Fortran backend when running via the shared library.
            Defaults to mirroring ``output_mode == "capture"``.
        optical_levels (Sequence[float] | None): Output levels (km) required
            when ``output_mode`` is ``"capture"``.
        optical_spectrum_index (int): Spectral index (1-based) to use when
            collecting captured optical depths (default ``1``).
        optical_match_tol (float): Absolute tolerance applied when matching
            requested levels to the captured profile grid (default ``1e-6`` km).

    Returns:
        RunResult: Status code, removed files, and optional in-memory payload.

    Raises:
        ValueError: Raised when mutually exclusive parameters are supplied or
            when capture output is requested without the required metadata.
        FileNotFoundError: Raised if the working directory or driver cannot be
            located.
        RuntimeError: Raised when the model reports a non-zero status code.
    """

    mode = output_mode.lower()
    if mode not in {"files", "capture"}:
        raise ValueError(
            f"Unsupported output_mode '{output_mode}'. Expected 'files' or 'capture'."
        )
    capture_requested = mode == "capture"
    if driver_lines is not None and driver_path is not None:
        raise ValueError("driver_lines and driver_path cannot be used together.")

    driver_path_obj: Path | None = None
    if driver_path is not None:
        driver_path_obj = Path(driver_path).expanduser()
        if not driver_path_obj.exists():
            raise FileNotFoundError(f"Driver table not found: {driver_path_obj}")
        if driver_path_obj.is_dir():
            raise IsADirectoryError(
                f"Driver table path must reference a file: {driver_path_obj}"
            )
        driver_path_obj = driver_path_obj.resolve()

    if directory is not None:
        root = Path(directory).expanduser()
    elif driver_path_obj is not None:
        root = driver_path_obj.parent
    else:
        root = Path.cwd()
    root = root.resolve()

    if not root.exists():
        raise FileNotFoundError(f"Working directory does not exist: {root}")
    if not root.is_dir():
        raise NotADirectoryError(f"Working directory is not a directory: {root}")

    buffer_driver_lines: Sequence[str] | None = None
    driver_copy_path: Path | None = None
    original_driver_bytes: bytes | None = None
    driver_copy_created = False

    if driver_lines is not None:
        buffer_driver_lines = list(driver_lines)
        driver_copy_path = root / "rfm.drv"
        if driver_copy_path.exists():
            original_driver_bytes = driver_copy_path.read_bytes()
        driver_copy_path.write_text(
            "\n".join(buffer_driver_lines) + "\n", encoding="utf-8"
        )
        driver_copy_created = True
    elif driver_path_obj is not None:
        buffer_driver_lines = _load_driver_lines(driver_path_obj)
        driver_copy_path = root / "rfm.drv"
        if driver_copy_path.resolve() != driver_path_obj:
            if driver_copy_path.exists():
                original_driver_bytes = driver_copy_path.read_bytes()
            shutil.copy2(driver_path_obj, driver_copy_path)
            driver_copy_created = True
        else:
            driver_copy_created = False
    else:
        driver_copy_path = root / "rfm.drv"
        if not driver_copy_path.exists():
            raise FileNotFoundError(f"Missing RFM driver file: {driver_copy_path}")
        buffer_driver_lines = _load_driver_lines(driver_copy_path)

    removed: List[Path] = []
    if clean_before:
        removed = clean_outputs(root, patterns)

    def _invoke() -> int:
        # Forward whichever optional pieces we have to the shared library.
        kwargs: dict[str, Any] = {}
        kwargs["driver_lines"] = []
        should_enable_capture = enable_capture
        if should_enable_capture is None:
            should_enable_capture = capture_requested
        if should_enable_capture:
            kwargs["enable_capture"] = True
        if run_id is not None:
            kwargs["run_id_in"] = run_id
        return rfm_py.rfm_run(**kwargs)

    status = 0
    try:
        if mode == "files":
            exe_candidates: list[Path] = [
                root / "rfm_original",
                root.parent / "rfm_original",
            ]
            executable: Path | None = next(
                (candidate for candidate in exe_candidates if candidate.exists()), None
            )
            if executable is None:
                raise FileNotFoundError(
                    f"rfm_original executable not found near working directory {root}"
                )
            with _temporary_cwd(root):
                completed = subprocess.run(
                    [str(executable.resolve())],
                    text=True,
                    capture_output=True,
                )
            status = completed.returncode
            if status != 0:
                message = textwrap.dedent(
                    f"""
                    RFM executable failed with status {status}
                    stdout:
                    {completed.stdout.strip()}
                    stderr:
                    {completed.stderr.strip()}
                    """
                ).strip()
                raise RuntimeError(message)
        else:
            with _temporary_cwd(root):
                status = _invoke()
    finally:
        if driver_copy_created and driver_copy_path is not None:
            if original_driver_bytes is None:
                try:
                    driver_copy_path.unlink()
                except OSError:
                    pass
            else:
                driver_copy_path.write_bytes(original_driver_bytes)
    if status != 0:
        error_message = None
        fetch_error = getattr(rfm_py, "_rfm_py_error", None)
        if callable(fetch_error):
            try:
                message = str(fetch_error()).strip()
            except Exception:
                message = ""
            if message:
                error_message = message
        if error_message:
            raise RuntimeError(f"RFM run failed with status {status}: {error_message}")
        raise RuntimeError(f"RFM run failed with status {status}")

    payload = None
    payload_df = None
    if capture_requested:
        if optical_levels is not None:
            payload_df = get_captured_optical_depths(
                optical_levels,
                spectrum_index=optical_spectrum_index,
                match_tol=optical_match_tol,
            )
            payload = payload_df
        else:
            payload = {}

    return RunResult(
        status=status, removed_files=removed, output=payload, output_df=payload_df
    )


def run_rfm(
    run_id: str | None = None,
    *,
    clean_before: bool = False,
    directory: Path | str | None = None,
    patterns: Sequence[str] | None = None,
    driver_lines: Sequence[str] | None = None,
    driver_path: Path | str | None = None,
    output_mode: Literal["files", "capture"] = "files",
    enable_capture: bool | None = None,
    optical_levels: Sequence[float] | None = None,
    optical_spectrum_index: int = 1,
    optical_match_tol: float = 1e-6,
) -> RunResult:
    """Public wrapper that runs the Fortran backend in a fresh subprocess."""

    if os.environ.get("RFM_DISABLE_SUBPROCESS") == "1":
        return _run_rfm_impl(
            run_id,
            clean_before=clean_before,
            directory=directory,
            patterns=patterns,
            driver_lines=driver_lines,
            driver_path=driver_path,
            output_mode=output_mode,
            enable_capture=enable_capture,
            optical_levels=optical_levels,
            optical_spectrum_index=optical_spectrum_index,
            optical_match_tol=optical_match_tol,
        )

    payload = {
        "args": (run_id,),
        "kwargs": dict(
            clean_before=clean_before,
            directory=directory,
            patterns=patterns,
            driver_lines=driver_lines,
            driver_path=driver_path,
            output_mode=output_mode,
            enable_capture=enable_capture,
            optical_levels=optical_levels,
            optical_spectrum_index=optical_spectrum_index,
            optical_match_tol=optical_match_tol,
        ),
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        args_path = Path(tmpdir) / "rfm_worker_args.pkl"
        result_path = Path(tmpdir) / "rfm_worker_result.pkl"
        with args_path.open("wb") as handle:
            pickle.dump(payload, handle)

        env = os.environ.copy()
        env["RFM_DISABLE_SUBPROCESS"] = "1"
        command = [
            sys.executable,
            "-c",
            (
                "from pathlib import Path; import sys; from srfm import rfm_helper; "
                "rfm_helper._run_rfm_worker_entry(Path(sys.argv[1]), Path(sys.argv[2]))"
            ),
            str(args_path),
            str(result_path),
        ]
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            cwd=str(Path.cwd()),
            env=env,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                "RFM worker subprocess failed to start:\n"
                f"stdout:\n{completed.stdout}\n\nstderr:\n{completed.stderr}"
            )
        if not result_path.exists():
            raise RuntimeError(
                "RFM worker subprocess did not produce a result payload."
            )

        with result_path.open("rb") as handle:
            status, data = pickle.load(handle)

    if status == "ok":
        return data
    raise RuntimeError(f"RFM subprocess failed:\n{data}")


def _load_driver_lines(driver_path: Path) -> list[str]:
    """Read a driver table from disk while preserving trailing whitespace.

    Args:
        driver_path (Path): Location of the ``rfm.drv`` file to load.

    Returns:
        list[str]: Driver lines excluding newline delimiters.
    """

    with driver_path.open("r", encoding="utf-8") as handle:
        return [line.rstrip("\r\n") for line in handle]


def get_captured_optical_depths(
    levels: Sequence[float],
    *,
    spectrum_index: int = 1,
    match_tol: float = 1e-6,
) -> pd.DataFrame:
    """Return the optical-depth dataframe captured by the Fortran backend.

    Args:
        levels (Sequence[float]): Altitude levels (km) defining layer bounds.
        spectrum_index (int): Spectral range index (1-based) to retrieve.
        match_tol (float): Absolute tolerance applied when matching requested
            levels to the captured grid.

    Returns:
        pandas.DataFrame: Layer metadata plus differential/integrated optical
        depths mirroring the legacy ``*.opt`` parsing.

    Raises:
        ValueError: Raised when the provided levels are invalid or cannot be
            matched to the captured profile.
        RuntimeError: Raised if the Fortran capture reports a failure.
    """

    if levels is None:
        raise ValueError("levels must contain at least two entries.")

    size_result = rfm_py.rfm_get_optical_grid_size(ispc=spectrum_index)
    if isinstance(size_result, tuple):
        status, n_levels, n_points = size_result
    else:
        status = size_result
        n_levels = n_points = 0
    if status != 0:
        raise RuntimeError(f"RFM optical-depth capture failed with status {status}.")

    status, altitudes, pressures, temperatures, wavenumbers, cumulative = (
        rfm_py.rfm_get_optical_grid(
            int(n_levels),
            int(n_points),
            ispc=spectrum_index,
        )
    )
    if status != 0:
        raise RuntimeError(f"RFM optical-depth capture failed with status {status}.")

    altitudes = np.array(altitudes, dtype=float, copy=True)
    pressures = np.array(pressures, dtype=float, copy=True)
    temperatures = np.array(temperatures, dtype=float, copy=True)
    wavenumbers = np.array(wavenumbers, dtype=float, copy=True)
    cumulative = np.array(cumulative, dtype=float, copy=True, order="C")

    if altitudes.ndim != 1 or cumulative.ndim != 2:
        raise ValueError("Unexpected dimensionality in captured optical-depth data.")

    requested_levels = np.atleast_1d(np.asarray(levels, dtype=float))
    if requested_levels.size < 2:
        raise ValueError("At least two output levels are required to form layers.")

    matched = np.zeros(requested_levels.shape, dtype=bool)
    selected_indices: list[int] = []
    for idx, value in enumerate(altitudes):
        matches = np.where(
            ~matched & np.isclose(value, requested_levels, atol=match_tol, rtol=0.0)
        )[0]
        if matches.size > 0:
            selected_indices.append(idx)
            matched[matches[0]] = True

    if not np.all(matched):
        missing = ", ".join(f"{val:g}" for val in requested_levels[~matched])
        raise ValueError(f"Requested levels not present in captured profile: {missing}")

    if len(selected_indices) < 2:
        raise ValueError(
            "Insufficient matching levels to compute layer optical depths."
        )

    selected_alt = altitudes[selected_indices]
    selected_pre = pressures[selected_indices]
    selected_tem = temperatures[selected_indices]
    selected_cumulative = cumulative[selected_indices, :]

    layer_delta = selected_cumulative[:-1, :] - selected_cumulative[1:, :]
    layer_delta = layer_delta[::-1, :]
    integrated = np.cumsum(layer_delta, axis=0)

    p_upper = selected_pre[1:][::-1]
    p_lower = selected_pre[:-1][::-1]
    p_avg = ((selected_pre[1:] - selected_pre[:-1]) / 2.0 + selected_pre[:-1])[::-1]

    h_upper = selected_alt[1:][::-1]
    h_lower = selected_alt[:-1][::-1]
    h_avg = ((selected_alt[1:] - selected_alt[:-1]) / 2.0 + selected_alt[:-1])[::-1]

    t_upper = selected_tem[1:][::-1]
    t_lower = selected_tem[:-1][::-1]
    t_avg = ((selected_tem[1:] - selected_tem[:-1]) / 2.0 + selected_tem[:-1])[::-1]

    for idx in range(t_avg.size - 1):
        if abs(t_avg[idx + 1] - t_avg[idx]) >= 10.0:
            warnings.warn(
                f"Temperature step between {h_avg[idx + 1]} and {h_avg[idx]} km "
                "larger than 10K, calculation may be inaccurate in DISORT.",
                stacklevel=2,
            )

    if layer_delta.shape[1] != wavenumbers.size:
        raise ValueError("Optical-depth grid and wavenumber axis are inconsistent.")

    dod_col_names = [f"dOD_{val:.4f}" for val in wavenumbers]
    iod_col_names = [f"iOD_{val:.4f}" for val in wavenumbers]

    layer_count = layer_delta.shape[0]
    prf_df = pd.DataFrame()
    prf_df["layer no."] = range(layer_count)
    prf_df.loc[:, "p_upper (mbar)"] = p_upper
    prf_df.loc[:, "p_lower (mbar)"] = p_lower
    prf_df.loc[:, "p_avg (mbar)"] = p_avg
    prf_df.loc[:, "h_upper (km)"] = h_upper
    prf_df.loc[:, "h_lower (km)"] = h_lower
    prf_df.loc[:, "h_avg (km)"] = h_avg
    prf_df.loc[:, "T_upper (K)"] = t_upper
    prf_df.loc[:, "T_lower (K)"] = t_lower
    prf_df.loc[:, "T_avg (K)"] = t_avg

    dod_df = pd.DataFrame(layer_delta, columns=dod_col_names)
    iod_df = pd.DataFrame(integrated, columns=iod_col_names)

    return pd.concat([prf_df, dod_df, iod_df], axis=1, join="outer")


@dataclass(frozen=True)
class DriverSectionInfo:
    """Metadata describing a section within ``rfm.drv``."""

    key: str
    mandatory: bool
    summary: str
    value_format: str
    aliases: tuple[str, ...] = ()


@dataclass(frozen=True)
class SpectralRange:
    """In-memory representation of a ``*SPC`` spectral range entry.

    Attributes:
        start (float): Lower bound of the spectral interval in cm^-1 (or GHz
            when the ``GHZ`` flag is set).
        stop (float): Upper bound of the interval in cm^-1.
        spacing (float): Sampling spacing (positive) or number of points per
            cm^-1 (>1).
        label (str | None): Optional label written ahead of the range.
    """

    start: float
    stop: float
    spacing: float
    label: str | None = None

    def as_record(self) -> str:
        body = f"{self.start:g} {self.stop:g} {self.spacing:g}"
        return f"{self.label} {body}" if self.label else body


@dataclass(frozen=True)
class SpectralFile:
    """Reference to a pre-defined spectral range file for ``*SPC``.

    Attributes:
        path (pathlib.Path): Path to the ``.spc`` file (coerced to
            :class:`~pathlib.Path`).
        label (str | None): Optional label written ahead of the filename.
    """

    path: Path
    label: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", Path(self.path))

    def as_record(self) -> str:
        body = str(self.path)
        return f"{self.label} {body}" if self.label else body


@dataclass(frozen=True)
class SectionLine:
    """Explicit list of tokens that form a single driver record.

    Attributes:
        parts (tuple[Any, ...]): Token sequence written verbatim into the
            driver table line.
    """

    parts: tuple[Any, ...]

    def __init__(self, *parts: Any) -> None:
        object.__setattr__(self, "parts", tuple(parts))

    def as_record(self) -> str:
        return " ".join(_format_token(part) for part in self.parts)


FLAG_CODES: tuple[str, ...] = (
    "ABS",
    "AVG",
    "BBT",
    "BFX",
    "BIN",
    "C32",
    "CHI",
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

OUTPUT_FLAG_CODES: tuple[str, ...] = (
    "ABS",
    "BIN",
    "BBT",
    "COO",
    "DBL",
    "GHZ",
    "JAC",
    "JTP",
    "LEV",
    "LOS",
    "NEW",
    "OPT",
    "PRF",
    "PTH",
    "RAD",
    "RJT",
    "SHH",
    "TAB",
    "TRA",
    "WID",
)

OUTPUT_FLAG_CODES_SET: set[str] = set(OUTPUT_FLAG_CODES)


DRIVER_SECTIONS: tuple[DriverSectionInfo, ...] = (
    DriverSectionInfo(
        "*HDR",
        True,
        "Free-form header written to the RFM log.",
        "One or more text records; the first non-blank line becomes the header.",
    ),
    DriverSectionInfo(
        "*FLG",
        True,
        "List of three-character option flags.",
        "Whitespace-separated flag codes (see ``FLAG_CODES``).",
    ),
    DriverSectionInfo(
        "*SPC",
        True,
        "Spectral ranges and resolutions to model.",
        "Records may specify ``start stop spacing`` triples, labelled ranges, or ``.spc`` filenames.",
    ),
    DriverSectionInfo(
        "*GAS",
        True,
        "Absorbing species to include.",
        "Gas names, ``.gas`` file names, or ``*`` wildcard optionally with qualifiers ``(i)`` or ``(i:j)``.",
    ),
    DriverSectionInfo(
        "*ATM",
        True,
        "Atmospheric profiles and single-valued overrides.",
        "One or more ``.atm`` files and/or ``PARAM=value`` assignments (eg ``TEM``, ``PRE``, ``HGT``).",
    ),
    DriverSectionInfo(
        "*TAN",
        True,
        "Ray geometry specification.",
        "Tangent altitudes, secants or ``UNITS=value`` for homogeneous paths; accepts filenames expanded via ``NXTFFL``.",
        aliases=("*SEC", "*ELE", "*GEO", "*HGT", "*LEN"),
    ),
    DriverSectionInfo(
        "*DIM",
        False,
        "Tabulation axes for ``TAB`` mode (replaces ``*TAN``).",
        "Pressure, temperature and VMR axis definitions as counts/limits or filenames; accepts ``PCG``/``PLV`` shortcuts.",
    ),
    DriverSectionInfo(
        "*CIA",
        False,
        "Collision-induced absorption datasets.",
        "One or more ``.cia`` filenames or a single filename template containing ``*``.",
    ),
    DriverSectionInfo(
        "*FIN",
        False,
        "Fine-mesh spectral sampling when convolving with an ILS.",
        "Single numeric value (spacing or points per cm-1; values <1 treated as cm-1 spacing or GHz when ``GHZ`` flag).",
    ),
    DriverSectionInfo(
        "*FOV",
        False,
        "Field-of-view definition.",
        "Single ``.fov`` filename.",
    ),
    DriverSectionInfo(
        "*GRD",
        False,
        "Irregular spectral grids.",
        "``.grd`` filenames or a template containing ``*``.",
    ),
    DriverSectionInfo(
        "*HIT",
        False,
        "HITRAN/ACE line-transition data sources.",
        "Records of ``BINFIL|PARFIL|HDBFIL = path`` or plain filenames (auto-typed).",
    ),
    DriverSectionInfo(
        "*ILS",
        False,
        "Instrument line-shape definitions.",
        "One or more ``.ils`` filenames.",
    ),
    DriverSectionInfo(
        "*JAC",
        False,
        "Jacobian targets and altitude grids.",
        "Records of ``target [z_min TAN z_max]`` or explicit altitude lists; accepts ``TAN`` shortcut.",
    ),
    DriverSectionInfo(
        "*LEV",
        False,
        "Intermediate altitude levels for output.",
        "Ascending list of altitude values (km).",
    ),
    DriverSectionInfo(
        "*LUT",
        False,
        "Absorption look-up tables.",
        "``.lut`` filenames or a template with ``*`` (one per gas/range).",
    ),
    DriverSectionInfo(
        "*NTE",
        False,
        "Non-LTE vibrational temperature datasets.",
        "``.nte`` filenames or a template containing ``*`` (optional ``(psi)`` suffix).",
    ),
    DriverSectionInfo(
        "*OBS",
        False,
        "Observer location for limb geometries.",
        "Parameters ``HGTOBS|PREOBS`` and optional ``PSIOBS`` via ``PARAM=value`` or positional fields.",
    ),
    DriverSectionInfo(
        "*OUT",
        False,
        "Override default output filenames and directory.",
        "Records ``PARAM = VALUE`` for ``OUTDIR``, ``ABSFIL``, ``BBTFIL``, ``COOFIL``, ``OPTFIL``, ``PRFFIL``, ``PTHFIL``, ``RADFIL``, ``RJTFIL``, ``TABFIL``, ``TRAFIL``, ``WIDFIL``.",
    ),
    DriverSectionInfo(
        "*PHY",
        False,
        "Adjust physical constants.",
        "``PARAM = value`` for ``CPKMOL``, ``GRAVTY``, ``WGTAIR``, ``TEMSPA``, ``RADCRV``, ``WNORFR`` or ``GHZRFR`` (converted to ``WNORFR``).",
    ),
    DriverSectionInfo(
        "*REJ",
        False,
        "Line-strength rejection thresholds.",
        "Records ``molecule/isotope strength`` or ``* strength`` (cm-1/(molecule cm-2)).",
    ),
    DriverSectionInfo(
        "*SFC",
        False,
        "Surface boundary conditions.",
        "Parameters ``TEMSFC``/``TEMREL``, ``EMSSFC`` (file or numeric), ``HGTSFC``/``PRESFC`` and optional ``RFLSFC``.",
    ),
    DriverSectionInfo(
        "*SHP",
        False,
        "Line-shape assignments per molecule.",
        "Records starting with ``VOI|LOR|DOP|CHI|VAN|VVW`` followed by gas names or ``*`` default.",
    ),
    DriverSectionInfo(
        "*SVD",
        False,
        "SVD-compressed absorption tables.",
        "``.svd`` filenames or template with ``*``.",
    ),
    DriverSectionInfo(
        "*XSC",
        False,
        "Continuum/foreign cross-section data.",
        "``.xsc`` filenames or template with ``*``.",
    ),
)


SECTION_WRITE_ORDER: tuple[str, ...] = tuple(
    info.key for info in DRIVER_SECTIONS if info.key not in ("*DIM",)
) + ("*DIM",)


def get_rfm_input_catalog() -> dict[str, Any]:
    """Compile metadata about supported driver sections and aliases.

    Returns:
        dict[str, Any]: Mapping with keys ``sections``, ``mandatory_sections``,
            ``optional_sections``, ``flag_codes``, and ``aliases``.
    """

    mandatory = tuple(info.key for info in DRIVER_SECTIONS if info.mandatory)
    optional = tuple(info.key for info in DRIVER_SECTIONS if not info.mandatory)
    alias_map = {alias: info.key for info in DRIVER_SECTIONS for alias in info.aliases}
    return {
        "sections": DRIVER_SECTIONS,
        "mandatory_sections": mandatory,
        "optional_sections": optional,
        "flag_codes": FLAG_CODES,
        "aliases": alias_map,
    }


def run_rfm_with_parameters(
    *,
    header: str | Sequence[str],
    flags: Sequence[str],
    spectral: Any,
    gases: Any,
    atmosphere: Any,
    tangent: Any | None = None,
    tab_dimensions: Any | None = None,
    cia: Any | None = None,
    fin: Any | None = None,
    fov: Any | None = None,
    grd: Any | None = None,
    hit: Any | None = None,
    ils: Any | None = None,
    jac: Any | None = None,
    lev: Any | None = None,
    lut: Any | None = None,
    nte: Any | None = None,
    obs: Any | None = None,
    out: Any | None = None,
    phy: Any | None = None,
    rej: Any | None = None,
    sfc: Any | None = None,
    shp: Any | None = None,
    svd: Any | None = None,
    xsc: Any | None = None,
    extra_sections: Mapping[str, Any] | None = None,
    run_id: str | None = None,
    clean_before: bool = False,
    directory: Path | str | None = None,
    patterns: Sequence[str] | None = None,
    output_mode: Literal["files", "capture"] = "files",
    enable_capture: bool | None = None,
    optical_levels: Sequence[float] | None = None,
    optical_spectrum_index: int = 1,
    optical_match_tol: float = 1e-6,
) -> RunResult:
    """Run RFM by constructing inputs directly from Python structures.

    Args:
        header (str | Sequence[str]): Text for the ``*HDR`` section.
        flags (Sequence[str]): Recognised three-letter flag codes
            (case-insensitive).
        spectral (Any): Entries for the ``*SPC`` section; accepts numbers,
            strings, :class:`SpectralRange`, :class:`SpectralFile`,
            :class:`SectionLine`, mappings, or nested sequences.
        gases (Any): Content for ``*GAS`` (molecule names, filenames,
            wildcards, etc.).
        atmosphere (Any): Content for ``*ATM``; either filenames or
            ``PARAM=value`` assignments.
        tangent (Any | None): Tangent path specification for ``*TAN`` (required
            unless the ``TAB`` flag is set).
        tab_dimensions (Any | None): ``*DIM`` records used when the ``TAB``
            flag is enabled.
        cia (Any | None): Content for ``*CIA``.
        fin (Any | None): Content for ``*FIN``.
        fov (Any | None): Content for ``*FOV``.
        grd (Any | None): Content for ``*GRD``.
        hit (Any | None): Content for ``*HIT``.
        ils (Any | None): Content for ``*ILS``.
        jac (Any | None): Content for ``*JAC``.
        lev (Any | None): Content for ``*LEV``.
        lut (Any | None): Content for ``*LUT``.
        nte (Any | None): Content for ``*NTE``.
        obs (Any | None): Content for ``*OBS``.
        out (Any | None): Content for ``*OUT``.
        phy (Any | None): Content for ``*PHY``.
        rej (Any | None): Content for ``*REJ``.
        sfc (Any | None): Content for ``*SFC``.
        shp (Any | None): Content for ``*SHP``.
        svd (Any | None): Content for ``*SVD``.
        xsc (Any | None): Content for ``*XSC``.
        extra_sections (Mapping[str, Any] | None): Additional ``*KEY`` names
            mapped to section data.
        run_id (str | None): Identifier appended to output filenames.
        clean_before (bool): When ``True`` remove artefacts matching
            ``patterns`` before running.
        directory (Path | str | None): Working directory containing assets and
            outputs.
        patterns (Sequence[str] | None): Filename glob patterns passed to
            :func:`clean_outputs`.
        output_mode (Literal["files", "capture"]): Pass ``"capture"`` to
            populate the returned :class:`RunResult` with the in-memory
            optical-depth dataframe; otherwise rely on file outputs.
        enable_capture (bool | None): Override whether the shared-library
            execution enables optical-depth capture; defaults to mirroring
            ``output_mode == "capture"``.
        optical_levels (Sequence[float] | None): Output levels (km) required
            when ``output_mode`` is ``"capture"``.
        optical_spectrum_index (int): Spectral index used when retrieving
            captured optical depths (default ``1``).
        optical_match_tol (float): Matching tolerance (km) applied to the
            captured profile grid (default ``1e-6``).

    Returns:
        RunResult: Summary of the execution and optional in-memory output.

    Raises:
        ValueError: If required arguments are missing or mutually exclusive.
        FileNotFoundError: If the requested working directory does not exist.
    """

    root = Path(directory) if directory is not None else Path.cwd()
    if not root.exists():
        raise FileNotFoundError(f"Directory does not exist: {root}")

    flag_list = [code.upper() for code in flags]
    if not flag_list:
        raise ValueError("At least one flag must be supplied for *FLG.")

    uses_tab = "TAB" in flag_list
    if uses_tab and tab_dimensions is None:
        raise ValueError("TAB flag requires tab_dimensions (for the *DIM section).")
    if not uses_tab and tangent is None:
        raise ValueError("Tangent geometry must be supplied when TAB flag is absent.")
    if uses_tab and tangent is not None:
        raise ValueError("Provide either tangent or tab_dimensions, not both.")
    if not isinstance(header, (str, Sequence)):
        raise ValueError("header must be a string or sequence of strings.")

    sections = _compose_driver_sections(
        header=header,
        flags=flag_list,
        spectral=spectral,
        gases=gases,
        atmosphere=atmosphere,
        tangent=tangent,
        tab_dimensions=tab_dimensions,
        cia=cia,
        fin=fin,
        fov=fov,
        grd=grd,
        hit=hit,
        ils=ils,
        jac=jac,
        lev=lev,
        lut=lut,
        nte=nte,
        obs=obs,
        out=out,
        phy=phy,
        rej=rej,
        sfc=sfc,
        shp=shp,
        svd=svd,
        xsc=xsc,
        extra_sections=extra_sections,
        uses_tab=uses_tab,
    )

    return run_rfm(
        run_id=run_id,
        clean_before=clean_before,
        directory=root,
        patterns=patterns,
        driver_lines=sections,
        output_mode=output_mode,
        enable_capture=enable_capture,
        optical_levels=optical_levels,
        optical_spectrum_index=optical_spectrum_index,
        optical_match_tol=optical_match_tol,
    )


def rfm_main(
    *,
    configuration: Mapping[str, Any] | None,
    driver_inputs: Mapping[str, Any],
    levels: Sequence[float],
    rfm_out_fldr: Path | str,
) -> RunResult:
    """Run the RFM model using dictionary-based configuration and inputs.

    Args:
        configuration (Mapping[str, Any] | None): Optional overrides that
            control execution, including ``driver_path``, ``generate_driver``,
            ``clean_before``, ``run_id``, ``patterns``, ``optical_spectrum_index``,
            ``optical_match_tol``, ``output_mode``, ``capture_files_content``, and
            ``verbose``. When supplied, ``output_mode`` must be ``"capture"`` or
            ``"files"``.
        driver_inputs (Mapping[str, Any]): Structured representation of driver
            sections mirroring :func:`run_rfm_with_parameters`.
        levels (Sequence[float]): Altitude levels (km) used when optical-depth
            capture is enabled via the ``OPT`` and ``LEV`` flags.
        rfm_out_fldr (Path | str): Directory in which the run should take place
            and where RFM writes its outputs.

    Returns:
        RunResult: The run status plus optional payloads. ``output`` holds the
            in-memory file map for capture runs (``None`` for file-mode runs),
            while ``output_df`` is populated with the optical-depth dataframe
            when ``OPT`` and ``LEV`` are active under capture mode.

    Raises:
        FileNotFoundError: If the requested output directory does not exist.
        NotADirectoryError: If ``rfm_out_fldr`` resolves to a file.
        ValueError: If mandatory sections are missing, incompatible geometry
            inputs are provided, optical-depth capture requirements are not
            satisfied, or ``output_mode`` is invalid for the provided flags.
        RuntimeError: If optical-depth capture is expected but a dataframe is
            not returned.

    Notes:
        ``capture_files_content`` controls whether capture mode records the
        generated ``*.asc``/log artefacts in memory (while deleting them on
        disk). When set to ``False`` the run still suppresses file creation but
        omits the in-memory payload, leaving only ``output_df`` (when
        applicable).
    """

    config = dict(configuration or {})
    driver_path_value = config.get("driver_path")
    generate_driver = config.get("generate_driver", False)
    clean_before = config.get("clean_before", True)
    run_id = config.get("run_id")
    patterns = config.get("patterns")
    optical_spectrum_index = config.get("optical_spectrum_index", 1)
    optical_match_tol = config.get("optical_match_tol", 1e-6)
    verbose = config.get("verbose", True)
    capture_files_content = bool(
        config.get(
            "capture_files_content",
            True,
        )
    )

    def _snapshot_files(directory: Path) -> dict[str, tuple[int, int]]:
        """Capture file metadata for change detection between runs.

        Args:
            directory (Path): Location containing the files of interest.

        Returns:
            dict[str, tuple[int, int]]: Mapping of filenames to (size, mtime_ns)
                tuples.
        """
        snapshot: dict[str, tuple[int, int]] = {}
        for entry in directory.iterdir():
            if entry.is_file():
                metadata = entry.stat()
                snapshot[entry.name] = (
                    int(metadata.st_size),
                    int(metadata.st_mtime_ns),
                )
        return snapshot

    def _read_file_payload(path: Path) -> str | bytes:
        """Read file contents as text, falling back to bytes when needed.

        Args:
            path (Path): File to load from disk.

        Returns:
            str | bytes: Textual payload when UTF-8 decodable, otherwise the raw
                bytes.
        """
        try:
            return path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            return path.read_bytes()

    levels_list = [float(val) for val in levels] if levels is not None else []

    run_directory = Path(rfm_out_fldr).expanduser().resolve()
    if not run_directory.exists():
        raise FileNotFoundError(f"Output directory does not exist: {run_directory}")
    if not run_directory.is_dir():
        raise NotADirectoryError(
            f"Output directory is not a directory: {run_directory}"
        )

    def _collect_file_payload(
        before_snapshot: dict[str, tuple[int, int]] | None,
        *,
        delete_after: bool = False,
    ) -> tuple[dict[str, str | bytes], dict[str, tuple[int, int]]]:
        ignore_names = {"rfm.drv"}
        after_snapshot = _snapshot_files(run_directory)
        payload: dict[str, str | bytes] = {}
        for name, meta in after_snapshot.items():
            if name in ignore_names:
                continue
            if (
                before_snapshot is None
                or name not in before_snapshot
                or before_snapshot[name] != meta
            ):
                path = run_directory / name
                payload[name] = _read_file_payload(path)
                if delete_after:
                    try:
                        path.unlink()
                    except OSError:
                        pass
        return payload, after_snapshot

    base_kwargs: dict[str, Any] = dict(driver_inputs)
    flags = base_kwargs.get("flags")
    if not flags:
        raise ValueError("driver_inputs must define the 'flags' section.")
    flag_list = [code.upper() for code in flags]
    base_kwargs["flags"] = flag_list

    has_opt = "OPT" in flag_list
    has_lev = "LEV" in flag_list
    capture_flags_active = has_opt and has_lev
    output_flags_present = any(flag in OUTPUT_FLAG_CODES_SET for flag in flag_list)

    requested_output_mode = config.get("output_mode")
    if requested_output_mode is None:
        output_mode = "capture" if capture_flags_active else "files"
    else:
        output_mode = str(requested_output_mode).lower()
        if output_mode not in {"files", "capture"}:
            raise ValueError(
                f"Unsupported output_mode '{requested_output_mode}'. Expected 'files' or 'capture'."
            )

    capture_requested = output_mode == "capture"
    optical_capture_active = capture_requested and capture_flags_active
    file_payload_required = capture_requested and output_flags_present

    if optical_capture_active and len(levels_list) < 2:
        raise ValueError(
            "At least two levels are required when OPT and LEV flags are enabled."
        )

    uses_tab = "TAB" in flag_list
    if levels_list:
        base_kwargs["lev"] = tuple(str(val) for val in levels_list)

    for required_key in ("header", "spectral", "gases", "atmosphere"):
        if required_key not in base_kwargs:
            raise ValueError(f"driver_inputs must define the '{required_key}' section.")

    if uses_tab and base_kwargs.get("tab_dimensions") is None:
        raise ValueError("TAB flag requires 'tab_dimensions' within driver_inputs.")
    if not uses_tab and base_kwargs.get("tangent") is None:
        raise ValueError(
            "driver_inputs must include 'tangent' unless TAB mode is used."
        )

    file_snapshot_before: dict[str, tuple[int, int]] | None = None
    if capture_requested:
        file_snapshot_before = _snapshot_files(run_directory)

    driver_path: Path | None = None
    if driver_path_value is not None:
        driver_path = Path(driver_path_value).expanduser().resolve()
        if generate_driver or not driver_path.exists():
            sections = _compose_driver_sections(
                header=base_kwargs["header"],
                flags=flag_list,
                spectral=base_kwargs["spectral"],
                gases=base_kwargs["gases"],
                atmosphere=base_kwargs["atmosphere"],
                tangent=base_kwargs.get("tangent"),
                tab_dimensions=base_kwargs.get("tab_dimensions"),
                cia=base_kwargs.get("cia"),
                fin=base_kwargs.get("fin"),
                fov=base_kwargs.get("fov"),
                grd=base_kwargs.get("grd"),
                hit=base_kwargs.get("hit"),
                ils=base_kwargs.get("ils"),
                jac=base_kwargs.get("jac"),
                lev=base_kwargs.get("lev"),
                lut=base_kwargs.get("lut"),
                nte=base_kwargs.get("nte"),
                obs=base_kwargs.get("obs"),
                out=base_kwargs.get("out"),
                phy=base_kwargs.get("phy"),
                rej=base_kwargs.get("rej"),
                sfc=base_kwargs.get("sfc"),
                shp=base_kwargs.get("shp"),
                svd=base_kwargs.get("svd"),
                xsc=base_kwargs.get("xsc"),
                extra_sections=base_kwargs.get("extra_sections"),
                uses_tab=uses_tab,
            )
            driver_path.parent.mkdir(parents=True, exist_ok=True)
            driver_path.write_text("\n".join(sections) + "\n", encoding="utf-8")

    def _execute_run(
        *,
        mode: Literal["files", "capture"],
        enable_capture_flag: bool | None,
        include_levels: bool,
        clean: bool,
    ) -> RunResult:
        run_kwargs: dict[str, Any] = dict(
            run_id=run_id,
            clean_before=clean,
            directory=run_directory,
            patterns=patterns,
            output_mode=mode,
            optical_spectrum_index=optical_spectrum_index,
            optical_match_tol=optical_match_tol,
        )
        if enable_capture_flag is not None:
            run_kwargs["enable_capture"] = enable_capture_flag
        if include_levels and levels_list:
            run_kwargs["optical_levels"] = levels_list
        if driver_path is not None:
            return run_rfm(driver_path=driver_path, **run_kwargs)
        return run_rfm_with_parameters(
            **base_kwargs,
            **run_kwargs,
        )

    run_result = _execute_run(
        mode=output_mode,
        enable_capture_flag=optical_capture_active,
        include_levels=optical_capture_active,
        clean=clean_before,
    )

    def _run_file_mode_capture() -> dict[str, str | bytes]:
        before = _snapshot_files(run_directory)
        _execute_run(
            mode="files",
            enable_capture_flag=False,
            include_levels=False,
            clean=False,
        )
        payload, _ = _collect_file_payload(before, delete_after=True)
        return payload

    if verbose:
        print(f"RFM run status: {run_result.status}")

    if optical_capture_active:
        status, n_levels, n_points = rfm_py.rfm_get_optical_grid_size(
            ispc=optical_spectrum_index
        )
        if verbose:
            print(
                f"Captured grid dimensions (status={status}): "
                f"{n_levels} levels, {n_points} spectral points"
            )

        df = run_result.output_df
        if df is None and isinstance(run_result.output, pd.DataFrame):
            df = run_result.output
        if df is None:
            raise RuntimeError(
                "Expected optical-depth dataframe when OPT and LEV are enabled."
            )
        run_result.output_df = df
        if verbose:
            dod_columns = [col for col in df.columns if col.startswith("dOD_")]
            if dod_columns:
                sample_vector = df[dod_columns].iloc[0].to_numpy(dtype=float)
                print(
                    f"\nSample layer optical depth vector "
                    f"(first layer, {sample_vector.size} points):"
                )
                print(np.array2string(sample_vector, threshold=10, max_line_width=120))
            print(df)
    else:
        run_result.output_df = None

    if capture_requested and capture_files_content:
        payload, _ = _collect_file_payload(
            file_snapshot_before,
            delete_after=True,
        )
        if optical_capture_active and file_payload_required:
            extra_payload = _run_file_mode_capture()
            payload.update(extra_payload)
        run_result.output = payload
        if verbose:
            if payload:
                print("Captured file outputs:")
                for name, content in payload.items():
                    size = len(content) if hasattr(content, "__len__") else 0
                    kind = "binary" if isinstance(content, bytes) else "text"
                    print(f"- {name} ({kind}, {size} bytes)")
            else:
                print("No output files were generated by the run.")
    elif capture_requested:
        _collect_file_payload(
            file_snapshot_before,
            delete_after=True,
        )
        if verbose:
            print(
                "capture_files_content disabled; outputs removed from disk without capture."
            )
        run_result.output = None
    else:
        if verbose:
            print("Outputs written to disk; run_result.output not populated.")
        run_result.output = None

    return run_result


def _compose_driver_sections(
    *,
    header: str | Sequence[str],
    flags: Sequence[str],
    spectral: Any,
    gases: Any,
    atmosphere: Any,
    tangent: Any | None,
    tab_dimensions: Any | None,
    cia: Any | None,
    fin: Any | None,
    fov: Any | None,
    grd: Any | None,
    hit: Any | None,
    ils: Any | None,
    jac: Any | None,
    lev: Any | None,
    lut: Any | None,
    nte: Any | None,
    obs: Any | None,
    out: Any | None,
    phy: Any | None,
    rej: Any | None,
    sfc: Any | None,
    shp: Any | None,
    svd: Any | None,
    xsc: Any | None,
    extra_sections: Mapping[str, Any] | None,
    uses_tab: bool,
) -> List[str]:
    """Materialise the ``rfm.drv`` driver table from structured inputs.

    Args:
        header (str | Sequence[str]): Records for the ``*HDR`` section.
        flags (Sequence[str]): Normalised RFM flag codes for ``*FLG``.
        spectral (Any): Payload describing spectral ranges for ``*SPC``.
        gases (Any): Entries for the ``*GAS`` section.
        atmosphere (Any): Atmospheric profile definitions for ``*ATM``.
        tangent (Any | None): Geometry definition for ``*TAN`` when ``uses_tab``
            is ``False``.
        tab_dimensions (Any | None): Tabulation axis definition for ``*DIM``
            when ``uses_tab`` is ``True``.
        cia (Any | None): Optional ``*CIA`` section records.
        fin (Any | None): Optional fine-mesh sampling (``*FIN``) entries.
        fov (Any | None): Optional field-of-view (``*FOV``) entries.
        grd (Any | None): Optional irregular spectral grid (``*GRD``) entries.
        hit (Any | None): Optional line-transition assignments (``*HIT``).
        ils (Any | None): Optional instrument line shape (``*ILS``) entries.
        jac (Any | None): Optional Jacobian target records (``*JAC``).
        lev (Any | None): Optional altitude levels (``*LEV``).
        lut (Any | None): Optional lookup table references (``*LUT``).
        nte (Any | None): Optional non-LTE dataset references (``*NTE``).
        obs (Any | None): Optional observer geometry entries (``*OBS``).
        out (Any | None): Optional output overrides (``*OUT``).
        phy (Any | None): Optional physical constant overrides (``*PHY``).
        rej (Any | None): Optional line rejection thresholds (``*REJ``).
        sfc (Any | None): Optional surface boundary conditions (``*SFC``).
        shp (Any | None): Optional molecule-specific line shapes (``*SHP``).
        svd (Any | None): Optional SVD table references (``*SVD``).
        xsc (Any | None): Optional cross-section data references (``*XSC``).
        extra_sections (Mapping[str, Any] | None): Additional driver sections
            keyed by their ``*KEY`` identifiers.
        uses_tab (bool): Indicates whether the ``TAB`` flag is active and the
            ``*DIM`` section should be emitted instead of ``*TAN``.

    Returns:
        list[str]: Lines that form a complete ``rfm.drv`` file terminated with
            ``*END``.

    Raises:
        ValueError: If required sections are empty, mutually exclusive geometry
            sections are provided together, or the same section appears more
            than once.
    """

    def ensure_non_empty(name: str, data: List[str]) -> List[str]:
        if not data:
            raise ValueError(f"{name} section requires at least one entry.")
        return data

    sections: list[str] = []
    seen_keys: set[str] = set()

    def append_section(key: str, records: List[str], required: bool = False) -> None:
        if required:
            ensure_non_empty(key, records)
        if not records:
            return
        if key in seen_keys:
            raise ValueError(f"Section {key} supplied more than once.")
        seen_keys.add(key)
        sections.append(key)
        sections.extend(f"  {record.strip()}" for record in records)

    append_section("*HDR", _normalize_section_data(header), required=True)
    append_section("*FLG", [" ".join(flags)], required=True)
    append_section(
        "*SPC", ensure_non_empty("*SPC", _normalize_section_data(spectral)), True
    )
    append_section(
        "*GAS", ensure_non_empty("*GAS", _normalize_section_data(gases)), True
    )
    append_section(
        "*ATM", ensure_non_empty("*ATM", _normalize_section_data(atmosphere)), True
    )

    if uses_tab:
        append_section(
            "*DIM", ensure_non_empty("*DIM", _normalize_section_data(tab_dimensions))
        )
    else:
        append_section(
            "*TAN", ensure_non_empty("*TAN", _normalize_section_data(tangent)), True
        )

    optional_map = {
        "*CIA": cia,
        "*FIN": fin,
        "*FOV": fov,
        "*GRD": grd,
        "*HIT": hit,
        "*ILS": ils,
        "*JAC": jac,
        "*LEV": lev,
        "*LUT": lut,
        "*NTE": nte,
        "*OBS": obs,
        "*OUT": out,
        "*PHY": phy,
        "*REJ": rej,
        "*SFC": sfc,
        "*SHP": shp,
        "*SVD": svd,
        "*XSC": xsc,
    }

    for key in (
        "*CIA",
        "*FIN",
        "*FOV",
        "*GRD",
        "*HIT",
        "*ILS",
        "*JAC",
        "*LEV",
        "*LUT",
        "*NTE",
        "*OBS",
        "*OUT",
        "*PHY",
        "*REJ",
        "*SFC",
        "*SHP",
        "*SVD",
        "*XSC",
    ):
        data = optional_map[key]
        if data is not None:
            append_section(key, _normalize_section_data(data))

    if extra_sections:
        for raw_key, payload in extra_sections.items():
            key = raw_key if raw_key.startswith("*") else f"*{raw_key.upper()}"
            append_section(key, _normalize_section_data(payload))

    sections.append("*END")
    return sections


def _normalize_section_data(data: Any) -> List[str]:
    """Convert flexible section payloads into driver table records.

    Args:
        data (Any): Arbitrary user-supplied representation of section content,
            including scalars, nested sequences, dataclass helpers, mappings,
            and filesystem paths.

    Returns:
        list[str]: Canonicalised records ready to be written under a section.
    """

    if data is None:
        return []

    if isinstance(data, (SpectralRange, SpectralFile, SectionLine)):
        return [data.as_record()]

    if isinstance(data, Mapping):
        lines: list[str] = []
        for key, value in data.items():
            lines.append(f"{key} = {_format_token(value)}")
        return lines

    if isinstance(data, Path):
        return [str(data)]

    if isinstance(data, str):
        return [data.strip()]

    if isinstance(data, Sequence):
        lines: list[str] = []
        for entry in data:
            if isinstance(entry, (SpectralRange, SpectralFile, SectionLine)):
                lines.append(entry.as_record())
            elif isinstance(entry, Mapping):
                lines.extend(_normalize_section_data(entry))
            elif isinstance(entry, Path):
                lines.append(str(entry))
            elif isinstance(entry, str):
                lines.append(entry.strip())
            elif isinstance(entry, Sequence):
                tokens = " ".join(_format_token(token) for token in entry)
                lines.append(tokens)
            else:
                lines.append(_format_token(entry))
        return lines

    return [_format_token(data)]


def _format_token(value: Any) -> str:
    """Render a single token suitable for inclusion in a driver record.

    Args:
        value (Any): Token supplied as part of a section entry.

    Returns:
        str: String representation matching the formatting rules expected by
            the RFM driver parser.
    """

    if isinstance(value, (SpectralRange, SpectralFile, SectionLine)):
        return value.as_record()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float):
        return f"{value:g}"
    if isinstance(value, (int,)):
        return str(value)
    return str(value)


@contextmanager
def _temporary_cwd(target: Path) -> Iterable[None]:
    """Temporarily change the working directory within a context manager.

    Args:
        target (pathlib.Path): Directory to switch into for the duration of the
            context.

    Yields:
        None: Control to the caller while the process operates inside
            ``target``.
    """

    original = Path.cwd()
    if original == target:
        yield
        return

    os.chdir(target)
    try:
        yield
    finally:
        os.chdir(original)


def _run_rfm_worker_entry(args_path: Path, result_path: Path) -> None:
    """Entry point used by the helper subprocess."""

    with Path(args_path).open("rb") as handle:
        payload = pickle.load(handle)

    try:
        result = _run_rfm_impl(*payload["args"], **payload["kwargs"])
    except Exception:
        outcome = ("error", traceback.format_exc())
    else:
        outcome = ("ok", result)

    with Path(result_path).open("wb") as handle:
        pickle.dump(outcome, handle)


def _cli() -> None:
    if len(sys.argv) == 4 and sys.argv[1] == "--rfm-worker":
        _run_rfm_worker_entry(Path(sys.argv[2]), Path(sys.argv[3]))
        return
    raise SystemExit("rfm_helper is not intended to be executed directly.")


if __name__ == "__main__":
    _cli()
