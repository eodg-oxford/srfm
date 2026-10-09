"""Used to access the ARIA database.

For details see https://eodg.atm.ox.ac.uk/ARIA/.

- Name: ARIA_module
- Parent package: srfm
- Author: Don Grainger
- Contributors: Antonin Knizek
- Date: 24 January 2025
"""

import os
import numpy as np
from importlib.resources import files, as_file

from ._aria_aliases import LEGACY_RI_FILENAMES


_GENERIC_RI_FILENAMES = {
    "ash": "eyjafjallajokull_ash_58.5%SiO2_Reed_2018.ri",
    "ice": "ice_266K_Warren_2008.ri",
    "sulphuric acid": "H2SO4_75%_300K_Palmer_1975.ri",
}


def get_ri_filepathname(input_string):
    """Resolve a composition or filename in the bundled ARIA database.

    Current and legacy filenames are matched exactly, including case. Legacy
    names select the corresponding renamed dataset, preserving the original
    sample, temperature, concentration, and crystal orientation.

    Args:
        input_string (str): Current or legacy ARIA basename, or one of
            ``"ash"``, ``"ice"``, and ``"sulphuric acid"``.

    Returns:
        str: Absolute path to the refractive-index file under ``srfm/data/ARIA``.

    Raises:
        TypeError: If the composition is not a string.
        FileNotFoundError: If no bundled dataset matches the composition.
    """

    if not isinstance(input_string, str):
        raise TypeError("composition must be a string.")

    requested_filename = _GENERIC_RI_FILENAMES.get(input_string, input_string)
    requested_filename = LEGACY_RI_FILENAMES.get(
        requested_filename, requested_filename
    )

    with as_file(files("srfm.data") / "ARIA") as path:

        # Recursively search for the file within the ARIA directory tree
        for directory, _, filenames in os.walk(path):
            if requested_filename in filenames:
                return os.path.join(directory, requested_filename)

    raise FileNotFoundError(
        f"Refractive-index file '{input_string}' was not found in bundled ARIA data."
    )


class ReadError(Exception):
    """Custom exception raised for errors in reading .ri files."""

    pass


class RI:
    """Class representing refractive index data."""

    expected_header_names = [
        "FORMAT",
        "DESCRIPTION",
        "DISTRIBUTEDBY",
        "SUBSTANCE",
        "SAMPLEFORM",
        "TEMPERATURE",
        "CONCENTRATION",
        "REFERENCE",
        "DOI",
        "SOURCE",
        "CONTACT",
        "COMMENT",
    ]
    expected_column_names = ["wavl", "wavn", "n", "dn", "k", "dk"]

    def __init__(self):
        """Initialize empty mappings for refractive-index headers and data."""
        self.header = {}
        self.data = {}

    def read(self, filepathname):
        """Reads and parses an .ri file into the object's attributes.

        Args:
            filepathname: Refractive index filepath.
        """

        self.header = {}
        self.data = {}
        try:
            with open(filepathname, "r", encoding="utf-8") as handle:
                lines = [line.strip() for line in handle]
        except (OSError, UnicodeError) as exc:
            raise ReadError(f"Could not read refractive-index file: {filepathname}") from exc

        while lines and not lines[-1]:
            lines.pop()

        header_lines: list[str] = []
        data_lines: list[str] = []
        data_started = False
        for line in lines:
            if not line:
                continue
            if line.startswith("#"):
                if data_started:
                    raise ReadError(
                        f"Incorrectly formatted file ({filepathname}): Header not contiguous."
                    )
                header_lines.append(line)
            else:
                data_started = True
                data_lines.append(line)

        if not header_lines:
            raise ReadError(f"Incorrectly formatted file ({filepathname}): No header.")
        if not data_lines:
            raise ReadError(f"Incorrectly formatted file ({filepathname}): No data.")

        for raw_line in header_lines:
            line = raw_line[1:].strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            tag_name, tag_content = (part.strip() for part in line.split("=", 1))
            tag_name = tag_name.upper()
            if tag_name not in self.expected_header_names:
                continue
            if tag_name in self.header:
                tag_content = self.header[tag_name] + " " + tag_content
            self.header[tag_name] = tag_content

        if "FORMAT" not in self.header:
            raise ReadError(
                f"Incorrectly formatted file ({filepathname}): No FORMAT tag in header."
            )

        column_labels = [item.lower() for item in self.header["FORMAT"].split()]
        if not column_labels:
            raise ReadError(
                f"Incorrectly formatted file ({filepathname}): Empty FORMAT tag."
            )
        unknown = [item for item in column_labels if item not in self.expected_column_names]
        if unknown:
            raise ReadError(
                f"Incorrectly formatted file ({filepathname}): Unknown FORMAT columns: "
                + ", ".join(unknown)
            )
        if len(column_labels) != len(set(column_labels)):
            raise ReadError(
                f"Incorrectly formatted file ({filepathname}): Duplicate FORMAT columns."
            )
        if "n" not in column_labels or "k" not in column_labels:
            raise ReadError(
                f"Incorrectly formatted file ({filepathname}): FORMAT requires n and k."
            )
        if "wavl" not in column_labels and "wavn" not in column_labels:
            raise ReadError(
                f"Incorrectly formatted file ({filepathname}): FORMAT requires wavl or wavn."
            )

        self.data = {label: [] for label in column_labels}
        for row_number, raw_line in enumerate(data_lines, start=1):
            fields = raw_line.split()
            if len(fields) != len(column_labels):
                raise ReadError(
                    f"Incorrectly formatted file ({filepathname}): Data row {row_number} "
                    f"has {len(fields)} columns; expected {len(column_labels)}."
                )
            try:
                values = [float(field) for field in fields]
            except ValueError as exc:
                raise ReadError(
                    f"Incorrectly formatted file ({filepathname}): Non-numeric data "
                    f"in row {row_number}."
                ) from exc
            for label, value in zip(column_labels, values):
                self.data[label].append(value)

        if "wavn" not in self.data:
            self.data["wavn"] = [
                10000.0 / value if value != 0 else float("nan")
                for value in self.data["wavl"]
            ]
        if "wavl" not in self.data:
            self.data["wavl"] = [
                10000.0 / value if value != 0 else float("nan")
                for value in self.data["wavn"]
            ]

    def select(self, wave=None, mode="wavelength", out_of_range="error"):
        """Return full-resolution indices or interpolate to a spectral grid.

        Args:
            wave (array-like, optional): Target spectral coordinates. If None,
                return the stored grid and indices without interpolation.
            mode (str): ``"wavelength"`` for micrometres or ``"wavenumber"``
                for inverse centimetres.
            out_of_range (str): ``"error"`` raises for coordinates outside the
                data range, ``"clip"`` uses the nearest endpoint, and ``"nan"``
                returns NaN outside the range while interpolating inside it.

        Returns:
            tuple: ``(grid, n, k)`` when wave is None, otherwise ``(n, k)``.

        Raises:
            ValueError: If the mode or range policy is invalid, the coordinate
                data is missing, or a requested coordinate violates ``"error"``.
        """

        if out_of_range not in {"error", "clip", "nan"}:
            raise ValueError(
                "Invalid value for out_of_range. Use 'error', 'clip', or 'nan'."
            )

        # Determine which data to use
        if mode == "wavelength":
            x_data = self.data.get("wavl")
            if x_data is None:
                raise ValueError("No wavelength (wavl) data available.")
        elif mode == "wavenumber":
            x_data = self.data.get("wavn")
            if x_data is None:
                raise ValueError("No wavenumber (wavn) data available.")
        else:
            raise ValueError("Invalid mode. Use 'wavelength' or 'wavenumber'.")

        # If no wave are provided, return full-resolution data
        if wave is None:
            return np.array(x_data), np.array(self.data["n"]), np.array(self.data["k"])

        wave = np.atleast_1d(np.asarray(wave, dtype=float))

        # Handle out-of-range values
        min_x, max_x = min(x_data), max(x_data)
        if min(wave) < min_x or max(wave) > max_x:
            if out_of_range == "error":
                raise ValueError(
                    f"Requested values are outside the valid range. "
                    f"Valid range: {min_x} to {max_x}."
                )
            elif out_of_range == "clip":
                wave = np.clip(wave, min_x, max_x)
            elif out_of_range == "nan":
                interpolated_n = np.full(len(wave), np.nan)
                interpolated_k = np.full(len(wave), np.nan)
                valid_indices = (wave >= min_x) & (wave <= max_x)
                if np.all(np.diff(x_data) > 0):  # Ascending order check
                    interpolated_n[valid_indices] = np.interp(
                        wave[valid_indices], x_data, self.data["n"]
                    )
                    interpolated_k[valid_indices] = np.interp(
                        wave[valid_indices], x_data, self.data["k"]
                    )
                else:  # Sort data for interpolation
                    sorted_indices = np.argsort(x_data)
                    x_data_sorted = np.array(x_data)[sorted_indices]
                    n_sorted = np.array(self.data["n"])[sorted_indices]
                    k_sorted = np.array(self.data["k"])[sorted_indices]
                    interpolated_n[valid_indices] = np.interp(
                        wave[valid_indices], x_data_sorted, n_sorted
                    )
                    interpolated_k[valid_indices] = np.interp(
                        wave[valid_indices], x_data_sorted, k_sorted
                    )
                return interpolated_n, interpolated_k

        # Perform interpolation if wave are provided
        if np.all(np.diff(x_data) > 0):  # Ascending order check
            interpolated_n = np.interp(wave, x_data, self.data["n"])
            interpolated_k = np.interp(wave, x_data, self.data["k"])
        else:  # Sort data for interpolation
            sorted_indices = np.argsort(x_data)
            x_data_sorted = np.array(x_data)[sorted_indices]
            n_sorted = np.array(self.data["n"])[sorted_indices]
            k_sorted = np.array(self.data["k"])[sorted_indices]
            interpolated_n = np.interp(wave, x_data_sorted, n_sorted)
            interpolated_k = np.interp(wave, x_data_sorted, k_sorted)
        return interpolated_n, interpolated_k

    def load_refractive_indices(
        self, composition, wave=None, mode="wavelength", out_of_range="error"
    ):
        """Load a bundled ARIA dataset and optionally interpolate its indices.

        Args:
            composition (str): Current or legacy ARIA basename, or ``"ash"``,
                ``"ice"``, or ``"sulphuric acid"``. See
                :func:`get_ri_filepathname` for filename compatibility.
            wave (array-like, optional): Target spectral coordinates. If None,
                return data at full resolution.
            mode (str): ``"wavelength"`` for micrometres or ``"wavenumber"``
                for inverse centimetres.
            out_of_range (str): Range policy: ``"error"``, ``"clip"``, or
                ``"nan"``. See :meth:`select`.

        Returns:
            tuple: ``(grid, n, k)`` when wave is None, otherwise interpolated
            ``(n, k)``. The extinction coefficient ``k`` retains ARIA's sign;
            the Mie loader converts it to the ``n - ik`` convention.

        Raises:
            TypeError: If composition is not a string.
            FileNotFoundError: If the dataset is absent from bundled ARIA.
            ReadError: If the file cannot be read or parsed.
            ValueError: If spectral selection fails.
        """

        filepathname = get_ri_filepathname(composition)

        self.read(filepathname)

        if wave is None:
            w, n, k = self.select(wave=wave, mode=mode)
            return w, n, k
        else:
            n, k = self.select(wave=wave, mode=mode, out_of_range=out_of_range)
            return n, k


def find_ri_files(ARIA_path):
    """Find refractive-index files recursively beneath a directory.

    Args:
        ARIA_path (str or path-like): Root directory to search.

    Returns:
        list[str]: Paths to files whose names end in ``.ri``. Paths retain
        the absolute or relative form of the supplied root directory.
    """
    refractive_index_paths = []
    for directory, _, filenames in os.walk(ARIA_path):
        for filename in filenames:
            if filename.endswith(".ri"):
                refractive_index_paths.append(os.path.join(directory, filename))
    return refractive_index_paths


def read_ri_file(filepathname, wave=None, mode="wavelength", out_of_range="error"):
    """Reads the refractive index data for a given ri file.

    Interpolates to the wave values if provided,
    or returns full-resolution data if wave is None.

    Args:
        filepathname: an ARIA filename
        wave (list or array, optional): The target wavelengths or wavenumbers to interpolate to. If None, returns data at full resolution.
        mode (str): 'wavelength' for wave in µm or 'wavenumber' for wave in cm⁻¹.
        out_of_range (str): Behavior for out-of-range values: 'error', 'clip', or 'nan'.

    Returns:
        Two or three arrays:
            - If wave is None: (x_data, n, k), where x_data is `wavl` or `wavn` depending on mode.
            - If wave is defined: (n, k), the interpolated real and imaginary parts of the refractive index.
    """
    #    from ARIA_module import RI  # Import within the function
    ri_object = RI()
    ri_object.read(filepathname)
    if "n" in ri_object.data and "k" in ri_object.data:
        if wave is None:
            w, n, k = ri_object.select(wave=wave, mode=mode)
            return w, n, k
        else:
            n, k = ri_object.select(wave=wave, mode=mode, out_of_range=out_of_range)
            return n, k
    else:
        print("Both N & k do not exist in this file")
