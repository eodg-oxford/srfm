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


def get_ri_filepathname(input_string):
    """Maps an input string to a file path.

    Args:
        input_string: either an ARIA filename or any of "ash", "ice", "sulphuric acid"

    Returns:
        path: absolute file path of the refractive indices file.

    """

    if not isinstance(input_string, str):
        raise TypeError("composition must be a string.")

    if input_string == "ash":
        input_string = "eyjafjallajokull-ash_Reed.ri"

    if input_string == "ice":
        input_string = "ICE_Warren_2008.ri"

    if input_string == "sulphuric acid":
        input_string = "H2SO4_75_Palmer_1975.ri"

    with as_file(files("srfm.data") / "ARIA") as path:

        # Recursively search for the file within the ARIA directory tree
        for root, dirs, fls in os.walk(path):
            if input_string in fls:
                return os.path.join(root, input_string)  # Return the absolute file path

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
        """Selects requested data.

        Args:
            wave: User input wave.
            mode: Output mode. Can be "wavelength" or "wavenumber".
                out_of_range: Defines how out-of-range values are handled.
                Could be:

                - error: Error is raised.
                - clip: Data is truncated.
                - nan: Data is interpolated.

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
        """Reads the refractive index data for a given ri file.

        Interpolates to the wave values if provided, or returns full-resolution data if wave is None.

        Args:
            composition (str):
                Any of:
                    - an ARIA file name
                    - an ARIA file group from the list below (NOT IMPLEMENTED)
                    - "ZASETSKY", keyword temperature needs to be set so the retruned reractive indexx of water is interpolated to that temperature
                    - a generic name from the list below:
                        - "ash"
                        - "ice"
                        - "sulphuric acid"
                        - "water". (NOT IMPLEMENTED)
            wave (list or array, optional): The target wavelengths or wavenumbers to interpolate to. If None, returns data at full resolution.
            mode (str): 'wavelength' for wave in µm or 'wavenumber' for wave in cm⁻¹.
            out_of_range (str): Behavior for out-of-range values: 'error', 'clip', or 'nan'.

        Returns:
            tuple: Two or three arrays:
                   - If wave is None: (x_data, n, k), where x_data is `wavl` or `wavn`.
                   - If wave is defined: (n, k), the interpolated real and imaginary parts of the refractive index.
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
    """Finds .ri files in requested path.

    Args:
        - ARIA_path: Requested path.

    Returns:
        list: Path to .ri files.
    """
    import os

    ri_files = []

    for root, dirs, files in os.walk(ARIA_path):
        for file in files:
            if file.endswith(".ri"):
                ri_files.append(os.path.join(root, file))
    return ri_files


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
