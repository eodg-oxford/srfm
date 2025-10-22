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
    
    if input_string == "ash":
        input_string = "eyjafjallajokull-ash_Reed.ri"
    
    if input_string == "ice":
        input_string = "ICE_Warren_2008.ri"   
    
    if input_string == "sulphuric acid":
        input_string = "H2SO4_75_Palmer_1975.ri"
    
    with as_file(files("srfm.data") / "ARIA" ) as path:
       
    # Recursively search for the file within the ARIA directory tree
        for root, dirs, fls in os.walk(path):
            if input_string in fls:
                return os.path.join(root, input_string)  # Return the absolute file path

    return f"Error: File '{input_string}' not found in ARIA directory tree." 


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

        with open(filepathname, "r") as f:
            t = f.readlines()
            t = [x.strip() for x in t]  # Strip whitespace from lines

            header_lines = 0
            data_lines = 0
            for line in t:
                if line.startswith("#"):
                    header_lines += 1
                    if data_lines > 0:
                        raise ReadError(
                            f"Incorrectly formatted file ({filepathname}): Header not contiguous."
                        )
                else:
                    data_lines += 1
            for i in range(1, data_lines):  # Ignore blank lines at the end
                if any(char.isdigit() for char in t[-i]):
                    break
                else:
                    data_lines -= 1

            if header_lines == 0:
                raise ReadError(
                    f"Incorrectly formatted file ({filepathname}): No header."
                )
            if data_lines == 0:
                raise ReadError(
                    f"Incorrectly formatted file ({filepathname}): No data."
                )

            # Parse headers
            for i in range(header_lines):
                line = t[i][1:]  # Strip leading '#'
                if line[0] != "#":
                    tag_name = line.split("=", 1)[0].strip().upper()
                    if tag_name not in self.expected_header_names:
                        print(
                            f'Unknown header tag "{tag_name}", so ignored (file: {filepathname})'
                        )
                        continue
                    try:
                        tag_content = line.split("=", 1)[1].strip()
                    except UnicodeDecodeError:
                        tag_content = line.split("=", 1)[1].strip()
                    if tag_name in self.header:
                        tag_content = self.header[tag_name] + " " + tag_content
                    self.header[tag_name] = tag_content

            if "FORMAT" not in self.header:
                raise ReadError(
                    f"Incorrectly formatted file ({filepathname}): No FORMAT tag in header."
                )

            column_labels = self.header["FORMAT"].split()
            column_labels = [x.strip().lower() for x in column_labels]
            for cl in column_labels:
                if cl not in self.expected_column_names:
                    print(
                        f'Unknown column name "{cl}", so ignored (file: {filepathname})'
                    )
                    continue
                self.data[cl] = []

            for l in range(header_lines, header_lines + data_lines):
                line = t[l].split()
                line = [x.strip() for x in line]
                for c in range(len(column_labels)):
                    if column_labels[c] in self.expected_column_names:
                        self.data[column_labels[c]].append(float(line[c]))

            # Add wavl & wavn columns if needed (wavl in micro-m, wavn in cm-1)
            if "wavn" not in self.data:
                self.data["wavn"] = [
                    float(10000) / x if x != 0 else float("nan")
                    for x in self.data["wavl"]
                ]
            if "wavl" not in self.data:
                self.data["wavl"] = [
                    float(10000) / x if x != 0 else float("nan")
                    for x in self.data["wavn"]
                ]
            col_lengths = [len(self.data[col]) for col in self.data]
            if len(set(col_lengths)) > 1:
                raise ReadError(
                    f"Incorrectly formatted file ({filepathname}): Data columns have different lengths."
                )

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
            else:
                raise ValueError(
                    "Invalid value for out_of_range. Use 'error', 'clip', or 'nan'."
                )

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
