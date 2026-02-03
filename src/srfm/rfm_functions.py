"""Provides functions that enable the user to work with RFM.

- Name: rfm_functions
- Parent package: srfm
- Author: Antonin Knizek
- Contributors:
- Date: 18 February 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import sys
import pandas as pd
from . import utilities as utils
import warnings


def read_output(filename):
    """Read output file from RFM into a dictionary.

    Args:
        filename (str): rfm output file filename.

    Returns:
        contents (dict): Dictionary containing RFM outputs. Contains the key "info",
            which is a dictionary with the description of other parameters.

    """
    f = open(filename, "r")
    f_lines = f.readlines()
    f.close()

    contents = {}
    contents["info"] = {}
    contents["info"][
        "header1"
    ] = "Spectrum type, ray nad and RFM version ID, format C80"
    contents["info"]["header2"] = "Text from *HDR section of Driver Table, format C80"
    contents["info"][
        "header3"
    ] = "Captions for next record, or additional info, format C80"
    contents["info"][
        "NPNT"
    ] = "No. spectral points in file, if<0, then unit = GHz, otherwise cm-1, format I"
    contents["info"]["WNO1"] = "Lower limit of spectrum, format D"
    contents["info"]["WNOD"] = "Spectral grid interval (0=irregular grid), format D"
    contents["info"]["WNO2"] = "Upper limit of spectrum, format D"
    contents["info"][
        "LABSPC"
    ] = "Spectral Label or type of spectrum, in single quotes, format C*"
    contents["info"]["WNO"] = "Wavenumber or Frequency of spectral point, format D"
    contents["info"]["SPC"] = "Spectral data value, format R/D"

    contents["header1"] = f_lines[0][1:].strip()
    contents["header2"] = f_lines[0][1:].strip()
    contents["header3"] = f_lines[0][1:].strip()

    flds = f_lines[3].strip().split()
    contents["NPNT"] = int(flds[0])
    contents["WNO1"] = float(flds[1])
    contents["WNOD"] = float(flds[2])
    contents["WNO2"] = float(flds[3])
    contents["LABSPC"] = (" ").join(flds[4:]).strip("'")
    if contents["WNOD"] > 0:  # signifies regular grid
        contents["WNO"] = (
            contents["WNO1"] + np.arange(contents["NPNT"]) * contents["WNOD"]
        )
        spc = [i.strip() for i in f_lines[4:]]
        spc = (" ").join(spc)
        spc = [float(ii) for ii in spc.split()]
        contents["SPC"] = np.asarray(spc)
    else:
        tot = [i.strip() for i in f_lines[4:]]
        tot = (" ").join(tot)
        tot = np.asarray([float(iii) for iii in tot.split()]).reshape(
            contents["NPNT"], 2
        )
        contents["WNO"] = tot[:, 0]
        contents["SPC"] = tot[:, 1]
    return contents


def read_output_prf(filename):  # read bits of internal profile output file prf.asc
    """This function reads the RFM output profile prf.asc.

    Args:
        filename (str): path, including filename, of the RFM prf.asc output file (and
            equivalents if that ever changes from the current RFM default).

    Returns:
        contents (dict): Dictionary that stores the read values. Keys are section
            headers from the prf.asc file.

    """
    print("Attempting to read prf.asc file.")
    f = open(filename, "r")
    f_lines = f.readlines()
    f.close()

    contents = {}  # dictionary to store that read values.

    # get section start lines
    sec_begin_lns = [f_lines.index(x) for x in f_lines if x.startswith("*")]
    for idx in range(len(sec_begin_lns) - 1):
        sec_lbl = f_lines[sec_begin_lns[idx]].strip(" \n*")
        if isinstance(sec_lbl, str) == False:
            raise TypeError("Section label must be a string.")

        sec_cont = [
            i.strip() for i in f_lines[sec_begin_lns[idx] + 1 : sec_begin_lns[idx + 1]]
        ]
        sec_cont = (" ").join(sec_cont)
        sec_cont = [float(ii) for ii in sec_cont.split()]
        contents[sec_lbl] = sec_cont

    print("Successfully read prf.asc file.")

    return contents


@utils.show_runtime
def get_rfm_optical_depths(fldr, levels):
    """Calculates layer optical depths from RFM output.

    This function is designed to calculate optical depths for layers in
    the atmosphere from RFM with the LEV and OPT flags set on.
    The basic idea is to load optical depth spectra and subtract the adjacent ones
    to get layer optical depth.

    Args:
        fldr (str): Path to RFM folder.
        levels (list): levels to derive the layers from. The RFM can be output at more
            levels than necessary (which happens when internal RFM levels don't match
            the required output levels, i.e. most of the time).

    Returns:
        df (Pandas.core.frame.Dataframe): Pandas dataframe with the calculated layer optical depths
            and atmospheric structure.

    Raises:
        ValueError: Raised when the pressure profile read from prf.asc and the number of
            files don't match. This is often invoked when some input parameters are
            changed, but RFM is not rerun.
        ValueError: Raised when wavenumber grids for adjacent levels don't match.
        ValueError: Raised when delta_OD (layer optical depth) and int_OD (integrated
            optical depth) don't match.
        ValueError: Raised when the length of the wavenumber grid and delta_OD array
            don't match.

    """

    # load relevant rfm output files, sort by modification time (i.e. altitude)
    fls = sorted(Path(fldr).iterdir())
    fls = [str(p.absolute()) for p in fls]
    fls = [
        i
        for i in fls
        if "up" in str(i)[str(i).rfind("/") :] and "opt" in str(i)[str(i).rfind("/") :]
    ]  # this is the list of OD output spectra

    prf = read_output_prf(f"{fldr}/prf.asc")  # read rfm output profile, type dict

    # create a subset of prf that matches the user-desired output levels
    lvls_idx = [prf["HGT [km]"].index(lvl) for lvl in prf["HGT [km]"] if lvl in levels]
    for pkey in prf.keys():
        prf[pkey] = [prf[pkey][i] for i in lvls_idx]

    if len(fls) != len(prf["PRE [mb]"]):
        raise ValueError(
            """length of pressure profile and the number of corresponding 
            optical depth files don't match. Have you perhaps changed some input 
            parameters of your model and not ran RFM?"""
        )

    delta_OD = []  # layers' optical depths

    u_p = []  # layers' top pressures [mbar]
    l_p = []  # layers' bottom pressures [mbar]
    avg_p = []  # layers' average pressures [mbar]

    u_h = []  # layers' top altitudes [km]
    l_h = []  # layers' bottom altitudes [km]
    avg_h = []  # layers' average altitudes [km]

    u_t = []  # layers' top temperatures [K]
    l_t = []  # layers' bottom temperatures [K]
    avg_t = []  # layers' average temperatures [K]

    # read layer altitude, pressure, optical depth
    for l in range(len(fls) - 1):
        upper_level = read_output(fls[l + 1])
        lower_level = read_output(fls[l])

        d_OD = lower_level["SPC"] - upper_level["SPC"]  # layer optical depths
        delta_OD.append(d_OD)

        u_p.append(prf["PRE [mb]"][l + 1])
        l_p.append(prf["PRE [mb]"][l])
        a_p = (prf["PRE [mb]"][l + 1] - prf["PRE [mb]"][l]) / 2 + prf["PRE [mb]"][
            l
        ]  # layer avg pressure
        avg_p.append(a_p)

        u_h.append(prf["HGT [km]"][l + 1])
        l_h.append(prf["HGT [km]"][l])
        a_h = (prf["HGT [km]"][l + 1] - prf["HGT [km]"][l]) / 2 + prf["HGT [km]"][
            l
        ]  # layer avg altitude
        avg_h.append(a_h)

        u_t.append(prf["TEM [K]"][l + 1])
        l_t.append(prf["TEM [K]"][l])
        a_t = (prf["TEM [K]"][l + 1] - prf["TEM [K]"][l]) / 2 + prf["TEM [K]"][
            l
        ]  # layer avg temperature
        avg_t.append(a_t)

    wnos_hi = upper_level["WNO"]
    wnos_lo = lower_level["WNO"]
    if not np.array_equal(wnos_hi, wnos_lo):
        raise ValueError("Wavenumber grids for adjacent levels don't match.")

    wnos = upper_level["WNO"]

    delta_OD = np.asarray(delta_OD)
    delta_OD = np.flipud(delta_OD)

    # calculate integrated optical depth at layer
    int_OD = []
    for ii in range(delta_OD.shape[0]):
        int_OD.append(sum(delta_OD[: ii + 1]))

    int_OD = np.asarray(int_OD)

    if len(delta_OD) != len(int_OD):
        raise ValueError("The length of delta_OD and int_OD don't match.")

    if delta_OD.shape[1] != len(wnos):
        raise ValueError("The length of delta_OD and wavenumber grid don't match.")

    # fold everything into a dataframe
    prf_df = pd.DataFrame()
    prf_df["layer no."] = range(len(avg_p))  # layer number, 0 = TOA
    prf_df.set_index("layer no.")

    pup = pd.DataFrame({"p_upper (mbar)": u_p[::-1]}, dtype="float")
    prf_df.loc[:, "p_upper (mbar)"] = pup

    plo = pd.DataFrame({"p_lower (mbar)": l_p[::-1]}, dtype="float")
    prf_df.loc[:, "p_lower (mbar)"] = plo

    pav = pd.DataFrame({"p_avg (mbar)": avg_p[::-1]}, dtype="float")
    prf_df.loc[:, "p_avg (mbar)"] = pav

    hup = pd.DataFrame({"h_upper (km)": u_h[::-1]}, dtype="float")
    prf_df.loc[:, "h_upper (km)"] = hup

    hlo = pd.DataFrame({"h_lower (km)": l_h[::-1]}, dtype="float")
    prf_df.loc[:, "h_lower (km)"] = hlo

    hav = pd.DataFrame({"h_avg (km)": avg_h[::-1]}, dtype="float")
    prf_df.loc[:, "h_avg (km)"] = hav

    tup = pd.DataFrame({"T_upper (K)": u_t[::-1]}, dtype="float")
    prf_df.loc[:, "T_upper (K)"] = tup

    tlo = pd.DataFrame({"T_lower (K)": l_t[::-1]}, dtype="float")
    prf_df.loc[:, "T_lower (K)"] = tlo

    tav = pd.DataFrame({"T_avg (K)": avg_t[::-1]}, dtype="float")
    prf_df.loc[:, "T_avg (K)"] = tav

    # check if there are too large temperature steps (10 K difference)
    for i in range(len(tav) - 1):
        if abs(tav["T_avg (K)"][i + 1] - tav["T_avg (K)"][i]) >= 10:
            warnings.warn(
                f"""Temperature step between {hav['h_avg (km)'][i+1]} and 
            {hav['h_avg (km)'][i]} larger than 10K, calculation may be 
            inaccurate in DISORT."""
            )
        else:
            pass

    dod_col_names = [f"dOD_{val:.4f}" for val in wnos]
    iod_col_names = [f"iOD_{val:.4f}" for val in wnos]

    dod_df = pd.DataFrame(delta_OD, columns=dod_col_names)
    iod_df = pd.DataFrame(int_OD, columns=iod_col_names)

    df = pd.concat([prf_df, dod_df, iod_df], axis=1, join="outer")

    return df


def compile_rfm(fldr):
    """Compiles RFM.

    Args:
        fldr (str): Path to RFM folder (absolute or relative).

    Raises:
        OSError: Raised when directory is invalid.

    """
    cwd = os.getcwd()
    try:
        os.chdir(f"{fldr}/source")
    except OSError:
        print("Invalid directory.")
        return
    try:
        os.system("make")
    except:
        print("Could not compile RFM, please check manually that RFM is fine.")
    os.chdir(cwd)

    return print("Successfully compiled rfm.")


def construct_rfm_driver_table(inp, fldr, force=True, **kwargs):
    """Constructs RFM driver table.

    The table is a text file saved to *fldr*.

    Args:
        inp (dict): Dictionary containing RFM driver table inputs.
        fldr (str): Path to RFM-containing folder.
        force (bool): If True, overwrites the current driver table. If False, save a
            copy of the old table first.

    """

    # check inputs
    if not isinstance(inp, dict):
        raise TypeError("Inp must be a dictionary.")
    if not os.path.isdir(fldr):
        raise FileNotFoundError(f"{fldr} does not exist.")

    # check for primary sections
    first_five_secs = ["HDR", "FLG", "SPC", "GAS", "ATM"]
    sixth_sec = ["TAN", "GEO", "ELE", "SEC", "HGT", "LEN", "DIM"]
    for i in first_five_secs:
        if i not in inp.keys():
            raise ValueError(f"Section {i} missing from RFM inputs.")

    ss = [i for i in sixth_sec if i in inp.keys()]  # idnetifies the 6th section
    if len(ss) == 0:
        raise ValueError(f"One of these sections must be present: {str(sixth_sec)}.")
    elif len(ss) == 1:
        pass
    else:
        raise ValueError(
            f"Exactly one of these sections must be present: {str(sixth_sec)}."
        )

    # open the file in its final directory
    if force == True:
        print("Overwriting the current rmf.drv file.")
        f = open(f"{fldr}/rfm.drv", "w")
    elif force == False:
        print("The current rfm.drv will be moved to old_rfm.drv")
        os.rename(f"{fldr}/rfm.drv", f"{fldr}/old_rfm.drv")
        f = open("{fldr}/rfm.drv", "w")
    else:
        raise ValueError("'force' accepts only True or False. Default is True.")

    # format all input strings to be written
    for key in inp.keys():
        inp[key] = utils.line_break_str(txt=inp[key], chars=200, delim=" ", indent=2)

    # write the first five sections
    for fs in first_five_secs:
        f.write(f"*{fs}\n")
        f.write(f"  {inp[fs]}\n")

    # write sixth section
    f.write(f"*{ss[0]}\n")
    f.write(f"  {inp[ss[0]]}\n")

    # write the remaining sections
    for kkey in inp.keys():
        if kkey not in first_five_secs and kkey not in sixth_sec:
            f.write(f"*{kkey}\n")
            f.write(f"  {inp[kkey]}\n")
    f.write("*END")

    f.close()

    print("rfm.drv has been successfully created.")
    return


def construct_rfm_output_levels_file(levels, fldr, fname="alts.lev", force=True):
    """Construct the levels file for RFM - specifies output levels.

    Mind that the LEV flag for RFM must be enabled, or else this file will be ignored.
    Mind that the filename must be passed to RFM in the LEV section of the RFM driver
    table.
    A text file is created and saved with the required *fname* in *fldr*.

    Args:
        levels (array-like): Required output levels. Doesn't have to be sorted.
        fldr (str): output folder path
        fname (str): Required output filename, Default is "alts.lev".
        force (bool): If True, overwrite current levels file, else saves a copy first.

    Raises:
        TypeError: Raised when *levels* has incorrect type.
        ValueError: Raised when *levels* is an array, but has more than one dimension.
        TypeError: Raised when *fname* is not a string.
        FileNotFoundError: Raised if *fldr* does not exist.
        ValueError: Raised when *force* is not bool.

    """

    if not isinstance(levels, (list, np.ndarray)):
        raise TypeError("levels must be a list or a numpy.ndarray")
    if isinstance(levels, np.ndarray):
        if levels.ndim != 1:
            raise ValueError("If levels is a np.ndarray, it must have one dimension")
    if not os.path.isdir(fldr):
        raise FileNotFoundError(f"{fldr} does not exist.")
    if not isinstance(fname, str):
        raise TypeError("Parameter fname must be a string.")

    # open the file in its final directory
    if force == True:
        print(f"Overwriting the current {fname} file.")
        f = open(f"{fldr}/{fname}", "w")
    elif force == False:
        print(f"The current {fname} will be moved to old_{fname}.")
        os.rename(f"{fldr}/{fname}", f"{fldr}/old_{fname}")
        f = open(f"{fldr}/{fname}", "w")
    else:
        raise ValueError("'force' accepts only True or False. Default is True.")

    f.write("!Output levels\n")

    for i, ii in enumerate(levels):
        if (i % 5) == 0 and i != 0:
            f.write(f"\n{levels[i]:.4f}    ")
        else:
            f.write(f"{levels[i]:.4f}    ")

    f.write("\n*END")
    f.close()
    return


def construct_rfm_grid_file(wvnm, filename="grid.spc", rfm_fldr="./srfm/RFM"):
    """Creates a grid.spc file for rfm.

    The file is saved in *rfm_fldr* with the required filename.

    RFM is run in the "irregular grid mode" where the SPC section in the driver table
    contains the filename of this file.

    Args:
        wvnm (array-like): Wavenumbers, the calculation grid.
        filename (str): The output filename. Default is "grid.spc".
        rfm_fldr (str): Path to RFM folder. Default is "./srfm/RFM".

    Raises:
        TypeError: Raised when any of the input values have incorrect dtype.
    """

    if not isinstance(wvnm, (np.ndarray, list, pd.core.series.Series)):
        raise TypeError("wvnm must be a np.ndarray, pandas.core.series.Series or list")
    if not isinstance(filename, str):
        raise TypeError("filename must be a string.")

    path = f"{rfm_fldr}/rfm_files/{filename}"
    with open(path, "w") as f:
        f.write("!\n")
        f.write("!\n")
        f.write("!\n")
        f.write(f"{len(wvnm)} {wvnm[0]:.4f} 0 {wvnm[-1]:.4f}\n")
        for i in wvnm:
            f.write(f"{i:.4f} 0\n")
    f.close()
    print(f"{path} successfully created.")
    return f"./rfm_files/{filename}"


def read_atm_file(filename):  # read bits of internal profile output file prf.asc
    """This function reads the RFM .atm file.

    Args:
        filename: Path to RFM-format .atm filename.

    Returns:
        contents (dir): Dictionary with the file contents.

    Raises:
        TypeError: Raised when a section label is not a string. Indicates that the file
            was not corectly read.

    """

    f = open(filename, "r")
    f_lines = f.readlines()
    f.close()

    contents = {}  # dictionary to store that read values.

    # get section start lines
    sec_begin_lns = [f_lines.index(x) for x in f_lines if x.startswith("*")]
    for idx in range(len(sec_begin_lns) - 1):
        sec_lbl = f_lines[sec_begin_lns[idx]].strip(" \n*")
        if isinstance(sec_lbl, str) == False:
            raise TypeError("Section label must be a string.")

        sec_cont = [
            i.strip() for i in f_lines[sec_begin_lns[idx] + 1 : sec_begin_lns[idx + 1]]
        ]
        sec_cont = (" ").join(sec_cont)
        sec_cont = [float(ii) for ii in sec_cont.split()]
        contents[sec_lbl] = sec_cont

    return contents


def write_atm_file(data, filename, header=None):
    """Write out a file with atmospheric profile in RFM's .atm file format.

    Saves an .atm file in the with the specified name.

    Args:
        data (dictionary): Data to be written to the file.
        filename (str): Filename for the output file (has to be a path).
        header (str): Header for the file, optional, each line starts with "!". Default
            is None.

    """
    with open(filename, "w") as f:
        if header is not None:
            if header.startswith("!"):
                f.write(header)
            else:
                f.write("! ")
                f.write(header)
        f.write("\n")
        f.write("  ")
        f.write(str(len(data["HGT [km]"])))
        f.write(" ! No. of levels in profiles.")
        f.write("\n")
        f.write("*HGT [km]\n")
        for i, ii in enumerate(data["HGT [km]"]):
            if (i % 5) == 0 and i != 0:
                f.write(f'\n{data["HGT [km]"][i]:.4e}    ')
            else:
                f.write(f'{data["HGT [km]"][i]:.4e}    ')
        f.write("\n")
        keys = list(data.keys())
        keys.remove("HGT [km]")
        for key in keys:
            f.write(f"*{key}\n")
            for i, ii in enumerate(data[key]):
                if (i % 5) == 0 and i != 0:
                    f.write(f"\n{data[key][i]:.4e}    ")
                else:
                    f.write(f"{data[key][i]:.4e}    ")
            f.write("\n")
        f.write("*END")
    f.close()
    return
