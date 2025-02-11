import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import sys
import pandas as pd
from . import utilities as utils


def read_output(filename):  # read output file from RFM into a dictionary
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
    """This function reads the RFM output profile prf.asc."""
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


def get_rfm_optical_depths(fldr):
    """
    This function is designed to calculate optical depths for layers in
    the atmosphere from RFM with the LEV and OPT flags set on.
    The basic idea is to load optical depth spectra and subtract the adjacent ones
    to get layer optical depth.
    Returns a pandas dataframe.
    """

    # load relevant rfm output files, sort by modification time (i.e. altitude)
    fls = sorted(Path(fldr).iterdir(), key=os.path.getmtime)
    fls = [str(p.absolute()) for p in fls]
    fls = [
        i for i in fls if "down" in i and "opt" in i
    ]  # this is the list of OD output spectra

    prf = read_output_prf(f"{fldr}/prf.asc")  # read rfm output profile, type dict
    if len(fls) != len(prf["PRE [mb]"]):
        raise ValueError(
            "length of pressure profile and the number of"
            + " corresponding optical depth files don't match."
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

        d_OD = upper_level["SPC"] - lower_level["SPC"]  # layer optical depths
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

    pup = pd.DataFrame({"p_upper (mbar)" : u_p[::-1]}, dtype = "float")
    prf_df.loc[:,"p_upper (mbar)"] = pup
    
    plo = pd.DataFrame({"p_lower (mbar)" : l_p[::-1]}, dtype = "float")
    prf_df.loc[:,"p_lower (mbar)"] = plo
    
    pav = pd.DataFrame({"p_avg (mbar)" : avg_p[::-1]}, dtype = "float")
    prf_df.loc[:,"p_avg (mbar)"] = pav


    hup = pd.DataFrame({"h_upper (km)" : u_h[::-1]}, dtype = "float")
    prf_df.loc[:,"h_upper (km)"] = hup
    
    hlo = pd.DataFrame({"h_lower (km)" : l_h[::-1]}, dtype = "float")
    prf_df.loc[:,"h_lower (km)"] = hlo
    
    hav = pd.DataFrame({"h_avg (km)" : avg_h[::-1]}, dtype = "float")
    prf_df.loc[:,"h_avg (km)"] = hav


    tup = pd.DataFrame({"T_upper (K)" : u_t[::-1]}, dtype = "float")
    prf_df.loc[:,"T_upper (K)"] = tup
    
    tlo = pd.DataFrame({"T_lower (K)" : l_t[::-1]}, dtype = "float")
    prf_df.loc[:,"T_lower (K)"] = tlo
    
    tav = pd.DataFrame({"T_avg (K)" : avg_t[::-1]}, dtype = "float")
    prf_df.loc[:,"T_avg (K)"] = tav


    dod_col_names = [f"dOD_{val:.4f}" for val in wnos]
    iod_col_names = [f"iOD_{val:.4f}" for val in wnos]
    
    dod_df = pd.DataFrame(delta_OD, columns = dod_col_names)
    iod_df = pd.DataFrame(int_OD, columns = iod_col_names)

    df = pd.concat([prf_df,dod_df,iod_df], axis=1, join='outer')
    
#    for num, val in enumerate(wnos):
#        dod = pd.DataFrame({f"dOD_{val:.2f}" : delta_OD[:, num][::-1]}, 
#                            dtype = "float"
#                        )
#        df.loc[:, f"dOD_{val:.2f}"] = dod
#        
#        iod = pd.DataFrame({f"iOD_{val:.2f}" : int_OD[:, num][::-1]}, 
#                            dtype = "float"
#                        )
#        df.loc[:, f"iOD_{val:.2f}"] = iod

    return df


def compile_rfm(fldr):
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
    
def construct_rfm_driver_table(inp, fldr, force=True,**kwargs):
    
    #check inputs
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
    
    ss = [i for i in sixth_sec if i in inp.keys()] # idnetifies the 6th section
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
        f.write(f"{inp[fs]}\n")
    
    # write sixth section
    f.write(f"*{ss[0]}\n")
    f.write(f"{inp[ss[0]]}\n")
    
    # write the remaining sections
    for kkey in inp.keys():
        if kkey not in first_five_secs and kkey not in sixth_sec:
            f.write(f"*{kkey}\n")
            f.write(f"{inp[kkey]}\n")
    f.write("*END")
    
    f.close()
    
    print("rfm.drv has been successfully created.")
    return
        
        
