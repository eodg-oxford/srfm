"""This code defines one function to run srfm.

    Note that the SRFM can be run interactively (i.e. used as a package, use the
    inside of run_srfm() as and example), or with a driver table.
    
    This module is intended to run the srfm with preprocessed iasi spectra.

- Name: iasi_main
- Parent package: srfm
- Author: Antonin Knizek
- Contributors:
- Date: 7 Nov 2025
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from . import utilities
from . import forward_model
from . import rfm_functions
from . import layer
import os
import sys
import datetime
import time
import warnings
from multiprocessing import Process, Manager
from bisect import bisect
import pickle
from importlib.resources import files, as_file
from pathlib import Path
import importlib.util
import uuid
from .RFM import rfm_py
from . import rfm_helper
from mergedeep import merge
from netCDF4 import Dataset
import json
import copy

@utilities.show_runtime
def run_srfm(inp):
    """Main function that runs srfm.

    Optimized for IASI data.

    Args:
        inp (obj): Instance of inputs.Inputs.
    
    Returns:
        model_SRFM (obj): Instance of forward_model.SRFM.

    """
    ########################################################################################
    # Assign some variables:
    ########################################################################################
    with as_file(files("srfm") / "RFM") as path:
        rfm_fldr = os.fspath(path)

    ########################################################################################
    # set iasi grid and final grid to interpolate to
    ########################################################################################
    fin_grid = np.linspace(
        inp.values["fin_wvnmlo"],
        inp.values["fin_wvnmhi"],
        int(
            (inp.values["fin_wvnmhi"] - inp.values["fin_wvnmlo"])
            / inp.values["fin_res"]
            + 1
        ),
    )
    iasi_grid = np.linspace(645, 2760, int((2760 - 645) / 0.25 + 1))

    ########################################################################################
    # pick an iasi spectral file (processed)
    ########################################################################################
    iasi_spc_fldr = inp.values["iasi_spc_fldr"]
    iasi_fl = inp.values["iasi_fl"]  # A for ascending, D for descending
    date = datetime.datetime.strptime(
        iasi_fl[iasi_fl.find("_") + 1 : iasi_fl.rfind("_")], "%Y%m%d"
    )
    year_day = date.timetuple().tm_yday

    ########################################################################################
    # load corresponding spectral file
    ########################################################################################
    with open(os.path.join(iasi_spc_fldr, iasi_fl), "rb") as f:
        iasi_data = pickle.load(f)
    f.close()

    # pick pixel (manually check the file first and pick a pixel
    px = inp.values[
        "px"
    ]  # the only pixel with 0% cloud cover on this day in the ascending orbit for hlat ssacc
    keystr = iasi_fl[: iasi_fl.find(".")] + f"_px{px}"

    iasi_out_spc = iasi_data["spec_bbt"][px, :]

    iasi_out_spc = np.interp(fin_grid, iasi_grid, iasi_out_spc)

    iasi_out_zen_deg = iasi_data["zen"][px]  # zenith angle [degrees]

    iasi_out_zen_rad = np.deg2rad(iasi_out_zen_deg)  # avg zenith angle [rad]
    iasi_zen_cos = np.cos(iasi_out_zen_rad).item()  # cosine of the zenith angle
    iasi_zen_sec = 1 / iasi_zen_cos  # secant of zenith angle
    # secant is for RFM, currently unused

    ########################################################################################
    # prepare ECMWF data (now from the spectral file)
    ########################################################################################
    z = iasi_data["ecmwf_z"][px][::-1]  # altitude grid

    slc_idx = next((index for index, value in enumerate(list(z)) if value > 0), None)

    # crop leading negative values
    z = z[slc_idx:]

    # load other ecmwf data
    T = iasi_data["ecmwf_T"][px][::-1][slc_idx:]
    p = iasi_data["ecmwf_p"][px][::-1][slc_idx:]
    o3 = iasi_data["ecmwf_O3"][px][::-1][slc_idx:]
    co = iasi_data["ecmwf_CO"][px][::-1][slc_idx:]
    n2o = iasi_data["ecmwf_N2O"][px][::-1][slc_idx:]
    co2 = iasi_data["ecmwf_CO2"][px][::-1][slc_idx:]
    ch4 = iasi_data["ecmwf_CH4"][px][::-1][slc_idx:]
    h2o = iasi_data["ecmwf_H2O"][px][::-1][slc_idx:]

    # read current rfm .atm file used to specify the atmosphere
    rfm_prf = rfm_functions.read_atm_file(inp.values["driver_inputs"]["atmosphere"][1])

    """The existing rfm data is bisected such that the levels which overlap with the
    ecmwf data levels are selected. On this part of the profile, the ecmwf data are used
    (interpolated to the .atm grid). The ecmwf data ends at ~75-80 km, so rfm day.atm 
    data are used above that point.    
    """
    b_idx = bisect(
        rfm_prf["HGT [km]"], z[-1]
    )  # bisect existing profile at end of ecmwf data

    z_lo = rfm_prf["HGT [km]"][:b_idx]

    # keep profiles where there's no ecmwf data (higher up in the atmosphere)
    new_T_hi = rfm_prf["TEM [K]"][b_idx:]
    new_p_hi = rfm_prf["PRE [mb]"][b_idx:]
    new_o3_hi = rfm_prf["O3 [ppmv]"][b_idx:]
    new_co_hi = rfm_prf["CO [ppmv]"][b_idx:]
    new_n2o_hi = rfm_prf["N2O [ppmv]"][b_idx:]
    new_co2_hi = rfm_prf["CO2 [ppmv]"][b_idx:]
    new_ch4_hi = rfm_prf["CH4 [ppmv]"][b_idx:]
    new_h2o_hi = rfm_prf["H2O [ppmv]"][b_idx:]

    # interpolate ecmwf data to the altitude grid from the .atm file, overlapping part
    new_T_lo = np.interp(z_lo, z, T)
    new_p_lo = np.interp(z_lo, z, p)
    new_o3_lo = np.interp(z_lo, z, o3)
    new_co_lo = np.interp(z_lo, z, co)
    new_n2o_lo = np.interp(z_lo, z, n2o)
    new_co2_lo = np.interp(z_lo, z, co2)
    new_ch4_lo = np.interp(z_lo, z, ch4)
    new_h2o_lo = np.interp(z_lo, z, h2o)

    # concatenate the two parts of the array
    new_T = np.concatenate((new_T_lo, np.array(new_T_hi)))
    new_p = np.concatenate((new_p_lo, np.array(new_p_hi)))
    new_o3 = np.concatenate((new_o3_lo, np.array(new_o3_hi)))
    new_co = np.concatenate((new_co_lo, np.array(new_co_hi)))
    new_n2o = np.concatenate((new_n2o_lo, np.array(new_n2o_hi)))
    new_co2 = np.concatenate((new_co2_lo, np.array(new_co2_hi)))
    new_ch4 = np.concatenate((new_ch4_lo, np.array(new_ch4_hi)))
    new_h2o = np.concatenate((new_h2o_lo, np.array(new_h2o_hi)))

    if inp.values["plot_profiles"] == True:
        # plot old and new profiles and save figure
        fig, axs = plt.subplots(2, 4, figsize=(11.7, 8.3), sharey=True)
        plt.rcParams.update({"font.size": 14})
        fig.supylabel("Altitude (km)")

        axs[0, 0].plot(rfm_prf["TEM [K]"], rfm_prf["HGT [km]"], label="rfm")
        axs[0, 0].plot(new_T, rfm_prf["HGT [km]"], label="ecmwf+rfm")
        axs[0, 0].set_xlabel("T (K)")
        axs[0, 0].legend()

        axs[0, 1].plot(rfm_prf["PRE [mb]"], rfm_prf["HGT [km]"], label="rfm")
        axs[0, 1].plot(new_p, rfm_prf["HGT [km]"], label="ecmwf+rfm")
        axs[0, 1].set_xlabel("p (mb)")
        axs[0, 1].legend()

        axs[0, 2].plot(rfm_prf["O3 [ppmv]"], rfm_prf["HGT [km]"], label="rfm")
        axs[0, 2].plot(new_o3, rfm_prf["HGT [km]"], label="ecmwf+rfm")
        axs[0, 2].set_xlabel(r"$x_{O_{3}}$ (ppmv)")
        axs[0, 2].legend()

        axs[0, 3].plot(rfm_prf["CO [ppmv]"], rfm_prf["HGT [km]"], label="rfm")
        axs[0, 3].plot(new_co, rfm_prf["HGT [km]"], label="ecmwf+rfm")
        axs[0, 3].set_xlabel(r"$x_{CO}$ (ppmv)")
        axs[0, 3].legend()

        axs[1, 0].plot(rfm_prf["N2O [ppmv]"], rfm_prf["HGT [km]"], label="rfm")
        axs[1, 0].plot(new_n2o, rfm_prf["HGT [km]"], label="ecmwf+rfm")
        axs[1, 0].set_xlabel(r"$x_{N_{2}O}$ (ppmv)")
        axs[1, 0].legend()

        axs[1, 1].plot(rfm_prf["CO2 [ppmv]"], rfm_prf["HGT [km]"], label="rfm")
        axs[1, 1].plot(new_co2, rfm_prf["HGT [km]"], label="ecmwf+rfm")
        axs[1, 1].set_xlabel(r"$x_{CO_{2}}$ (ppmv)")
        axs[1, 1].legend()

        axs[1, 2].plot(rfm_prf["CH4 [ppmv]"], rfm_prf["HGT [km]"], label="rfm")
        axs[1, 2].plot(new_ch4, rfm_prf["HGT [km]"], label="ecmwf+rfm")
        axs[1, 2].set_xlabel(r"$x_{CH_{4}}$ (ppmv)")
        axs[1, 2].legend()

        axs[1, 3].plot(rfm_prf["H2O [ppmv]"], rfm_prf["HGT [km]"], label="rfm")
        axs[1, 3].plot(new_h2o, rfm_prf["HGT [km]"], label="ecmwf+rfm")
        axs[1, 3].set_xlabel(r"$x_{H_{2}O}$ (ppmv)")
        axs[1, 3].set_xscale = "log"
        axs[1, 3].legend()

        # plt.suptitle(f"Comparison between original rfm and interpolated ecmwf+rfm profiles\n{keystr[1:keystr.rfind('.')]}")
        plt.suptitle(
            f"Comparison between original rfm and interpolated ecmwf+rfm profiles\n{keystr}"
        )

        plt.tight_layout()
        plt.savefig(f"{inp.values['results_fldr']}/{keystr}_prf.png")
        plt.close()

    ########################################################################################
    # work out where to insert a liquid water cloud
    ########################################################################################

    # (exclude region above 2200 cm-1 where reflected solar radiation starts to matter)
    lim = (
        fin_grid[-1] if fin_grid[-1] < 2200 else 2200
    )  # determine where to cut off the
    # spectrum for looking for the maximum observed bbt

    # index of maximum observed bbt in the truncated spectrum, which is in theory the
    # observed blackbody temperature + instrument noise
    max_obs_bbt_idx = iasi_out_spc[0 : np.absolute(iasi_grid - lim).argmin()].argmax()

    obs_bbt = iasi_out_spc[max_obs_bbt_idx]  # maximum observed temperature

    obs_bbt_wvnm = fin_grid[max_obs_bbt_idx]  # wavenumber corresponding to the
    # maximum in the spectrum

    # check for ambiguity (should not happen, this is just check, not a fix)
    if not isinstance(obs_bbt_wvnm, float):
        raise RuntimeError(
            """Ambiguity encountered. Multiple wavenumbers matched to the
        maximum brightness temperature value. Can't assign noise equivalent delta
        temperature."""
        )

    # concept of the following few lines:
    # 1. Assumption 1: The main cause of noise in the spectrum is instrument noise.
    # 2. iasi_nedt.txt is a file that contains iasi NEDT @280K, see the file for details
    # 3. I interpolate the nedt to obs_bbt_wvnm
    # 4. Since this is noise, I take the value*1/2 and subtract from obs_bbt
    # 5. In this way I have hopefully removed the effect of noise from my reading.
    # 6. Assumption 2: obs_bbt does not coincide with a strong absorption line (ie is
    #    truly from a black/grey body.
    # 7. Assumption 3: NEDT is taken just as the instrument noise at 280K, no
    #    interpolation to other temperatures (generally in this case, the difference
    #    should be <20 K anyway.
    # 8. Assumption 3: Scene specific NEDT is not used (see Vincent, A, Retrieval of
    #    trace gases using IASI, thesis, Univ Oxford (EODG), 2016
    # 9. Assumption 4: The noise is random and taking a mean, or halving it, makes sense.

    # open iasi noise equivalent delta temperature file
    nedt = np.loadtxt(inp.values["nedt"], skiprows=3)  # [[wvnm,nedt]]

    # interpolate to obs_bbt_wvnm
    obs_bbt_nedt = np.interp(obs_bbt_wvnm, nedt[:, 0], nedt[:, 1])

    # subtract half of the noise value from the temperature value
    obs_bbt = obs_bbt - obs_bbt_nedt / 2

    # find tropopause (min in TEM profile between 0 and 30 km, 30 is arbitrary)
    tropopause_idx = rfm_prf["TEM [K]"].index(min(rfm_prf["TEM [K]"][0:30]))

    # get the altitude from a the ecmwf profile that corresponds to this temperature
    # the profile is linearly interpolated.
    corr_alt = np.interp(
        obs_bbt,
        rfm_prf["TEM [K]"][1:tropopause_idx][::-1],  # the slice [1:10] is to avoid
        # the fact that the temperature profile isn't monotonic with
        # altitude. Using the first ten values should both be
        # sufficient to capture clouds and small enough to avoid
        # encountering the cold trap temperature.
        # The profile is inverted for interpolation purposes, because
        # np.interp requires strictly increasing values in abscissa.
        rfm_prf["HGT [km]"][1:tropopause_idx][::-1],
    )
    corr_alt = round(corr_alt, 3)
    corr_alt = corr_alt.item()

    ########################################################################################
    # update rfm profiles with data from previous two sections
    ########################################################################################

    # update first temperature in the profile according to the maximum observed temperature
    # minus nedt. The logic being that if there's no cloud, this correspond to surface
    # temperature. This is necessary since the spectrum is sensitive to surface temperature
    # and ecmwf data may be inaccurate in this. If there is a cloud, then the surface
    # temperature does not matter anyway.
    # new_T[0] = iasi_out_spc[max_obs_bbt_idx]

    # update rfm_prf
    rfm_prf["TEM [K]"] = new_T
    rfm_prf["PRE [mb]"] = new_p
    rfm_prf["O3 [ppmv]"] = new_o3
    rfm_prf["CO [ppmv]"] = new_co
    rfm_prf["N2O [ppmv]"] = new_n2o
    rfm_prf["CO2 [ppmv]"] = new_co2
    rfm_prf["CH4 [ppmv]"] = new_ch4
    rfm_prf["H2O [ppmv]"] = new_h2o

    # set profiles of species to retrieve
    if "g_rtv" in inp.values.keys():
        for key in inp.values["g_rtv"]:
            rfm_prf[key] = inp.values["g_rtv"][key]

    # save output profile
    rfm_functions.write_atm_file(
        data=rfm_prf,
        filename=f"{inp.values['results_fldr']}/{keystr}.atm",
        header=f"""!Combination of (ecmwf) and MID-LATITUDE DAY PROFILES FOR MIPAS
        !by Anu Dudhia.\n""",
    )

    ########################################################################################
    # specify spectral calculation grid
    ########################################################################################
    spec_res = inp.values["spc_res"]  # model spectral resolution,[spec_units]
    low_spc = inp.values["spc_wvnmlo"]  # model start wavenumber (lower), [spec_units]
    upp_spc = inp.values["spc_wvnmhi"]  # model end wavenumber (upper), [spec_units]
    spec_units = inp.values["spc_units"]  # accepted values "cm-1", "um", "nm"

    RFM_wvnm, wvls = utilities.calc_grids(low_spc, upp_spc, spec_res, spec_units)

    rfm_grid_fname = rfm_functions.construct_rfm_grid_file(
        RFM_wvnm, filename="grid.spc", rfm_fldr=rfm_fldr
    )
    ########################################################################################
    # define an atmospheric scattering layers
    ########################################################################################
    scat_lyrs = (
        {}
    )  # dictionary of scattering layers, key - layer name, value - Layer() object

    # define Layer properties
    scat_lyrs_inputs = {}

    if "scat_lyrs_inputs" in inp.values.keys():

        # calculate MieLayer optical properties
        for lyr in inp.values["scat_lyrs_inputs"].keys():
            scat_lyrs_inputs[lyr] = (
                inp.values["scat_lyrs_inputs"][lyr]
                | inp.values["scat_lyrs_inputs"][lyr]
            )
            scat_lyrs[lyr] = layer.MieLayer()
            scat_lyrs[lyr].set_input_from_dict(
                scat_lyrs_inputs[lyr]
            )  # sets input for scattering layer
            scat_lyrs[
                lyr
            ].calculate_op()  # calculates layer optical properties, may run in parallel

    ########################################################################################
    # prepare atmospheric layer structure
    ########################################################################################

    # define some requested output levels
    levels = [
        0.0,
        1.0,
        2.0,
        3.0,
        4.0,
        5.0,
        6.0,
        7.0,
        8.0,
        9.0,
        10.0,
        11.0,
        12.0,
        13.0,
        14.0,
        15.0,
        16.0,
        17.0,
        18.0,
        19.0,
        20.0,
        21.0,
        22.0,
        23.0,
        24.0,
        25.0,
        27.5,
        30.0,
        32.5,
        35.0,
        37.5,
        40.0,
        42.5,
        45.0,
        47.5,
        50.0,
        55.0,
        57.5,
        60.0,
        62.5,
        65.0,
        67.5,
        70.0,
        75.0,
        80.0,
        82.5,
        85.0,
        87.5,
        90.0,
        95.0,
        100.0,
        101.0,
        102.0,
        103.0,
        104.0,
        105.0,
        106.0,
        107.0,
        108.0,
        109.0,
        110.0,
        110.5,
        111.0,
        111.5,
        112.0,
        112.5,
        113.0,
        113.5,
        114.0,
        114.5,
        115.0,
        115.5,
        116.0,
        116.5,
        117.0,
        117.5,
        118.0,
        118.5,
        119.0,
        119.5,
        120.0,
    ]

    # define tracker array for atmospheric structure
    track_lev = [None for i in levels]

    # add upper and lower particle layer boundaries, delete any levels "within" the layer
    for lyr in scat_lyrs.keys():
        levels, track_lev = utilities.add_lyr_from_Layer(
            lev=levels, track_lev=track_lev, new_lyr=scat_lyrs[lyr]
        )

    # convert the tracking levels array to a tracking layers array
    track_lyr = utilities.track_lev_to_track_lyr(track_lev)
    track_lyr = track_lyr[::-1]

    ########################################################################################
    # prepare and call RFM
    ########################################################################################
    # RFM global config
    rfm_config = inp.values["rfm_config"]

    # RFM driver table
    driver_inputs = inp.values["driver_inputs"]
    driver_inputs["tangent"] = (str(iasi_zen_sec),)
    driver_inputs["lev"] = tuple(str(val) for val in levels)
    driver_inputs["atmosphere"] = list(driver_inputs["atmosphere"])
    driver_inputs["atmosphere"][1] = f"{inp.values['results_fldr']}/{keystr}.atm"
    driver_inputs["atmosphere"] = tuple(driver_inputs["atmosphere"])

    # initialize RFM model class
    model_RFM = forward_model.RFM()

    # print current status:
    print(model_RFM.status)

    # run rfm
    rfm_run_result = rfm_helper.rfm_main(
        configuration=rfm_config,
        driver_inputs=driver_inputs,
        levels=levels,
        rfm_out_fldr=inp.values["results_fldr"],
    )
    # Store the full RunResult for debugging/metadata while keeping the legacy dataframe API.
    model_RFM.rfm_run_result = rfm_run_result
    if rfm_run_result.output_df is None:
        raise RuntimeError("RFM capture did not return an optical-depth dataframe.")
    model_RFM.rfm_output = rfm_run_result.output_df.copy()

    # print current status:
    print(model_RFM.status)

    # rebuild the spectral grid directly from the captured columns
    cols = [i for i in model_RFM.rfm_output.columns if i.startswith("dOD_")]
    if not cols:
        raise RuntimeError(
            "RFM output does not include any differential optical-depth columns."
        )
    try:
        RFM_wvnm = np.array([float(col.split("_", 1)[1]) for col in cols], dtype=float)
    except ValueError as exc:
        raise RuntimeError(
            "Failed to parse wavenumbers from RFM output columns."
        ) from exc
    wvls = (1.0 / RFM_wvnm) * 1e4

    # add output from optical properties calculation, TODO MOVE UP?
    for lyr in scat_lyrs.keys():
        scat_lyrs[lyr].add_op_calc_output()

        # interpolate layer optical properties
        scat_lyrs[lyr].regrid(wvls, track_diff=False)
        scat_lyrs[lyr].calc_tau()

    ########################################################################################
    # prepare DISORT common variables
    ########################################################################################

    # initialize DISORT model class
    model_DISORT = forward_model.DISORT()

    # check if number of columns and wavelengths match
    if len(cols) != len(wvls):
        raise ValueError(
            f"Number of RFM and scattering wavelengths don't match "
            f"({len(cols)} spectral columns vs {len(wvls)} scattering wavelengths)."
        )

    # set disort_input parameters common to all loop iterations
    # these need to be set first:
    nmom = inp.values["nmom"]
    for lyr in scat_lyrs.keys():
        if (scat_lyrs[lyr].legendre_coefficient.shape[1] - 1) > nmom:
            nmom = scat_lyrs[lyr].legendre_coefficient.shape[1] - 1

    model_DISORT.set_maxcmu(inp.values["maxcmu"])

    model_DISORT.set_maxmom(nmom)
    if nmom < model_DISORT.disort_input["maxcmu"]:
        model_DISORT.set_maxmom(model_DISORT.disort_input["maxcmu"])

    model_DISORT.set_maxumu(inp.values["maxumu"])
    model_DISORT.set_maxphi(inp.values["maxphi"])
    model_DISORT.set_maxulv(inp.values["maxulv"])

    # now the rest
    model_DISORT.set_usrang(inp.values["usrang"])
    model_DISORT.set_usrtau(inp.values["usrtau"])
    model_DISORT.set_ibcnd(inp.values["ibcnd"])
    model_DISORT.set_onlyfl(inp.values["onlyfl"])
    # model_DISORT.set_prnt([True, True, True, False, False])
    model_DISORT.set_prnt(inp.values["prnt"])
    model_DISORT.set_plank(inp.values["planck"])
    model_DISORT.set_lamber(inp.values["lamber"])
    model_DISORT.set_deltamplus(inp.values["deltamplus"])
    model_DISORT.set_do_pseudo_sphere(inp.values["do_pseudo_sphere"])
    model_DISORT.set_utau(inp.values["utau"])

    model_DISORT.set_fisot(inp.values["fisot"])
    model_DISORT.set_albedo(inp.values["albedo"])

    model_DISORT.set_temis(inp.values["temis"])
    model_DISORT.set_earth_radius(6371)
    model_DISORT.set_rhoq(
        np.zeros(
            shape=(
                int(model_DISORT.disort_input["maxcmu"] / 2),
                int(model_DISORT.disort_input["maxcmu"] / 2 + 1),
                int(model_DISORT.disort_input["maxcmu"]),
            )
        )
    )
    model_DISORT.set_rhou(
        np.zeros(
            shape=(
                model_DISORT.disort_input["maxumu"],
                int(model_DISORT.disort_input["maxcmu"] / 2 + 1),
                model_DISORT.disort_input["maxcmu"],
            )
        )
    )
    model_DISORT.set_rho_accurate(
        np.zeros(
            shape=(
                model_DISORT.disort_input["maxumu"],
                model_DISORT.disort_input["maxphi"],
            )
        )
    )
    model_DISORT.set_bemst(
        np.zeros(shape=(int(model_DISORT.disort_input["maxcmu"] / 2)))
    )
    model_DISORT.set_emust(np.zeros(shape=(model_DISORT.disort_input["maxumu"])))
    model_DISORT.set_accur(0)
    ########################################################################################
    # prepare solar spectrum
    ########################################################################################

    # load solar spectrum from file
    solar_spc, solar_spc_wvnm = utilities.load_solar_spectrum_Gueymard20018()
    solar_spc = np.interp(
        RFM_wvnm, solar_spc_wvnm[::-1], solar_spc[::-1]
    )  # interpolate to RFM_wvnm (the calculation grid)

    # scale with year day (different Sun-Earth distance throughout the year
    # the original specturm is for 1 AU
    solar_spc = utilities.scale_solar_spectrum(solar_spc, year_day)

    # get incoming solar beam polar angle for DISORT
    solar_zen_deg = iasi_data["sza"][px]  # solar zenith angle [degrees]
    solar_zen_rad = np.deg2rad(solar_zen_deg)  # solar zenith angle [rad]
    solar_zen_cos = np.cos(
        solar_zen_rad
    ).item()  # cosine of the solar zenith angle, UMU0

    model_DISORT.set_umu0(solar_zen_cos)
    model_DISORT.set_phi0(iasi_data["saa"][px])
    model_DISORT.set_umu([iasi_zen_cos])
    model_DISORT.set_phi([iasi_data["azi"][px]])

    # initialize disort input arrays for output variables from a single run
    model_DISORT.initialize_disort_output_arrays()

    ########################################################################################
    # set DISORT DISOBRDF variables
    ########################################################################################

    # TBD

    ########################################################################################
    # prepare SRFM common variables
    ########################################################################################
    model_SRFM = forward_model.SRFM()
    model_SRFM.set_wvnm(RFM_wvnm)
    model_SRFM.set_wvls(wvls)
    model_SRFM.initialize_srfm_output_arrays_from_disort(model_DISORT)

    # track progress
    pct = [
        1,
        2,
        5,
        10,
        20,
        30,
        40,
        50,
        60,
        70,
        80,
        90,
        100,
    ]  # percent completed to report
    pct_val = [RFM_wvnm[int(i * len(RFM_wvnm) / 100 - 1)] for i in pct]

    ########################################################################################
    # set dynamic variables and run DISORT
    ########################################################################################
    # loop over columns, dynamically set disort input variables in each loop
    for wvl_idx, (wvnm, wvl, col) in enumerate(zip(RFM_wvnm, wvls, cols)):

        # track progress
        if wvnm in pct_val:
            val_idx = pct_val.index(wvnm)
            print(f"Running main DISORT loop. {pct[val_idx]}% done...")

        #        model_DISORT.set_header(f"Now starting calculation for {col} cm-1.")
        model_DISORT.set_header(
            "NO HEADER"
        )  # the string "NO HEADER" will cause no printout
        model_DISORT.set_wvnm(wvnm)
        model_DISORT.set_wvl(wvl)

        # get layer optical depths from gas absorption
        tau_g = model_RFM.rfm_output[col].to_numpy()

        # layer optical depths from Rayleigh scattering
        #    tau_R = np.zeros(shape=(tau_g.shape))
        tau_R = utilities.calc_Rayleigh_opt_depths(
            ps=model_RFM.rfm_output["p_lower (mbar)"].iloc[-1],
            pu=model_RFM.rfm_output["p_upper (mbar)"],
            pl=model_RFM.rfm_output["p_lower (mbar)"],
            l=wvnm,
        )

        # particle layer optical depths (from particle scattering)
        tau_p = np.zeros(shape=(len(tau_g)))
        for lyr in scat_lyrs.keys():
            tau_p[track_lyr.index(lyr)] = scat_lyrs[lyr].tau[wvl_idx]

        # particle layer single scatter albedo
        w_p = np.zeros(shape=(len(tau_g)))
        for lyr in scat_lyrs.keys():
            w_p[track_lyr.index(lyr)] = scat_lyrs[lyr].ssalb[wvl_idx]

        dtauc_tot = utilities.calc_tot_dtauc(tau_g=tau_g, tau_R=tau_R, tau_p=tau_p)

        # truncate optical depths
        threshold_od = 1e-8  # threshold at which to truncate optical depths
        idx = next(
            (
                index
                for index, value in enumerate(list(dtauc_tot))
                if value > threshold_od
            ),
            None,
        )
        dtauc_tot = dtauc_tot[idx:]
        tau_g = tau_g[idx:]
        tau_R = tau_R[idx:]
        tau_p = tau_p[idx:]
        w_p = w_p[idx:]
        track_lyr_local = track_lyr[idx:]

        # set layer optical depths
        model_DISORT.set_dtauc_manually(dtauc=dtauc_tot)
        model_DISORT.disort_input["dtauc"][model_DISORT.disort_input["dtauc"] < 0] = 0

        # set maxcly based on current length of tau_g
        model_DISORT.set_maxcly(len(dtauc_tot))

        ## set some maxcly-dependent DISORT input variables

        # set and truncate temper
        model_DISORT.set_temper_from_rfm(model_RFM)
        model_DISORT.disort_input["temper"] = model_DISORT.disort_input["temper"][idx:]

        model_DISORT.set_btemp(model_DISORT.disort_input["temper"][-1])
        model_DISORT.set_ttemp(model_DISORT.disort_input["temper"][0])
        model_DISORT.set_h_lyr(
            np.zeros(shape=(model_DISORT.disort_input["maxcly"] + 1))
        )

        # set single scatter albedo
        model_DISORT.set_ssalb(tau_g=tau_g, tau_R=tau_R, tau_p=tau_p, w_p=w_p)

        # calculate phase function moments for Rayleigh and particle scattering from DISORT
        pmom_R = model_DISORT.calc_pmom(iphas=2, prec=inp.values["disort_precision"])

        # set phase function moments for particle scattering from Mie code
        pmom_p = np.zeros(
            (
                model_DISORT.disort_input["maxmom"] + 1,
                model_DISORT.disort_input["maxcly"],
            )
        )
        for lyr in scat_lyrs.keys():
            if scat_lyrs[lyr].tau[wvl_idx] < threshold_od:
                print(
                    f"""Scattering layer optical depth was < {threshold_od} and was
                truncated from the optical depths profile. No particle scattering at this
                 wavelength."""
                )
            else:
                pmom_p[
                    : len(scat_lyrs[lyr].legendre_coefficient[wvl_idx, :]),
                    track_lyr_local.index(lyr),
                ] = scat_lyrs[lyr].legendre_coefficient[wvl_idx, :]
                Legendre_precision = (
                    1 / pmom_p[0, track_lyr_local.index(lyr)]
                )  # unused, can be used to track expansion precision
                if abs(Legendre_precision - 1) > 1e-5:
                    raise RuntimeError("Something is wrong with the phase function")
                pmom_p[0, track_lyr_local.index(lyr)] = 1.0

        # calculate the weighted sum of phase function moments
        model_DISORT.set_pmom(
            pmom_R=pmom_R, tau_R=tau_R, w_p=w_p, tau_p=tau_p, pmom_p=pmom_p
        )
        #    print(model_DISORT.disort_input["pmom"]).shape

        # set wavenumber range for DISORT (for Planck function)
        model_DISORT.set_wvnm_range(wvnm - 0.5, wvnm + 0.5)

        # set incoming beam of (solar) radiation
        model_DISORT.set_fbeam(solar_spc[wvl_idx])
        #    model_DISORT.set_fbeam(0.1)

        # run disort input tests
        model_DISORT.test_disort_input_format()
        model_DISORT.test_disort_input_integrity()

        # call DISOBRDF
        #    model_DISORT.run_disobrdf(prec=inp.values["disort_precision"],
        #                              debug=False,
        #                              brdf_type=2, # Cox-Munk
        #                              brdf_arg=[1,1.34,False,0], # wind speed, water refractive index, do_shadow
        #                              nmug=200) # number of quadrature angles
        # run disort
        model_DISORT.run_disort(prec=inp.values["disort_precision"])

        # store result in SRFM()
        model_SRFM.store_disort_result(model_DISORT, wvl_idx)

        # print current status:m
    #    print(model_DISORT.status)
    print("Main DISORT loop finished.")

    # convolve final radiance spectrum with iasi instrument line shape
    model_SRFM.convolve_with_iasi(inp.values["ils"])

    # interpolate resulting bbt and radiances to final grid
    model_SRFM.interp(fin_grid)

    # calculate brightness temperature for the final spectrum
    model_SRFM.calc_bbt()

    ########################################################################################
    # (optional) save spectrum to file
    ########################################################################################
    if inp.values["out_mode"] == "txt":
        if inp.values["bbt"] == True:
            if isinstance(inp.values["bbt_out_fname"], str):
                np.savetxt(
                    f"{inp.values['results_fldr']}/{inp.values['bbt_out_fname']}.txt",
                    np.column_stack((model_SRFM.wvnm, model_SRFM.bbt[:, 0, 0, 0])),
                )       
            else:         
                np.savetxt(
                    f"{inp.values['results_fldr']}/bbt.txt",
                    np.column_stack((model_SRFM.wvnm, model_SRFM.bbt[:, 0, 0, 0])),
                )
        if inp.values["rad"] == True:
            if isinstance(inp.values["rad_out_fname"], str):
                np.savetxt(
                    f"{inp.values['results_fldr']}/{inp.values['rad_out_fname']}.txt",
                    np.column_stack((model_SRFM.wvnm, model_SRFM.uu[:, 0, 0, 0])),
                )       
            else:         
                np.savetxt(
                    f"{inp.values['results_fldr']}/rad.txt",
                    np.column_stack((model_SRFM.wvnm, model_SRFM.uu[:, 0, 0, 0])),
                )
        else:
            warnings.warn("driver table doesn't specify output bbt or rad.")
            pass
        
    elif inp.values["out_mode"] == "netcdf":
        if inp.values["bbt"] == True:
            if isinstance(inp.values["bbt_out_fname"], str):
                out_nm = f"{inp.values['results_fldr']}/{inp.values['bbt_out_fname']}.nc" 
            else:
                out_nm = f"{inp.values['results_fldr']}/bbt.nc"
            
            with Dataset(out_nm, 'w', format='NETCDF4') as nc_file:
                nc_file.description = f"SRFM output."
                nc_file.history = f"Created {datetime.datetime.now().strftime('%Y-%m-%d')}"
                nc_file.createDimension("wavenumber", fin_grid.shape[0]) # determines the output spectrum shape
                bbt = nc_file.createVariable("bbt", "f8", ("wavenumber",), zlib=True, complevel=4)
                
                bbt.units = "K"
                bbt.long_name = "Brightness temperature"
                bbt[:] = model_SRFM.bbt[:, 0, 0, 0]
                
                # store inputs as well
                metadata = copy.deepcopy(inp.values)
                metadata["driver_inputs"]["spectral"] = str(metadata["driver_inputs"]["spectral"])
                _inp = json.dumps(metadata)
                bbt.srfm_params = _inp
        
        if inp.values["rad"] == True:
            if isinstance(inp.values["rad_out_fname"], str):
                out_nm = f"{inp.values['results_fldr']}/{inp.values['rad_out_fname']}.nc" 
            else:
                out_nm = f"{inp.values['results_fldr']}/rad.nc"
            
            with Dataset(out_nm, 'w', format='NETCDF4') as nc_file:
                nc_file.description = f"SRFM output."
                nc_file.history = f"Created {datetime.datetime.now().strftime('%Y-%m-%d')}"
                nc_file.createDimension("wavenumber", fin_grid.shape[0]) # determines the output spectrum shape
                rad = nc_file.createVariable("rad", "f8", ("wavenumber",), zlib=True, complevel=4)
                
                rad.units = "W m-2 sr-1 cm"
                rad.long_name = "Radiance"
                rad[:] = model_SRFM.uu[:, 0, 0, 0]
                
                # store inputs as well
                metadata = copy.deepcopy(inp.values)
                metadata["driver_inputs"]["spectral"] = str(metadata["driver_inputs"]["spectral"])
                _inp = json.dumps(metadata)
                bbt.srfm_params = _inp
            
    elif inp.values["out_mode"] == None:
        pass
        

    if inp.values["base_plots"] == True:
        ########################################################################################
        # (optional) create base plots
        ########################################################################################
        plot = False  # create a plot?
        plot_srfm = True  # plot disort output?
        plot_rfm = False  # plot rfm output?
        plot_residual = False  # plot difference between rfm and diosrt?
        y_type = (
            "bbt"  # plot in radiances (W m-2 sr-1 cm) or brightness temperatures [K]
        )
        x_type = "cm-1"  # plot vs. wavenumbers [cm-1] or wavelengths [um] or [nm]

        if plot == True:
            plt.figure(figsize=(11.7, 8.4))
            plt.rcParams.update({"font.size": 12})
            plt.cla()

            # determine x label:
            if x_type == "cm-1":
                x_lbl = r"Wavenumbers (cm$^{-1}$)"
            elif x_type == "um":
                x_lbl = r"Wavelength ($\mu$m)"
            elif x_type == "nm":
                x_lbl = "Wavelength (nm)"

            # determine y label:
            if y_type == "bbt":
                y_lbl = "Brightness temperature (K)"
            elif y_type == "rad":
                y_lbl = r"Radiance (W m$^{-2}$ sr$^{-1}$ cm)"

            # plot RFM
            if plot_rfm == True:
                # determine x:
                if x_type == "cm-1":
                    x = RFM_wvnm
                elif x_type == "um":
                    x = wvls
                elif x_type == "nm":
                    x = wvls * 1e3

                filename = f"{rfm_fldr}/rad_001000.asc"
                data = rfm_functions.read_output(filename)
                RFM_rad = data["SPC"] * 1e-5

                if y_type == "bbt":
                    RFM_bbt = utilities.convert_spectral_radiance_to_bbt(
                        RFM_rad, RFM_wvnm
                    )
                    plt.plot(
                        x, RFM_bbt, label=f"no scattering", alpha=1, c="tab:orange"
                    )
                elif y_type == "rad":
                    plt.plot(
                        x, RFM_rad, label=f"no scattering", alpha=0.9, c="tab:orange"
                    )

            # plot SRFM
            if plot_srfm == True:
                # determine x:
                if x_type == "cm-1":
                    x = model_SRFM.wvnm
                elif x_type == "um":
                    x = model_SRFM.wvls
                elif x_type == "nm":
                    x = model_SRFM.wvls * 1e3

                if y_type == "bbt":
                    y = model_SRFM.bbt[:, 0, 0, 0]
                elif y_type == "rad":
                    y = model_SRFM.uu[:, 0, 0, 0]

                plt.plot(x, y, label=f"scattering", c="tab:blue")

            # plot residual
            if plot_residual == True:
                if plot_srfm != True and plot_rfm != True:
                    print(f"Can't plot residual without data.")
                else:
                    plt.ion()
                    if y_type == "bbt":
                        plt.plot(x, y - RFM_bbt, label="residual", c="tab:green")
                    elif y_type == "rad":
                        plt.plot(x, y - RFM_rad, label="residual", c="tab:green")

            # common
            plt.xlabel(x_lbl)
            plt.ylabel(y_lbl)

            plt.legend()
            plt.title(f"{keystr}")
            #        plt.show()
            plt.savefig(f"inp.values['results_fldr']/{keystr}_spc.png")
            plt.close()

        ########################################################################################
        # (optional) create plots vs iasi spectra
        ########################################################################################
        plot = True  # create a plot?
        plot_srfm = True  # plot srfm output?
        plot_iasi = True  # plot iasi spectrum?
        y_type = (
            "bbt"  # plot in radiances (W m-2 sr-1 cm) or brightness temperatures [K]
        )
        x_type = "cm-1"  # plot vs. wavenumbers [cm-1] or wavelengths [um] or [nm]

        if plot == True:
            plt.figure(figsize=(11.7, 8.4))
            plt.rcParams.update({"font.size": 12})
            plt.cla()

            x_lbl = r"Wavenumbers (cm$^{-1}$)"
            y_lbl = "Brightness temperature (K)"

            # plot SRFM
            if plot_srfm == True:
                x = model_SRFM.wvnm
                y = model_SRFM.bbt[:, 0, 0, 0]

                plt.plot(x, y, label=f"srfm", c="tab:blue")

            # plot iasi
            if plot_iasi == True:
                x = fin_grid
                y = iasi_out_spc

                plt.plot(x, y, label=f"iasi", c="tab:orange")

            # common
            plt.xlabel(x_lbl)
            plt.ylabel(y_lbl)

            plt.legend()
            plt.title(f"{keystr}")
            #        plt.show()
            plt.savefig(f"{inp.values['results_fldr']}/{keystr}_spc_vs_iasi.png")
            if inp.values["show_plots"] == True:
                plt.show()

        ########################################################################################
        # (optional) plot difference between srfm and iasi spectra
        ########################################################################################
        plot = True  # create a plot?
        y_type = (
            "bbt"  # plot in radiances (W m-2 sr-1 cm) or brightness temperatures [K]
        )
        x_type = "cm-1"  # plot vs. wavenumbers [cm-1] or wavelengths [um] or [nm]

        if plot == True:
            plt.figure(figsize=(11.7, 8.4))
            plt.rcParams.update({"font.size": 12})
            plt.cla()

            x_lbl = r"Wavenumbers (cm$^{-1}$)"
            y_lbl = "Brightness temperature (K)"

            x = model_SRFM.wvnm
            y = model_SRFM.bbt[:, 0, 0, 0] - iasi_out_spc

            plt.plot(x, y, label=f"srfm - iasi", c="tab:red")

            # common
            plt.xlabel(x_lbl)
            plt.ylabel(y_lbl)

            plt.legend()
            plt.title(f"{keystr}")
            #        plt.show()
            plt.savefig(f"{inp.values['results_fldr']}/{keystr}_diff.png")
            plt.close()

        return model_SRFM
