"""This code defines one function to run srfm.

    Note that the SRFM can be run interactively (i.e. used as a package, use the
    inside of run_srfm() as and example), or with a driver table.

- Name: main
- Parent package: srfm
- Author: Antonin Knizek
- Contributors:
- Date: 25 Nov 2025
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
from .RFM import rfm_py
from . import rfm_helper
from mergedeep import merge
from netCDF4 import Dataset
import json
import copy

@utilities.show_runtime
def run_srfm(inp):
    """Main function that runs srfm.
    
    Generic srfm  run.

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
    # set final grid to interpolate to
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
    
    if "date" in inp.values:
        if not isintance(inp.values["date"], tuple) or not len(inp.values["date"]) == 3:
            raise ValueError("date must be a tuple of len 3 (see datetime.datetime.")
        else:
            date = inp.values["date"]
    else:
        date = datetime.datetime(2025,3,23) # default date, approx spring equinox, avg Earth-Sun dist
    year_day = date.timetuple().tm_yday
    
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
    levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0,
        14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 27.5,
        30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 55.0, 57.5, 60.0, 62.5,
        65.0, 67.5, 70.0, 75.0, 80.0, 82.5, 85.0, 87.5, 90.0, 95.0, 100.0, 101.0, 102.0, 
        103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 110.5, 111.0, 111.5,
        112.0, 112.5, 113.0, 113.5, 114.0, 114.5, 115.0, 115.5, 116.0, 116.5, 117.0,
        117.5, 118.0, 118.5, 119.0, 119.5, 120.0,
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
    # precalculate angles    
    zen_rad = np.deg2rad(inp.values["zen"]) # zenith angle to rad
    zen_cos = np.cos(zen_rad).item() # cosine of zenith angle
    zen_sec = 1/ zen_cos # secant of zenith angle
    
    # RFM global config
    rfm_config = inp.values["rfm_config"]

    # RFM driver table
    driver_inputs = inp.values["driver_inputs"]
    driver_inputs["lev"] = tuple(str(val) for val in levels)
    driver_inputs["tangent"] = (str(zen_sec),)

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
    model_RFM.status = "RFM completed"
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
    if inp.values["sun"] == True:
        # load solar spectrum from file
        solar_spc, solar_spc_wvnm = utilities.load_solar_spectrum_Gueymard20018()
        solar_spc = np.interp(
            RFM_wvnm, solar_spc_wvnm[::-1], solar_spc[::-1]
        )  # interpolate to RFM_wvnm (the calculation grid)

        # scale with year day (different Sun-Earth distance throughout the year
        # the original specturm is for 1 AU
        solar_spc = utilities.scale_solar_spectrum(solar_spc, year_day)

        # get incoming solar beam polar angle for DISORT
        solar_zen_deg = inp.values["sza"]  # solar zenith angle [degrees]
        solar_zen_rad = np.deg2rad(solar_zen_deg)  # solar zenith angle [rad]
        solar_zen_cos = np.cos(
            solar_zen_rad
        ).item()  # cosine of the solar zenith angle, UMU0

        model_DISORT.set_umu0(solar_zen_cos)
        model_DISORT.set_phi0(inp.values["saa"])
    else:
        model_DISORT.set_umu0(1)
        model_DISORT.set_phi0(0)

    # Angles
    model_DISORT.set_phi([inp.values["azi"]])
    model_DISORT.set_umu([zen_cos])
    
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
    pct = [1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]  # percent completed to report
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
#        model_DISORT.test_disort_input_format()
#        model_DISORT.test_disort_input_integrity()
        # These tests are currently disabled, pending review. It turns out they test
        # inputs, inlc. arrays, element by element, taking about 25% of the total 
        # runtime.
        # TODO: options:
        # 1. Remove the tests, potentially unsafe, or rather may be less explanatory
        # when the run fails.
        # 2. Split the tests (i.e. test elements that do not change between runs
        # separately and then test within the loop only the elements that change from
        # iteration to iteration.
        # 3. Move the test behind a debug flag (add the flag to the driver table).
        # 4. Replace element-wise checks with vectorized dtype assertions
        # (e.g. for dtauc do np.asarray(dtauc).dtype / np.issubdtype... 

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
    if inp.values["convolve_iasi"] == True:
        model_SRFM.convolve_with_iasi(inp.values["iasi_ils"])

    # interpolate resulting bbt and radiances to final grid
    model_SRFM.interp(fin_grid)

    # calculate brightness temperature for the final spectrum
    model_SRFM.calc_bbt()

    ########################################################################################
    # (optional) save spectrum to file(s)
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
        y_type = (
            "bbt"  # plot in radiances (W m-2 sr-1 cm) or brightness temperatures [K]
        )
        x_type = "cm-1"  # plot vs. wavenumbers [cm-1] or wavelengths [um] or [nm]

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


        # determine x:
        if x_type == "cm-1":
            x = model_SRFM.wvnm
        elif x_type == "um":
            x = model_SRFM.wvls
        elif x_type == "nm":
            x = model_SRFM.wvls * 1e3
        
        # determine y:
        if y_type == "bbt":
            y = model_SRFM.bbt[:, 0, 0, 0]
        elif y_type == "rad":
            y = model_SRFM.uu[:, 0, 0, 0]

        plt.plot(x, y, label=f"SRFM", c="tab:blue")

        # common
        plt.xlabel(x_lbl)
        plt.ylabel(y_lbl)

        plt.legend()
        plt.savefig(f"{inp.values['results_fldr']}/base_plot.png")
        if inp.values["show_plots"] == True:
            plt.show()
        else:
            plt.close()

        return model_SRFM
