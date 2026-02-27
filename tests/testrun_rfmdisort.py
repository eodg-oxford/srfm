import numpy as np
from srfm import *
import matplotlib.pyplot as plt
import sys
import datetime
import time
import warnings
from multiprocessing import Process, Manager
from oxford_colours.oxford_colours import clrs as oxclrs

start_time = time.monotonic()

########################################################################################
# Assign some variables:
########################################################################################

rfm_fldr = "./srfm/RFM"  # where rfm is
aria_fldr = "/network/group/aopp/eodg/RGG009_GRAINGER_EODGCOMN/ARIA/"  # where ARIA is

prec = "double"  # choose disort precision, accepted values "single" or "double"

multiprocess = True  # if True, parallelizes some calculations,
# WARNING: EXPERIMENTAL, MAY OR MAY NOT WORK, IF UNSURE, SET FALSE

########################################################################################
# specify spectral calculation grid
########################################################################################
spec_res = 0.005  # model spectral resolution,[spec_units]
low_spc = 637  # model start wavenumber (lower), [spec_units]
upp_spc = 2763  # model end wavenumber (upper), [spec_units]
spec_units = "cm-1"  # accepted values "cm-1", "um", "nm"

RFM_wvnm, wvls = utilities.calc_grids(low_spc, upp_spc, spec_res, spec_units)

rfm_grid_fname = rfm_functions.construct_rfm_grid_file(
    RFM_wvnm, filename="grid.spc", rfm_fldr=rfm_fldr
)
########################################################################################
# define an atmospheric scattering layer
########################################################################################
scat_lyrs = (
    {}
)  # dictionary of scattering layers, key - layer name, value - Layer() object

scat_lyrs_inputs = {
    "Mie_1": {
        "name": "Mie_1",  # layer identifier/name
        "low_spc": low_spc,  # spectral claculation grid lower limit
        "upp_spc": upp_spc,  # spectral calculation grid upper limit
        "res": 1,  # specral calculation grid resolution
        "spec_units": spec_units,  # spectral calculation grid units
        "mass_loading": 3.23,  # column particle loading, [g m-2]
        "rho": "mineral",  # particle density, can be number or one of permitted strings
        "n": None,  # total particle concentration [cm-3]
        "r": 2,  # mean particle radius [um]
        "s": 1.5,  # spread, for the lognormal distribution
        "s_a_den": None,  # surf. area density, will be calculated in size dist.
        "v_den": None,  # volume density, will be calculated in size dist.
        "dist_type": "log_normal",  # choose size distribution type
        "comp": "ash",  # define particle type (by composition)
        "center_alt": 3.5,  # avg particle layer altitude (center of layer altitude), [km]
        "thick": 0.1,  # particle layer thickness, [km]
        "alt_upp": None,  # particle layer upper boundary altitude, [km]
        "alt_low": None,  # particle layer lower boundary altitude, [km],
        "radii": 200,  # number of radii in size distribution quadrature
        "eta": 1e-6,  # value n(r) at which the size distribution upper and lower limits are set
        "phase_quad_N": 181,  # number of quadrature points for the phase function (no. of angles)
        "phase_quad_type": "L",  # type of phase function quadrature
        "radii_quad_type": "T",  # type of radii size distribution quadrature
        "leg_coeffs": True,  # toggle legendre expansion coefficients for the phase function
        "leg_coeffs_type": "normalised",  # type of Legendre coeffs, normalised or regular
        "aria": aria_fldr,  # where ARIA is
        "multiprocess": True,  # type of Legendre polynomial expansion coefficients
    },
    #        "Mie_2" : {
    #            "name" : "Mie_2", # layer identifier/name
    #            "low_spc" : low_spc, # spectral claculation grid lower limit
    #            "upp_spc" : upp_spc, # spectral calculation grid upper limit
    #            "res" : 1, # specral calculation grid resolution
    #            "spec_units" : spec_units, # spectral calculation grid units
    #            "mass_loading" : 2.23, #column particle loading, [g m-2]
    #            "rho" : "glass", # particle density, can be number or one of permitted strings
    #            "n" : None, # total particle concentration [cm-3]
    #            "r" : 2,  # mean particle radius [um]
    #            "s" : 1.5, # spread, for the lognormal distribution
    #            "s_a_den" : None, # surf. area density, will be calculated in size dist.
    #            "v_den" : None, # volume density, will be calculated in size dist.
    #            "dist_type" : "log_normal", # choose size distribution type
    #            "comp" : "ice", # define particle type (by composition)
    #            "center_alt" : 8, # avg particle layer altitude (center of layer altitude), [km]
    #            "thick" : 1, # particle layer thickness, [km]
    #            "alt_upp" : None, # particle layer upper boundary altitude, [km]
    #            "alt_low" : None, # particle layer lower boundary altitude, [km],
    #            "radii" : 181, # number of radii in size distribution quadrature
    #            "eta" : 1e-6, # value n(r) at which the size distribution upper and lower limits are set
    #            "phase_quad_N" : 200, # number of quadrature points for the phase function (no. of angles)
    #            "phase_quad_type" : "L", # type of phase function quadrature
    #            "radii_quad_type" : "T", # type of radii size distribution quadrature
    #            "leg_coeffs" : True, # toggle legendre expansion coefficients for the phase function
    #            "leg_coeffs_type" : "normalised", # type of Legendre coeffs, normalised or regular
    #            "aria" : aria_fldr, # where ARIA is
    #            "multiprocess" : True, # type of Legendre polynomial expansion coefficients
    #        },
}

for lyr in scat_lyrs_inputs.keys():
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

# write output levels file for RFM
rfm_out_lvl_fname = "alts.lev"
rfm_functions.construct_rfm_output_levels_file(
    levels=levels, fldr=rfm_fldr, fname=f"{rfm_out_lvl_fname}"
)

########################################################################################
# prepare RFM driver table
########################################################################################

# input values in the format (key:value) section name:value
rfm_inp = {}

# primary sections (mandatory), see the documentation for alternatives
rfm_inp["HDR"] = f"{str(datetime.date.today())} test run"  # RFM header
rfm_inp["FLG"] = "OPT NAD SFC PRF LEV RAD DBL"  # RFM flags
rfm_inp["SPC"] = f"{rfm_grid_fname}"  # RFM spectral settings
rfm_inp["GAS"] = (
    """N2 O2 CO2 O3 H2O CH4 N2O HNO3 CO NO2 N2O5 ClO HOCl ClONO2 NO HNO4 HCN NH3 F11 F12 F14 F22 CCl4 COF2 H2O2 C2H2 C2H6 OCS SO2 SF6"""
)
# RFM chemical species
rfm_inp["ATM"] = "./rfm_files/hgt_std.atm ./rfm_files/day.atm"  # RFM vertical grids
rfm_inp["SEC"] = "1.0"  # RFM geometry

# secondary sections, contain information required by any of the primary sections
rfm_inp["LEV"] = f"./rfm_files/{rfm_out_lvl_fname}"  # RFM required output levels
rfm_inp["ILS"] = "./rfm_files/iasi.ils"  # RFM instrument profile for convolution

# optional sections, change defaults or identify spectroscopic data files
rfm_inp["XSC"] = "/network/aopp/matin/eodg/crun/eodg/rfm/rfm_files/xsc/*.xsc"
# RFM xsc files
rfm_inp["HIT"] = "/network/aopp/matin/eodg/crun/eodg/rfm/rfm_files/bin/hitran_2012.bin"
# RFM hitran database
# rfm_inp["REJ"] = "*   1.0E-5"

# construct rfm.drv table
rfm_functions.construct_rfm_driver_table(inp=rfm_inp, fldr=rfm_fldr)

########################################################################################
# run RFM
########################################################################################

# compile rfm
rfm_functions.compile_rfm(rfm_fldr)

# initialize RFM model class
model_RFM = forward_model.RFM()

# print current status:
print(model_RFM.status)

# run rfm
model_RFM.run_rfm(rfm_fldr)

# print current status:
print(model_RFM.status)

# add output from optical properties calculation (placed here, because if calculations
# run in parallel processes, here is the place they join the main process.)
for lyr in scat_lyrs.keys():
    scat_lyrs[lyr].add_op_calc_output()

    # interpolate layer optical properties
    scat_lyrs[lyr].regrid(wvls, track_diff=True)
    scat_lyrs[lyr].calc_tau()

## add rfm opt output to model_RFM
model_RFM.add_rfm_opt_output(rfm_fldr, levels)

# determine cols (wavelength channels) to loop over
cols = [i for i in model_RFM.rfm_output.columns if i.startswith("dOD")]

# print current status:
print(model_RFM.status)

########################################################################################
# prepare DISORT common variables
########################################################################################

# initialize DISORT model class
model_DISORT = forward_model.DISORT()

# check if number of columns and wavelengths match
if len(cols) != len(wvls):
    raise ValueError("Number of RFM and scattering wavelengths don't match.")

# set disort_input parameters common to all loop iterations
# these need to be set first:
nmom = 0
for lyr in scat_lyrs.keys():
    if (scat_lyrs[lyr].legendre_coefficient.shape[1] - 1) > nmom:
        nmom = scat_lyrs[lyr].legendre_coefficient.shape[1] - 1

model_DISORT.set_maxcmu(16)

model_DISORT.set_maxmom(nmom)
if nmom < model_DISORT.disort_input["maxcmu"]:
    model_DISORT.set_maxmom(model_DISORT.disort_input["maxcmu"])


model_DISORT.set_maxumu(4)
model_DISORT.set_maxphi(2)
model_DISORT.set_maxulv(3)

# now the rest
model_DISORT.set_usrang(True)
model_DISORT.set_usrtau(True)
model_DISORT.set_ibcnd(0)
model_DISORT.set_onlyfl(False)
# model_DISORT.set_prnt([True, True, True, False, False])
model_DISORT.set_prnt([False, False, False, False, False])
model_DISORT.set_plank(True)
model_DISORT.set_lamber(True)
model_DISORT.set_deltamplus(False)
model_DISORT.set_do_pseudo_sphere(False)
model_DISORT.set_utau([0, 0.1, 0.2])
model_DISORT.set_umu0(0.1)
model_DISORT.set_phi0(0)
model_DISORT.set_umu([-1, -0.1, 0.1, 1])
model_DISORT.set_phi([0, 1])
model_DISORT.set_fbeam(0)
model_DISORT.set_fisot(0)
model_DISORT.set_albedo(0)


model_DISORT.set_temis(1)
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
        shape=(model_DISORT.disort_input["maxumu"], model_DISORT.disort_input["maxphi"])
    )
)
model_DISORT.set_bemst(np.zeros(shape=(int(model_DISORT.disort_input["maxcmu"] / 2))))
model_DISORT.set_emust(np.zeros(shape=(model_DISORT.disort_input["maxumu"])))
model_DISORT.set_accur(0)


# initialize disort input arrays for output variables from a single run
model_DISORT.initialize_disort_output_arrays()

########################################################################################
# prepare SRFM common variables
########################################################################################
model_SRFM = forward_model.SRFM()
model_SRFM.set_wvnm(RFM_wvnm)
model_SRFM.set_wvls(wvls)
model_SRFM.initialize_srfm_output_arrays_from_disort(model_DISORT)

# loop over columns, dynamically set disort input variables in each loop
for wvl_idx, (wvnm, wvl, col) in enumerate(zip(RFM_wvnm, wvls, cols)):
    model_DISORT.set_header(f"Now starting calculation for {col} cm-1.")
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
        (index for index, value in enumerate(list(dtauc_tot)) if value > threshold_od),
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
    model_DISORT.set_h_lyr(np.zeros(shape=(model_DISORT.disort_input["maxcly"] + 1)))

    # set single scatter albedo
    model_DISORT.set_ssalb(tau_g=tau_g, tau_R=tau_R, tau_p=tau_p, w_p=w_p)

    # calculate phase function moments for Rayleigh and particle scattering from DISORT
    pmom_R = model_DISORT.calc_pmom(iphas=2, prec=prec)

    # set phase function moments for particle scattering from Mie code
    pmom_p = np.zeros(
        (model_DISORT.disort_input["maxmom"] + 1, model_DISORT.disort_input["maxcly"])
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
            Legendre_precision = 1 / pmom_p[0, track_lyr_local.index(lyr)]
            pmom_p[0, track_lyr_local.index(lyr)] = 1.0

    # calculate the weighted sum of phase function moments
    model_DISORT.set_pmom(
        pmom_R=pmom_R, tau_R=tau_R, w_p=w_p, tau_p=tau_p, pmom_p=pmom_p
    )
    #    print(model_DISORT.disort_input["pmom"]).shape

    # set wavenumber range for DISORT (for Planck function)
    model_DISORT.set_wvnm_range(wvnm - 0.5, wvnm + 0.5)

    # run disort input tests
    model_DISORT.test_disort_input_format()
    model_DISORT.test_disort_input_integrity()

    # run disort
    model_DISORT.run_disort(prec=prec)

    # store result in SRFM()
    model_SRFM.store_disort_result(model_DISORT, wvl_idx)

    # print current status:
#    print(model_DISORT.status)
print("Main DISORT loop finished.")

model_SRFM.convolve_with_iasi(f"{rfm_fldr}/rfm_files/iasi.ils")

model_SRFM.calc_bbt()

model_SRFM.interp(np.linspace(637, 2763, int((2763 - 637) / 0.25 + 1)))

########################################################################################
# (optional) create plots
########################################################################################
plot = True  # create a plot?
plot_srfm = True  # plot disort output?
plot_rfm = True  # plot rfm output?
plot_residual = False  # plot difference between rfm and diosrt?
y_type = "bbt"  # plot in radiances (W m-2 sr-1 cm) or brightness temperatures [K]
x_type = "cm-1"  # plot vs. wavenumbers [cm-1] or wavelengths [um] or [nm]

if plot == True:
    # plt.rcParams.update({'font.size': 22})
    plt.ion()
    plt.cla()

    # determine x:
    if x_type == "cm-1":
        x = RFM_wvnm
        x_lbl = r"Wavenumbers (cm$^{-1}$)"
    elif x_type == "um":
        x = wvls
        x_lbl = r"Wavelength ($\mu$m)"
    elif x_type == "nm":
        x = wvls * 1e3
        x_lbl = "Wavelength (nm)"

    # determine y:
    if y_type == "bbt":
        y_lbl = "Brightness temperature (K)"
    elif y_type == "rad":
        y_lbl = r"Radiance (W m$^{-2}$ sr$^{-1}$ cm)"

    # plot RFM
    if plot_rfm == True:
        filename = f"{rfm_fldr}/rad_001000.asc"
        data = rfm_functions.read_output(filename)
        RFM_rad = data["SPC"] * 1e-5

        if y_type == "bbt":
            RFM_bbt = utilities.convert_spectral_radiance_to_bbt(RFM_rad, RFM_wvnm)
            plt.plot(
                x, RFM_bbt, label=f"no scattering", alpha=1, c=oxclrs["Oxford blue"]
            )
        elif y_type == "rad":
            plt.plot(
                x, RFM_rad, label=f"no scattering", alpha=0.9, c=oxclrs["Oxford blue"]
            )

    # plot SRFM
    if plot_srfm == True:
        if y_type == "bbt":
            y = model_SRFM.bbt[:, 3, 0, 0]
        elif y_type == "rad":
            y = model_SRFM.uu[:, 3, 0, 0]

        plt.plot(x, y, label=f"ash layer 3-4 km", c=oxclrs["Oxford pink"])

    # plot residual
    if plot_residual == True:
        if plot_srfm != True and plot_rfm != True:
            print(f"Can't plot residual without data.")
        else:
            plt.ion()
            if y_type == "bbt":
                plt.plot(
                    x, y - RFM_bbt, label="residual", c=oxclrs["Oxford lime green"]
                )
            elif y_type == "rad":
                plt.plot(
                    x, y - RFM_rad, label="residual", c=oxclrs["Oxford lime green"]
                )

    # common
    plt.xlabel(x_lbl)
    plt.ylabel(y_lbl)

    #    plt.plot(model_SRFM.wvnm,model_SRFM.bbt_unconvolved[:,3,0,0],label="unconvolved",c="tab:red")
    #    plt.plot(model_SRFM.wvnm,model_SRFM.bbt[:,3,0,0],label="convolved",c="tab:green")

    plt.legend()
    plt.show()

end_time = time.monotonic()
print(f"Total run time: {datetime.timedelta(seconds=end_time - start_time)}")
