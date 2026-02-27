import numpy as np
from srfm import *
import matplotlib.pyplot as plt
import sys
import datetime
import time

start_time = time.monotonic()


###################################
# Assign some variables:
###################################

rfm_fldr = "./srfm/RFM"
disort_fldr = "./srfm/DISORT"

spec_res = 1  # model spectral resolution, [cm-1]
low_wvn = 760  # model start wavenumber (lower), [cm-1]
upp_wvn = 2500  # model end wavenumber (upper), [cm-1]

###################################
# prepare RFM driver table
###################################

# input values in the format (key:value) section name:value
rfm_inp = {}

# primary sections (mandatory), see the documentation for alternatives
rfm_inp["HDR"] = f"{str(datetime.date.today())} test run"
rfm_inp["FLG"] = "OPT NAD SFC PRF LEV BBT"
rfm_inp["SPC"] = f"{low_wvn} {upp_wvn} {spec_res}"
rfm_inp["GAS"] = (
    "N2 O2 CO2 O3 H2O CH4 N2O HNO3 CO NO2 N2O5 ClO HOCl ClONO2 NO HNO4 HCN NH3 F11 F12 F14 F22 CCl4 COF2 H2O2 C2H2 C2H6 OCS SO2 SF6"
)
rfm_inp["ATM"] = "./rfm_files/hgt_std.atm ./rfm_files/day.atm"
rfm_inp["SEC"] = "1.0"

# secondary sections, contain information required by any of the primary sections
rfm_inp["LEV"] = "./rfm_files/alts.lev"
rfm_inp["ILS"] = "./rfm_files/iasi.ils"

# optional sections, change defaults or identify spectroscopic data files
rfm_inp["XSC"] = "/network/aopp/matin/eodg/crun/eodg/rfm/rfm_files/xsc/*.xsc"
rfm_inp["HIT"] = "/network/aopp/matin/eodg/crun/eodg/rfm/rfm_files/bin/hitran_2012.bin"

# construct rfm.drv table
rfm_functions.construct_rfm_driver_table(inp=rfm_inp, fldr=rfm_fldr)

#######################################
# run RFM
#######################################

# compile rfm
rfm_functions.compile_rfm(rfm_fldr)

# initialize RFM model class
model_RFM = forward_model.RFM()

# print current status:
print(model_RFM.status)

# run rfm
model_RFM.run_rfm(rfm_fldr)

# add rfm opt output to model_RFM
model_RFM.add_rfm_opt_output(rfm_fldr)

# determine cols (wavelength channels) to loop over
cols = [i for i in model_RFM.rfm_output.columns if i.startswith("dOD")]

# print current status:
print(model_RFM.status)

#######################################
# prepare particle scattering properties
#######################################

# define particle properties
n = 10  # total particle concentration [cm-3]
r = 2  # mean particle radius [um]
s = 1.5  # spread, for the lognormal distribution
dist_type = "log_normal"  # choose size distribution type
comp = "ash"  # define particle type (by composition)
h_p = 24  # avg particle layer height, [km]

# create particle size distribution
sd = size_distribution.create_distribution(dist_type=dist_type, n=n, r=r, s=s)

# set-up the scattering calculation
angles = np.linspace(0, 180, 181)  # angles to calculate scattering at, [degrees]
leg_coeffs = True  # toggle legendre expansion coefficients for the phase function
RFM_wvnm = np.linspace(
    low_wvn, upp_wvn, int((upp_wvn - low_wvn) / spec_res + 1), endpoint=True
)
wvls = (1 / RFM_wvnm) * 1e4  # convert RFM wavenumber to wavelength in [um]

# calculate the optical properties, where
# e - extinction coefficient
# w - single scatter albedo
# p - phase function
# l - coefficients of the Legendre expansion of the phase function
e, w, p, l, c, cw = optical_properties.ewp_hs(
    wvls, comp, sd, legendre_coefficients_flag=leg_coeffs
)
# print(l.shape)

#######################################
# prepare and run DISORT
#######################################

# initialize DISORT model class
model_DISORT = forward_model.DISORT()


# check if number of columns and wavelengths match
if len(cols) != len(wvls):
    raise ValueError("Number of RFM and scattering wavelengths don't match.")

# set disort_input parameters common to all loop iterations
# these need to be set first:

model_DISORT.set_maxmom(l.shape[1] - 1)
model_DISORT.set_maxcmu(16)
model_DISORT.set_maxumu(1)
model_DISORT.set_maxphi(1)
model_DISORT.set_maxulv(1)

# now the rest
model_DISORT.set_usrang(True)
model_DISORT.set_usrtau(True)
model_DISORT.set_ibcnd(0)
model_DISORT.set_onlyfl(False)
model_DISORT.set_prnt([True, True, True, False, True])
model_DISORT.set_plank(True)
model_DISORT.set_lamber(True)
model_DISORT.set_deltamplus(False)
model_DISORT.set_do_pseudo_sphere(False)
model_DISORT.set_utau([0])
model_DISORT.set_umu0(0.1)
model_DISORT.set_phi0(0)
model_DISORT.set_umu([1])
model_DISORT.set_phi([0])
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
model_DISORT.set_header("")

# initialize disort input arrays for output variables
model_DISORT.initialize_disort_output_arrays(
    model_DISORT.disort_input["maxumu"],
    model_DISORT.disort_input["maxphi"],
    model_DISORT.disort_input["maxphi"],
)

# loop over columns, dynamically set disort input variables in each loop
for wvl_idx, (wvl, col) in enumerate(zip(wvls, cols)):
    print(col)

    ## set number of computational layers (truncate zeros in opt. depths)

    # get layer optical depths from gas absorption
    tau_g = model_RFM.rfm_output[col].to_numpy()
    #    model_DISORT.set_temper(np.linspace(250,280,len(tau_g)+1))

    # truncate tau_g
    idx = next((index for index, value in enumerate(list(tau_g)) if value > 1e-4), None)
    tau_g = tau_g[idx:]

    # set maxcly based on current length of tau_g
    model_DISORT.set_maxcly(len(tau_g))

    ## set some maxcly-dependent DISORT input variables

    # set and truncate temper
    model_DISORT.set_temper_from_rfm(model_RFM)
    model_DISORT.disort_input["temper"] = model_DISORT.disort_input["temper"][idx:]

    model_DISORT.set_btemp(model_DISORT.disort_input["temper"][-1])
    model_DISORT.set_ttemp(model_DISORT.disort_input["temper"][0])
    model_DISORT.set_h_lyr(np.zeros(shape=(model_DISORT.disort_input["maxcly"] + 1)))

    # get wavenumber from columns name
    wvnm = float(col[col.rfind("_") + 1 :])
    # layer optical depths from Rayleigh scattering
    #    tau_R = np.zeros(shape=(model_DISORT.disort_input["maxcly"]))
    tau_R = utilities.calc_Rayleigh_opt_depths(
        ps=model_RFM.rfm_output["p_lower (mbar)"].iloc[-1],
        pu=model_RFM.rfm_output["p_upper (mbar)"],
        pl=model_RFM.rfm_output["p_lower (mbar)"],
        l=wvnm,
    )
    tau_R = tau_R[idx:]

    # determine layer closest to aerosol height
    hgts = list(model_RFM.rfm_output["h_avg (km)"])[idx:]
    h_diff = [abs(h_p - i) for i in hgts]
    h_idx = h_diff.index(min(h_diff))

    # layer optical depths from aerosol scattering
    tau_p = np.zeros(shape=(model_DISORT.disort_input["maxcly"]))
    tau_p[h_idx] = e[wvl_idx] * (
        model_RFM.rfm_output["h_upper (km)"][h_idx]
        - model_RFM.rfm_output["h_lower (km)"][h_idx]
    )
    #    tau_p = tau_p[idx:]

    # layer aerosol single scatter albedo, currently not set
    w_p = np.zeros(shape=(model_DISORT.disort_input["maxcly"]))
    #    print(w_p.shape)
    w_p[h_idx] = w[wvl_idx]
    #    w_p = w_p[idx:]

    # set layer optical depths
    model_DISORT.set_dtauc(tau_g=tau_g, tau_R=tau_R, tau_p=tau_p)

    # set single scatter albedo
    model_DISORT.set_ssalb(tau_g=tau_g, tau_R=tau_R, tau_p=tau_p, w_p=w_p)

    # calculate phase function moments for Rayleigh and particle scattering
    pmom_R = model_DISORT.calc_pmom(iphas=2)
    #    pmom_p = model_DISORT.calc_pmom(iphas=6)
    print(model_DISORT.disort_input["maxmom"] + 1)
    print(model_DISORT.disort_input["maxcly"])
    pmom_p = np.zeros(
        (model_DISORT.disort_input["maxmom"] + 1, model_DISORT.disort_input["maxcly"])
    )
    #    pmom_p[:,h_idx] = l[wvl_idx]

    # set phase function moments
    model_DISORT.set_pmom(
        pmom_R=pmom_R, tau_R=tau_R, w_p=w_p, tau_p=tau_p, pmom_p=pmom_p
    )
    #    print(model_DISORT.disort_input["pmom"])
    model_DISORT.set_wvnm_range(wvnm - 0.5, wvnm + 0.5)

    # run disort input tests
    model_DISORT.test_disort_input_format()
    model_DISORT.test_disort_input_integrity()

    # initialize dictionary to store results
    model_DISORT.disort_out[col] = {}
    model_DISORT.disort_input["dtauc"][model_DISORT.disort_input["dtauc"] < 0] = 0

    # run disort itself
    model_DISORT.store_disort(col, model_DISORT.run_disort(), bbt=True)

    # print current status:
    print(model_DISORT.status)


################################################################################
plot_disort = True
plot_rfm = False
plot_residual = False

if plot_disort == True:
    plt.ion()
    #    plt.cla()

    # plot the resulting brightness temperatures at the TOA in the forward direction
    wavenumbers = [
        model_DISORT.disort_out[i]["wavenumber (cm-1)"]
        for i in model_DISORT.disort_out.keys()
    ]
    y = [
        model_DISORT.disort_out[i]["uu_bbt"][0].item()
        for i in model_DISORT.disort_out.keys()
    ]

    plt.plot(wavenumbers, y, label="opt. depths truncated below 1e-4, res. 0.1 cm-1")
    plt.xlabel(r"Wavenumbers (cm$^{-1}$)")
    plt.ylabel("Brightness temperature (K)")
    #    plt.ylabel("Radiance (W m-2 sr-1 cm)")

    plt.legend()
    plt.show()

if plot_rfm == True:
    plt.ion()
    #    plt.cla()
    filename = "./oxharp/RFM/rad_01000.asc"
    data = rfm_functions.read_output(filename)
    plt.plot(
        data["WNO"],
        data["SPC"] * 1e-5,
        label=f"RFM_{filename}, res. 0.1 cm-1",
        alpha=0.5,
    )  # data['SPC']*1e-5 if plotting radiances to convert to W m-2 sr-1 cm
    #    if data['NPNT'] > 0:
    #        plt.xlabel(r'Wavenumber (cm$^{-1}$)')
    #    else:
    #        plt.xlabel('Frequency (GHz)')
    #    plt.ylabel(data['LABSPC'])
    plt.legend()
    plt.show()

# plt.title("tau_g,temper=rfm; tau_R=0,ssalb=True;plank=true,rfm;res=0.25cm-1")
# plt.savefig('2024_12_06_test_2_bbt.pdf')

if plot_residual == True:
    if plot_disort != True and plot_rfm != True:
        print(f"Can't plot residual difference without data.")
    else:
        plt.ion()
        plt.plot(wavenumbers, y - (data["SPC"] * 1e-5), label="residual")
        plt.legend()
        plt.show()

end_time = time.monotonic()
print(datetime.timedelta(seconds=end_time - start_time))
