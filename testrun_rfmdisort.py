import numpy as np
from srfm import *
import matplotlib.pyplot as plt
import sys
import datetime
import time
import warnings
from multiprocessing import Process, Manager

start_time = time.monotonic()


###################################
# Assign some variables:
###################################

rfm_fldr = "./srfm/RFM" # where rfm is
aria_fldr = "/network/group/aopp/eodg/RGG009_GRAINGER_EODGCOMN/ARIA/" # where ARIA is

spec_res = 0.1 # model spectral resolution, [cm-1]
low_wvn = 750 # model start wavenumber (lower), [cm-1]
upp_wvn = 1250 # model end wavenumber (upper), [cm-1]

prec = "double" # choose disort precision, accepted values "single" or "double"

multiprocess = True # if True, parallelizes some calculations, 
#WARNING: EXPERIMENTAL, MAY OR MAY NOT WORK, IF UNSURE, SET FALSE

###################################
# prepare RFM driver table
###################################

# input values in the format (key:value) section name:value
rfm_inp = {}

# primary sections (mandatory), see the documentation for alternatives
rfm_inp["HDR"] = f"{str(datetime.date.today())} test run" # RFM header
rfm_inp["FLG"] = "OPT NAD SFC PRF LEV BBT RAD REJ DBL" # RFM flags
rfm_inp["SPC"] = f"{low_wvn} {upp_wvn} {spec_res}" # RFM spectral settings
rfm_inp["GAS"] = "N2 O2 CO2 O3 H2O CH4 N2O HNO3 CO NO2 N2O5 ClO HOCl ClONO2 NO HNO4 HCN NH3 F11 F12 F14 F22 CCl4 COF2 H2O2 C2H2 C2H6 OCS SO2 SF6" # RFM chemical species
rfm_inp["ATM"] = "./rfm_files/hgt_std.atm ./rfm_files/day.atm" # RFM vertical grids
rfm_inp["SEC"] = "1.0" # RFM geometry

# secondary sections, contain information required by any of the primary sections
rfm_inp["LEV"] = "./rfm_files/alts.lev" # RFM required output levels
rfm_inp["ILS"] = "./rfm_files/iasi.ils" # RFM instrument profile for convolution

# optional sections, change defaults or identify spectroscopic data files
rfm_inp["XSC"] = "/network/aopp/matin/eodg/crun/eodg/rfm/rfm_files/xsc/*.xsc" # RFM xsc files
rfm_inp["HIT"] = "/network/aopp/matin/eodg/crun/eodg/rfm/rfm_files/bin/hitran_2012.bin" # RFM hitran database
rfm_inp["REJ"] = "*   1.0E-5"

#construct rfm.drv table
rfm_functions.construct_rfm_driver_table(inp=rfm_inp,fldr=rfm_fldr)

#######################################
# prepare particle scattering properties
#######################################

# define particle properties
particle_loading = 3.23 #column particle loading, [g m-2]
#n = 1e4 # total particle concentration [cm-3]
r = 2 # mean particle radius [um]
s = 1.5 # spread, for the lognormal distribution
s_a_den = None # surf. area density, will be calculated in size dist.
v_den = None # volume density, will be calculated in size dist.
dist_type = "log_normal" # choose size distribution type
comp = "ash" # define particle type (by composition)
p_lyr_a_avg = 3.5 # avg particle layer altitude (center of layer altitude), [km]

p_lyr_thick = 1 # particle layer thickness, [km]
n = utilities.number_conc_from_particle_loading(l = particle_loading,
                                                rho = "glass",
                                                thick = p_lyr_thick,
                                                r = r)

# calculate particle layer upper and lower boundary altitude
p_lyr_u, p_lyr_l = utilities.calc_layer_extent(p_lyr_a_avg, p_lyr_thick)
#alternatively can specify layer in terms of lower and upper boundary altitude by 
# directly setting p_lyr_u and p_lyr_l here

# create particle size distribution
sd = size_distribution.create_distribution(dist_type = dist_type,
                                           n = n,
                                           r = r,
                                           s = s,
                                           surface_area_density = s_a_den,
                                           volume_density = v_den
    )

# set-up the scattering calculation
radii = 200 # number of radii in size distribution quadrature
eta = 1e-6 # value n(r) at which the size distribution upper and lower limits are set
phase_quad_N = 181 # number of quadrature points for the phase function (no. of angles)
phase_quad_type = "L" # type of phase function quadrature
radii_quad_type = "T" # type of radii size distribution quadrature
leg_coeffs = True # toggle legendre expansion coefficients for the phase function
leg_coeffs_type = "normalised" # type of Legendre polynomial expansion coefficients

RFM_wvnm = np.linspace(low_wvn,
                       upp_wvn,
                       int((upp_wvn-low_wvn) / spec_res + 1),
                       endpoint=True
                      ) # mirrors the RFM wavenumber grid
wvls = (1/RFM_wvnm)*1e4 # convert RFM wavenumber to wavelength in [um]

ewp_res = 1 # resolution of the grid to calculate optical properties at
ewp_wvnm = np.linspace(low_wvn,
                       upp_wvn,
                       int((upp_wvn-low_wvn) / ewp_res + 1),
                       endpoint=True
                     ) # grid for the optical properties
ewp_grid = (1/ewp_wvnm)*1e4 # convert ewo_wvnm wavenumbers to wavelength in [um]

#calculate optical properties either standard or parallelized, where
#   op_dict["beta_ext"] - extinction coefficient
#   op_dict["ss_alb"] - single scatter albedo
#   op_dict["phase_function"] - phase function    
#   op_dict["legendre_coefficient"] - coefficients of the Legendre expansion of the phase function, optional
if multiprocess == True:
    op_proxy_dict = Manager().dict()
    op_process = Process(target = optical_properties.ewp_hs,
                         kwargs = {"wavelengths" : ewp_grid, 
                                 "composition" : comp, 
                                 "distribution" : sd,
                                 "legendre_coefficients_flag" : leg_coeffs,
                                 "legendre_coefficients_type" : leg_coeffs_type,
                                 "radii" : radii,
                                 "eta" : eta,
                                 "phase_quad_N" : phase_quad_N,
                                 "phase_quad_type" : phase_quad_type,
                                 "radii_quad_type" : radii_quad_type,
                                 "aria" : aria_fldr,
                                 "return_dict"  :  op_proxy_dict,
                                 "multiprocess" : True
                             }
                        )
    op_process.start()
    
else:
    op_dict = optical_properties.ewp_hs(wavelengths=ewp_grid, 
                                        composition=comp, 
                                        distribution=sd,
                                        legendre_coefficients_flag=leg_coeffs,
                                        legendre_coefficients_type=leg_coeffs_type,
                                        radii=radii,
                                        eta=eta,
                                        phase_quad_N=phase_quad_N,
                                        phase_quad_type=phase_quad_type,
                                        radii_quad_type=radii_quad_type,
                                        aria=aria_fldr
                                    )

#######################################
# prepare atmospheric layer structure
#######################################

#define some requested output levels
init_lev = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0,
            14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0,
            27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 55.0, 60.0,
            65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 101.0, 102.0, 103.0, 104.0,
            105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0,
            116.0, 117.0, 118.0, 119.0, 120.0
            ]

# add upper and lower particle layer boundaries, delete any levels "within" the layer
levels, ilyr, inv_ilyr = utilities.add_lyr(old_lev = init_lev, u = p_lyr_u, l = p_lyr_l)

# write output levels file for RFM
rfm_functions.construct_rfm_output_levels_file(levels=levels,
                                               fldr=rfm_fldr
                                              )

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

# join output of the optical properties calculation (if multiprocessing)
if multiprocess == True:
    op_process.join()
    op_dict = {}
    op_dict.update(op_proxy_dict)
    
# regrid op_dict from the ewp_grid to wvls
#input("waiting for input:")
#old_dict = op_dict.copy()
op_dict, diff_dict = optical_properties.regrid(op_dict, wvls, track_diff=True,diff_type="abs")
#input("waiting for input:")
## add rfm opt output to model_RFM
model_RFM.add_rfm_opt_output(rfm_fldr, levels)

# determine cols (wavelength channels) to loop over
cols = [i for i in model_RFM.rfm_output.columns if i.startswith("dOD")]

# print current status:
print(model_RFM.status)
        
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
model_DISORT.set_maxmom(op_dict["legendre_coefficient"].shape[1]-1)
model_DISORT.set_maxcmu(16)
model_DISORT.set_maxumu(1)
model_DISORT.set_maxphi(1)
model_DISORT.set_maxulv(1)

# now the rest
model_DISORT.set_usrang(True)
model_DISORT.set_usrtau(True)
model_DISORT.set_ibcnd(0)
model_DISORT.set_onlyfl(False)
#model_DISORT.set_prnt([True, True, True, False, False])
model_DISORT.set_prnt([False, False, False, False, False])
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


# initialize disort input arrays for output variables
model_DISORT.initialize_disort_output_arrays(
    model_DISORT.disort_input["maxumu"],
    model_DISORT.disort_input["maxphi"],
    model_DISORT.disort_input["maxphi"],
)

# loop over columns, dynamically set disort input variables in each loop
for wvl_idx, (wvl,col) in enumerate(zip(wvls,cols)):
    model_DISORT.set_header(f"Now starting calculation for {col} cm-1.")

    # get wavenumber from columns name
    wvnm = RFM_wvnm[wvl_idx]
    
    # get layer optical depths from gas absorption
    tau_g = model_RFM.rfm_output[col].to_numpy()
        
    # layer optical depths from Rayleigh scattering
#    tau_R = np.zeros(shape=(tau_g.shape))
    tau_R=utilities.calc_Rayleigh_opt_depths(
                                      ps = model_RFM.rfm_output['p_lower (mbar)'].iloc[-1],
                                      pu = model_RFM.rfm_output['p_upper (mbar)'],
                                      pl = model_RFM.rfm_output['p_lower (mbar)'],
                                      l = wvnm
                                      )
    
    #particle layer optical depths (from particle scattering)
    tau_p = np.zeros( shape=( len(tau_g) ) )
    tau_p[inv_ilyr] = op_dict["beta_ext"][wvl_idx] * 1e3 * ( p_lyr_u - p_lyr_l )
    # factor 1e3 because of unit conversion (beta_ext [m-1], whereas RFM uses [km-1]
    
    #particle layer single scatter albedo
    w_p = np.zeros( shape=( len(tau_g) ) )
    w_p[inv_ilyr] = op_dict["ssalb"][wvl_idx]
    
    dtauc_tot = utilities.calc_tot_dtauc(tau_g=tau_g,
                                         tau_R=tau_R,
                                         tau_p=tau_p
                                        )
                                        
    # truncate optical depths
    threshold_od = 1e-7 # threshold at which to truncate optical depths
    idx = next(
        (index for index,value in enumerate(list(dtauc_tot)) if value > threshold_od),None)
    dtauc_tot = dtauc_tot[idx:]
    tau_g = tau_g[idx:]
    tau_R = tau_R[idx:]
    tau_p = tau_p[idx:]
    w_p = w_p[idx:]
    
    # set layer optical depths
    model_DISORT.set_dtauc_manually(dtauc=dtauc_tot)
    
    # set maxcly based on current length of tau_g
    model_DISORT.set_maxcly(len(dtauc_tot))

    ## set some maxcly-dependent DISORT input variables
    
    # set and truncate temper
    model_DISORT.set_temper_from_rfm(model_RFM)
    model_DISORT.disort_input['temper'] = model_DISORT.disort_input['temper'][idx:]

    model_DISORT.set_btemp(model_DISORT.disort_input["temper"][-1])
    model_DISORT.set_ttemp(model_DISORT.disort_input["temper"][0])
    model_DISORT.set_h_lyr(np.zeros(shape=(model_DISORT.disort_input["maxcly"] + 1)))
                   
    # set single scatter albedo
    model_DISORT.set_ssalb(tau_g=tau_g,
                           tau_R=tau_R,
                           tau_p=tau_p,
                           w_p=w_p)
    
    # calculate phase function moments for Rayleigh and particle scattering from DISORT
    pmom_R = model_DISORT.calc_pmom(iphas=2,prec=prec)


    #set phase function moments for particle scattering from Mie code
    pmom_p = np.zeros((model_DISORT.disort_input["maxmom"] + 1,
                       model_DISORT.disort_input["maxcly"])
                     )
    if ilyr > (pmom_p.shape[1]-1):
        print(f"""Scattering layer optical depth was < {threshold_od} and was
        truncated from the optical depths profile. No particle scattering at this
         wavelength.""")
    else:
        pmom_p[:,-(ilyr+1)] = op_dict["legendre_coefficient"][wvl_idx,:]
        Legendre_precision = 1/pmom_p[:,-(ilyr+1)][0] # determine the Legendre expansion precision
        pmom_p[:,-(ilyr+1)][0] = 1.0 # default the first coefficient to 0

    # calcualte the weighted sum of phase function moments      
    model_DISORT.set_pmom(pmom_R=pmom_R,
                          tau_R=tau_R,
                          w_p=w_p,
                          tau_p=tau_p,
                          pmom_p=pmom_p)
#    print(model_DISORT.disort_input["pmom"]).shape  

    # set wavenumber range for DISORT (for Planck function)
    model_DISORT.set_wvnm_range(wvnm - 0.5, wvnm + 0.5)

    # run disort input tests
    model_DISORT.test_disort_input_format()
    model_DISORT.test_disort_input_integrity()

    # initialize dictionary to store results
    model_DISORT.disort_out[col] = {}
    model_DISORT.disort_input['dtauc'][model_DISORT.disort_input['dtauc']<0] = 0
    
    # run disort itself
    model_DISORT.store_disort(col, model_DISORT.run_disort(prec=prec), bbt=True)
    
    # print current status:
#    print(model_DISORT.status)
print("Main DISORT loop finished.")

################################################################################
plot_disort = True # plot disort output?
plot_rfm = False # plot rfm output?
plot_residual = False # plot difference between rfm and diosrt?
bbt_or_rad = "bbt" # plot in radiances (W m-2 sr-1 cm) or brightness temperatures [K]

#plt.rcParams.update({'font.size': 22})

if plot_rfm == True:
    plt.ion()
    #    plt.cla()
    if bbt_or_rad == "bbt":
        filename = f"{rfm_fldr}/bbt_001000.asc"
        data = rfm_functions.read_output(filename)
        plt.plot(data["WNO"], data["SPC"], label=f"no scattering",alpha=0.9,c="tab:orange")
    elif bbt_or_rad == "rad":
        filename = f"{rfm_fldr}/rad_001000.asc"
        data = rfm_functions.read_output(filename)
        plt.plot(data["WNO"], data["SPC"]*1e-5, label=f"no scattering",alpha=0.9,c="tab:orange")
    plt.legend()
    plt.show()

if plot_disort == True:
    plt.ion()
    #    plt.cla()

    # plot the resulting brightness temperatures at the TOA in the forward direction
    wavenumbers = [
        model_DISORT.disort_out[i]["wavenumber (cm-1)"]
        for i in model_DISORT.disort_out.keys()
    ]
    
    if bbt_or_rad == "bbt": 
        y = [
            model_DISORT.disort_out[i]["uu_bbt"][0].item()
            for i in model_DISORT.disort_out.keys()
        ]
        plt.ylabel('Brightness temperature (K)')
    elif bbt_or_rad == "rad":
        y = [
            model_DISORT.disort_out[i]["uu"][0].item()
            for i in model_DISORT.disort_out.keys()
        ]
        plt.ylabel("Radiance (W m-2 sr-1 cm)")
        
    plt.plot(wavenumbers, y, label=f"ash layer: {p_lyr_a_avg - 0.5} - {p_lyr_a_avg + 0.5} km", c="tab:blue")
    plt.xlabel(r"Wavenumbers (cm$^{-1}$)")

    plt.legend()
    plt.show()

#plt.title("tau_g,temper=rfm; tau_R=0,ssalb=True;plank=true,rfm;res=0.25cm-1")
#plt.savefig('2025_02_17_test.pdf')

if plot_residual == True:
    if plot_disort != True and plot_rfm != True:    
        print(f"Can't plot residual difference without data.")
    else:
        plt.ion()
        if bbt_or_rad == "bbt":
            plt.plot(wavenumbers, y-(data["SPC"]), label = "residual")
        elif bbt_or_rad == "rad":
            plt.plot(wavenumbers, y-(data["SPC"]*1e-5), label = "residual")
        plt.legend()
        plt.show()
    
end_time = time.monotonic()
print(f"Total run time: {datetime.timedelta(seconds=end_time - start_time)}")
