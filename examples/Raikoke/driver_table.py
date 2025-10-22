import os
from srfm import *

### Ancillary information
# dict, can contain any information that needs to be passed to the forward model,
# but that does not in any way constitute a part of the retrieval

b = {
    ## Measurement (IASI L1c preprocessed file, angles, etc.
    "iasi_spc_fldr" : "/network/group/aopp/eodg/RGG008_GRAINGER_IASIVOLC/knizek/raikoke_ssaccs/spectra_w_iasi_L2", # folder with preprocessed IASI L1c spectra (.pkl files) from my preprocessor
    "iasi_fl" : "hlat_20190701_A.pkl", # which iasi file to use from my preprocessed files, corresponds to one pre-processed orbit file
    "px" : 8, # which pixel to use
    
    ## Paths
    "results_fldr" : f"{abs_path}/results", # where to store results from one calculation
    "hitbin" : "/network/aopp/matin/eodg/crun/eodg/rfm/rfm_files/bin/hitran_2024_10_8.bin", # hitran database bin file path for RFM
    "xsc" : "/network/aopp/matin/eodg/crun/eodg/rfm/rfm_files/xsc/*.xsc", # xsc files file paths for RFM (molecular absorption cross-sections)
    "ils" : f"{abs_path}/iasi.ils", # instrument line shape file
    "nedt" : f"{abs_path}/iasi_nedt.txt", # iasi noise equivalent delta temperature file
    
    ## Grids
    # note that calculation output wavenumber gridd are separate objects, i.e. you can run a high resolution calculation and then interpolate to lower resolution (e.g. to match IASI gridpoints)
    "fin_wvnmlo" : fin_wvnmlo,                                    # lower wavenumber for the final spectrum output, cm-1, int/float
    "fin_wvnmhi" : fin_wvnmhi,                                   # upper wavenumber for the final spectrum output, cm-1, int/float
    "fin_res" : fin_res,                                      # resolution of the final spectrum output, cm-1

    "spc_res" : spc_res,                                      # main calculation spectral resolution
    "spc_wvnmlo" : fin_wvnmlo - 3,                         # main calculation upper limit, if convolving with ils, must be in cm-1 and at least 2 cm-1 greater than fin_wvnmhi
    "spc_wvnmhi" : fin_wvnmhi + 3,                         # main calculation upper limit, if convolving with ils, must be in cm-1 and at least 2 cm-1 greater than fin_wvnmhi
    "spc_units" : spc_units,                                  # main calculation units, cm-1, um or nm
    
    ## Atmospheric profiles
    "atm" : f"{abs_path}/day.atm", # absolute or relative path to .atm file which contains atmospheric profiles, uses the RFM .atm file format
    
    ## DISORT variables
    "fisot" : 0, # isotropic radiation at the top of the atmosphere
    "albedo" : 0, # bottom boundary albedo
    "temis" : 1, # top boundary emissivity
    "earth_radius" : 6371, # Earth radius
    "nmom" : 3, # maximum amount of phase function moments, varies wildly. Try to use as small number as possible (higher number means more RAM usage, min. 3 (for Rayleigh scattering)
    "maxcmu" : 16, # number of computational streams
    "maxumu" : 1, # number of polar computational angles
    "maxphi" : 1, # number of computational azimuthal angles
    "maxulv" : 1, # number of output optical depths (layers)
    "usrang" : True, # user-defined output polar angles?
    "usrtau" : True, # user-defined output optical depths?
    "ibcnd" : 0, # boundary conditions case
    "onlyfl" : False, # return only fluxes?
    "prnt" : [False, False, False, False, False], # terminal print flags for DISORT
    "planck" : True, # toggle internal blackbody radiation
    "lamber"  : True, # bottom boundary Lambertian reflectivity
    "deltamplus" : True, # Delta-M-plus approximation to phase function
    "do_pseudo_sphere" : False, # spherical correction to atmosphere 
    "utau" : [0] # user output optical depths, list
    
    ## Other
    "disort_precision" : "double", # sets DISORT precision in Fortran, can be "double" or "single"
    
    ## Scattering layers
    # each layer is defined by a set of parameters. Use as many layers as you want. Do not overlap layers in altitude!
    # for ancillary information, keep values which are not retrieved
    "scat_lyrs_inputs" : { 
        "Sulphuric_acid_1" : {
                ##################################
                #see https://doi.org/10.1175/1520-0469(1950)007<0054:VALWCI>2.0.CO;2
                #see https://acp.copernicus.org/articles/10/7197/2010/#:~:text=The%20real%20geometrical%20thickness%20of%20optically%20thin,diffusive%20cloud%20tops%20than%20at%20higher%20latitudes.
                # avg_wc = 0.2 # g m-3 avg stratus cloud.
                # thickness of the cloud is 1.5 km, so this gives a mass loading as
                # 0.2*1500 = 300 g m-2                
                "name" : "Sulphuric_acid_1", # layer identifier/name
                "low_spc" : spc_wvnmlo, # spectral claculation grid lower limit 
                "upp_spc" : spc_wvnmhi, # spectral calculation grid upper limit
                "res" : 1, # specral calculation grid resolution
                "spec_units" : spc_units, # spectral calculation grid units
                "rho" : 1670, # kg m-3, particle density (water here), can be number or one of permitted strings
                "n" : None, # total particle concentration [cm-3]
                # see https://adele.faculty.ucdavis.edu/research/projects/cloud-dsd/#:~:text=Clouds%20are%20composed%20of%20countless,of%20processes%20related%20to%20clouds.
                "r" : 0.4,  # mean particle radius [um]
                "s" : 1.5, # spread, for the lognormal distribution
                "s_a_den" : None, # surf. area density, will be calculated in size dist.
                "v_den" : None, # volume density, will be calculated in size dist.
                "dist_type" : "log_normal", # choose size distribution type
                "comp" : "sulphuric acid", # define particle type (by composition)
        #                "comp" : "ash", # define particle type (by composition)
                "center_alt" : None, # avg particle layer altitude (center of layer altitude), [km]
                "thick" : None, # particle layer thickness, [km]
                "alt_upp" : 17, # particle layer upper boundary altitude, [km]
                "alt_low" : 16.1, # particle layer lower boundary altitude, [km],
                "radii" : 200, # number of radii in size distribution quadrature
                "eta" : 1e-6, # value n(r) at which the size distribution upper and lower limits are set
                "phase_quad_N" : 181, # number of quadrature points for the phase function (no. of angles)
                "phase_quad_type" : "L", # type of phase function quadrature
                "radii_quad_type" : "T", # type of radii size distribution quadrature
                "leg_coeffs" : True, # toggle legendre expansion coefficients for the phase function
                "leg_coeffs_type" : "normalised", # type of Legendre coeffs, normalised or regular
                "multiprocess" : False, # type of Legendre polynomial expansion coefficients     
            },
            "Ash_1" : {
                "name" : "Ash_1", # layer identifier/name
                "low_spc" : spc_wvnmlo, # spectral claculation grid lower limit 
                "upp_spc" : spc_wvnmhi, # spectral calculation grid upper limit
                "res" : 1, # specral calculation grid resolution
                "spec_units" : spc_units, # spectral calculation grid units
                "rho" : 2300, # particle density, can be number or one of permitted strings
                "n" : None, # total particle concentration [cm-3]
                "s" : 1.5, # spread, for the lognormal distribution
                "s_a_den" : None, # surf. area density, will be calculated in size dist.
                "v_den" : None, # volume density, will be calculated in size dist.
                "dist_type" : "log_normal", # choose size distribution type
                "comp" : "tongariro-ash_Reed.ri", # define particle type (by composition)
                "center_alt" : None, # avg particle layer altitude (center of layer altitude), [km]
                "thick" : None, # particle layer thickness, [km]
                "alt_upp" : 16, # particle layer upper boundary altitude, [km]
                "alt_low" : 15, # particle layer lower boundary altitude, [km],
                "radii" : 181, # number of radii in size distribution quadrature
                "eta" : 1e-6, # value n(r) at which the size distribution upper and lower limits are set
                "phase_quad_N" : 200, # number of quadrature points for the phase function (no. of angles)
                "phase_quad_type" : "L", # type of phase function quadrature
                "radii_quad_type" : "T", # type of radii size distribution quadrature
                "leg_coeffs" : True, # toggle legendre expansion coefficients for the phase function
                "leg_coeffs_type" : "normalised", # type of Legendre coeffs, normalised or regular
                "multiprocess" : False, # type of Legendre polynomial expansion coefficients     
            },
            "Water_cloud_1" : {
                "name" : "Water_cloud_1", # layer identifier/name
                "low_spc" : spc_wvnmlo, # spectral claculation grid lower limit 
                "upp_spc" : spc_wvnmhi, # spectral calculation grid upper limit
                "res" : 1, # specral calculation grid resolution
                "spec_units" : spc_units, # spectral calculation grid units
                "rho" : 997, # particle density, can be number or one of permitted strings, water density 997 kg m-3
                "n" : None, # total particle concentration [cm-3]
                "s" : 1.5, # spread, for the lognormal distribution
                "s_a_den" : None, # surf. area density, will be calculated in size dist.
                "v_den" : None, # volume density, will be calculated in size dist.
                "dist_type" : "log_normal", # choose size distribution type
                "comp" : "H2O_263K_Rowe_2020.ri", # define particle type (by composition)
                "center_alt" : None, # avg particle layer altitude (center of layer altitude), [km]
                "thick" : None, # particle layer thickness, [km]
                "radii" : 181, # number of radii in size distribution quadrature
                "eta" : 1e-6, # value n(r) at which the size distribution upper and lower limits are set
                "phase_quad_N" : 200, # number of quadrature points for the phase function (no. of angles)
                "phase_quad_type" : "L", # type of phase function quadrature
                "radii_quad_type" : "T", # type of radii size distribution quadrature
                "leg_coeffs" : True, # toggle legendre expansion coefficients for the phase function
                "leg_coeffs_type" : "normalised", # type of Legendre coeffs, normalised or regular
                "multiprocess" : False, # type of Legendre polynomial expansion coefficients     
            },
        },
    "plot_profiles" : True # bool, create a plot comparison or not
    "base_plots" : True # bool, toggle basic plotting outputs
    "savetxt" : True # bool, if True, save output spectrum to txt file
    
}


## Prepare profiles of species to retrieve (now H2O, SO2):
rfm_prf = rfm_functions.read_atm_file(b["atm"])

abs_path = os.getcwd() # current folder absolute path

## Grids
# note that calculation output wavenumber gridd are separate objects, i.e. you can run a high resolution calculation and then interpolate to lower resolution (e.g. to match IASI gridpoints)
fin_wvnmlo = 645,                                    # lower wavenumber for the final spectrum output, cm-1, int/float
fin_wvnmhi = 1500,                                   # upper wavenumber for the final spectrum output, cm-1, int/float
fin_res = 0.25,                                      # resolution of the final spectrum output, cm-1

spc_res = 0.05,                                      # main calculation spectral resolution
spc_wvnmlo = fin_wvnmlo - 3,                         # main calculation upper limit, if convolving with ils, must be in cm-1 and at least 2 cm-1 greater than fin_wvnmhi
spc_wvnmhi = fin_wvnmhi + 3,                         # main calculation upper limit, if convolving with ils, must be in cm-1 and at least 2 cm-1 greater than fin_wvnmhi
spc_units = "cm-1",                                  # main calculation units, cm-1, um or nm

### state vector, ordered, contains:

x = {
    ## Scattering layers
    # for state vector; keep dict items which are retrieved
    "scat_lyrs_inputs" : { 
        "Sulphuric_acid_1" : {
            "mass_loading" : 5, #column particle loading, [g m-2]
            "r" : 0.4,  # mean particle radius [um]    
        },
        "Ash_1" : {
            "mass_loading" : 0.21, #column particle loading, [g m-2]
            "r" : 0.1,  # mean particle radius [um]    
        },
        "Water_cloud_1" : {
            "mass_loading" : 100, #column particle loading, [g m-2] #cloud is 100 m thick, 1 g m-3 mass loading
            "r" : 15,  # mean particle radius [um]
            "alt_upp" : 3, # particle layer upper boundary altitude, [km]
            "alt_low" : 2, # particle layer lower boundary altitude, [km], # cloud 100 m thick 
        },
    },
    
    # atmospheric profiles of gases to retrieve
    "g_rtv" : {
        "SO2 [ppmv]" : rfm_prf["SO2 [ppmv]"],
        "H2O [ppmv]" : rfm_prf["H2O [ppmv]"]
        }
    
}



