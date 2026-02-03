import os
from srfm import *

# Define common grids and helper constants first so they can be reused below
ABS_PATH = os.getcwd()

# Final grid to output results at, units cm-1
FIN_WVNMLO = 920.0 # min
FIN_WVNMHI = 925.0 # max
FIN_RES = 0.25 # resolution

# Computational grid
SPC_RES = 0.001 # resolution
SPC_WVNMLO = FIN_WVNMLO - 3.0 # min
SPC_WVNMHI = FIN_WVNMHI + 3.0 # max
SPC_UNITS = "cm-1" # units (cm-1, nm, um)

inputs = {
    ## Retrieval configuration (optional)
    "fwd_model" : "SRFM", # which forward model to use, currently only SRFM
    "instrument" : "IASI", # which instrument are data from, currently only IASI
    
    ## Measurement configuration (optional)
    "iasi_spc_fldr": "/network/group/aopp/eodg/RGG008_GRAINGER_IASIVOLC/knizek/raikoke_ssaccs/spectra_w_iasi_L2",
    "iasi_fl": "hlat_20190701_A.pkl",
    "px": 8,

    ## General    
    "plot_profiles": True, # (optional) is using srfm.iasi_main
    "base_plots": True, # if srfm.main, then plots output spectrum, if srfm.iasi_main, plots output spectrum and comparisons
    "out_mode": "netcdf", # file format to save to, netcdf, txt or None
    "show_plots": True, # show plots on screen
    "results_fldr": os.path.join(ABS_PATH, "results"),
    "rad": False, # True - save radiances
    "bbt": True, # True - save brightness temperatures
    "rad_out_fname": None, # manually set radiances output filename
    "bbt_out_fname": None, # manually set bbt output filename
    
    # insturment line shape (ILS)
    "convolve_iasi": False, # if True, convolves spectrum with IASI ILS
    "iasi_ils": os.path.join(ABS_PATH, "iasi.ils"), # (optional) path to ILS file

    ## Grids
    "fin_wvnmlo": FIN_WVNMLO,
    "fin_wvnmhi": FIN_WVNMHI,
    "fin_res": FIN_RES,
    "spc_res": SPC_RES,
    "spc_wvnmlo": SPC_WVNMLO,
    "spc_wvnmhi": SPC_WVNMHI,
    "spc_units": SPC_UNITS,
    
    ## RFM configuration:
    
    # RFM global config
    "rfm_config" : {
        "output_mode": "capture",            # "capture" or "files"
        "driver_path": None,               # e.g. Path("drvs/example.drv")
        "generate_driver": False,          # Set True to overwrite driver_path contents
        "verbose": False,
        "capture_files_content": False # True - store output files content in result
        # the logic is to save memory if processed outputs, such as dataframes are requested
        # and not capture the raw files as well.
    },
    
    # inputs that mirror the standard RFM driver table https://eodg.atm.ox.ac.uk/RFM/index.html
    "driver_inputs":
        dict(
        # mandatory sections
        header="SRFM_standard_run",
        flags=("OPT", "NAD", "SFC", "PRF", "LEV", "DBL", "CHI", "MIX"),
        spectral=[rfm_helper.SpectralRange(SPC_WVNMLO, SPC_WVNMHI, SPC_RES)],
        gases=("N2", "O2", "CO2", "O3", "H2O", "CH4", "N2O", "HNO3", "CO", "NO2", 
            "N2O5", "ClO", "HOCl", "ClONO2", "NO", "HNO4", "HCN", "NH3", "F11", "F12",
            "F14", "F22", "CCl4", "COF2", "H2O2", "C2H2", "C2H6", "OCS", "SO2"),
        atmosphere=(
            os.path.join(ABS_PATH,"hgt_std.atm"),
            os.path.join(ABS_PATH,"day.atm"), # make sure that the actual profiles are the second item in this tuple
            # TODO change the reliance on positions
        ),

        hit=("/network/aopp/matin/eodg/crun/eodg/rfm/rfm_files/bin/hitran_2012.bin",),
        xsc=("/network/aopp/matin/eodg/crun/eodg/rfm/rfm_files/xsc/*.xsc",),
    ),

    ## DISORT configuration
    "fisot": 0.0, # isotropic illumination at the top of the atmosphere
    "albedo": 0.0, # bottom boundary albedo
    "temis": 1.0, # top boundary emissivity
    "earth_radius": 6371.0, # Earth radius (km0
    "nmom": 3, # number of phase function moments
    "maxcmu": 16, # number of computational streams
    "maxumu": 1, # number of user output polar angles
    "maxphi": 1, # number of user azimuth angles
    "maxulv": 1, # number of user optical depths
    "usrang": True, # return output at user angles?
    "usrtau": True, # return output at user optical depths?
    "ibcnd": 0, # boundary conditions
    "onlyfl": False, # return only fluxes?
    "prnt": [False, False, False, False, False], # what gets prints to terminal
    "planck": True, # include internal Planck function?
    "lamber": True, # Lambertian reflector surface
    "deltamplus": True, # Delta-M+ approximation to phase functions
    "do_pseudo_sphere": False, # bent surface?
    "utau": [0.0], # user optical depths for output
    "disort_precision": "double", # Fortran precision
    "header": "NO HEADER", # header for terminal printing, "NO HEADER" == supressed.

    ## Scattering configuration
    # scattering layers are named and are as keys in this dict, refer to docs for 
    # specific parameters
    "scat_lyrs_inputs": {
        "Sulphuric_acid_1": {
            "name": "Sulphuric_acid_1",
            "low_spc": SPC_WVNMLO,
            "upp_spc": SPC_WVNMHI,
            "res": 3.0,
            "spec_units": SPC_UNITS,
            "rho": 1670.0, # scatterer density
            "n": None, # number concentration
            "s": 1.75, # size distribution spread
            "s_a_den": None, # surface area density
            "v_den": None, # volume density
            "dist_type": "log_normal", # size distribution type
            "comp": "sulphuric acid", # refractive index
            "center_alt": 14.0, # scattering layer center altitude
            "thick": 1.5, # scattering layer thickness
            "alt_upp": None, # scattering layer upper boundary
            "alt_low": None, # scattering layer lower boundary
            "radii": 181, # number of particle radii in particle size distribution
            "eta": 1e-6, # size distribution cut-off
            "phase_quad_N": 200, # number of quadrature points in the phase function
            "phase_quad_type": "L", # phase function quadrature type
            "radii_quad_type": "T", # size distribution quadrature type
            "leg_coeffs": True, # return Legendre expansion coefficients
            "leg_coeffs_type": "normalised", # what type of coefficients
            "multiprocess": False, # attempt to parallelize calculations, EXPERIMENTAL
            "mass_loading": 0.2, # scatterer mass loading (g m-2)
            "r": 0.4 # particle mean radius
        },
        "Ash_1": {
            "name": "Ash_1",
            "low_spc": SPC_WVNMLO,
            "upp_spc": SPC_WVNMHI,
            "res": 1.0,
            "spec_units": SPC_UNITS,
            "rho": 2300.0,
            "n": None,
            "s": 1.8,
            "s_a_den": None,
            "v_den": None,
            "dist_type": "log_normal",
            "comp": "ash",
            "center_alt": 9.0,
            "thick": 1.0,
            "alt_upp": None,
            "alt_low": None,
            "radii": 181,
            "eta": 1e-6,
            "phase_quad_N": 200,
            "phase_quad_type": "L",
            "radii_quad_type": "T",
            "leg_coeffs": True,
            "leg_coeffs_type": "normalised",
            "multiprocess": False,
            "mass_loading": 0.21,
            "r": 0.1,
        },
        "Water_cloud_1": {
            "name": "Water_cloud_1",
            "low_spc": SPC_WVNMLO,
            "upp_spc": SPC_WVNMHI,
            "res": 1.0,
            "spec_units": SPC_UNITS,
            "rho": 997.0,
            "n": None,
            "s": 1.5,
            "s_a_den": None,
            "v_den": None,
            "dist_type": "log_normal",
            "comp": "H2O_263K_Rowe_2020.ri",
            "center_alt": None,
            "thick": None,
            "radii": 181,
            "eta": 1e-6,
            "phase_quad_N": 200,
            "phase_quad_type": "L",
            "radii_quad_type": "T",
            "leg_coeffs": True,
            "leg_coeffs_type": "normalised",
            "multiprocess": False,
            "mass_loading": 100.0,
            "r": 15.0,
            "alt_upp": 3.0,
            "alt_low": 2.0,
        },
    },
    
    ## solar reflection
    "sun": True, # if True, include solar reflection
    "sza": 0, # solar zenith angle, 0-180, 0 for directly overhead, >90 for night (not included)
    "saa": 0, # solar azimuth angle, 0-360
    
    ## Angles
    "zen": 0, # satellite zenith angle, 0-180, 0 for directly overhead, >90 for night (not included)
    "azi": 0, # satellite azimuth angle, 0-360
}
