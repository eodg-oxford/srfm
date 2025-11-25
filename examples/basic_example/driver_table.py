import os
#import sys
#from pathlib import Path

#DEFAULT_SRFM_DIR = Path("/home/k/knizek/Documents/srfm_main/src")
#SRFM_HOME = Path(os.environ.get("SRFM_HOME", DEFAULT_SRFM_DIR))
#if SRFM_HOME.is_dir() and str(SRFM_HOME) not in sys.path:
#    sys.path.insert(0, str(SRFM_HOME))

from srfm import *

# Define common grids and helper constants first so they can be reused below
ABS_PATH = os.getcwd()

FIN_WVNMLO = 645.0
FIN_WVNMHI = 1500.0
FIN_RES = 0.25

SPC_RES = 0.05
SPC_WVNMLO = FIN_WVNMLO - 3.0
SPC_WVNMHI = FIN_WVNMHI + 3.0
SPC_UNITS = "cm-1"

inputs = {
    ## Retrieval configuration
    "fwd_model" : "SRFM", # which forward model to use, currently only SRFM
    "instrument" : "IASI", # which instrument are data from, currently only IASI
    
    ## Measurement configuration
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
    
    # inputs that mirror the standard RFM driver table
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
    "fisot": 0.0,
    "albedo": 0.0,
    "temis": 1.0,
    "earth_radius": 6371.0,
    "nmom": 3,
    "maxcmu": 16,
    "maxumu": 1,
    "maxphi": 1,
    "maxulv": 1,
    "usrang": True,
    "usrtau": True,
    "ibcnd": 0,
    "onlyfl": False,
    "prnt": [False, False, False, False, False],
    "planck": True,
    "lamber": True,
    "deltamplus": True,
    "do_pseudo_sphere": False,
    "utau": [0.0],
    "disort_precision": "double",
    "header": "NO HEADER", # header for terminal printing, "NO HEADER" == supressed.

    ## Scattering configuration
    "scat_lyrs_inputs": {
        "Sulphuric_acid_1": {
            "name": "Sulphuric_acid_1",
            "low_spc": SPC_WVNMLO,
            "upp_spc": SPC_WVNMHI,
            "res": 3.0,
            "spec_units": SPC_UNITS,
            "rho": 1670.0,
            "n": None,
            "s": 1.75,
            "s_a_den": None,
            "v_den": None,
            "dist_type": "log_normal",
            "comp": "sulphuric acid",
            "center_alt": 14.0,
            "thick": 1.5,
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
            "mass_loading": 0.2,
            "r": 0.4
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
