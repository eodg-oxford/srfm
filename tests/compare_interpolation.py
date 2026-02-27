"""
This code calculates optical properties (extinction coefficient, single scatter albedo,
phase function and optionally its Legendre expansion at two different grids.
It then interpolates the calculation from the second grid to the first grid and
calculates the ensuing difference.
Optionally also creates a plot.
"""

import numpy as np
from srfm import *
import matplotlib.pyplot as plt
import warnings

###################################
# Assign some variables:
###################################
aria_fldr = "/network/group/aopp/eodg/RGG009_GRAINGER_EODGCOMN/ARIA/"  # where ARIA is

spec_res = 0.001  # model spectral resolution, [cm-1]
low_wvn = 750  # model start wavenumber (lower), [cm-1]
upp_wvn = 1250  # model end wavenumber (upper), [cm-1]

multiprocess = False  # don't change here

#######################################
# prepare particle scattering properties
#######################################

# define particle properties
particle_loading = 3.23  # column particle loading, [g m-2]
# n = 1e4 # total particle concentration [cm-3]
r = 2  # mean particle radius [um]
s = 1.5  # spread, for the lognormal distribution
s_a_den = None  # surf. area density, will be calculated in size dist.
v_den = None  # volume density, will be calculated in size dist.
dist_type = "log_normal"  # choose size distribution type
comp = "ash"  # define particle type (by composition)
p_lyr_a_avg = 3.5  # avg particle layer altitude (center of layer altitude), [km]

p_lyr_thick = 1  # particle layer thickness, [km]
n = utilities.number_conc_from_particle_loading(
    l=particle_loading, rho="glass", thick=p_lyr_thick, r=r
)

# create particle size distribution
sd = size_distribution.create_distribution(
    dist_type=dist_type,
    n=n,
    r=r,
    s=s,
    surface_area_density=s_a_den,
    volume_density=v_den,
)

# set-up the scattering calculation
radii = 200  # number of radii in size distribution quadrature
eta = 1e-6  # value n(r) at which the size distribution upper and lower limits are set
phase_quad_N = 181  # number of quadrature points for the phase function (no. of angles)
phase_quad_type = "L"  # type of phase function quadrature
radii_quad_type = "T"  # type of radii size distribution quadrature
leg_coeffs = True  # toggle legendre expansion coefficients for the phase function
leg_coeffs_type = "normalised"  # type of Legendre polynomial expansion coefficients

RFM_wvnm = np.linspace(
    low_wvn, upp_wvn, int((upp_wvn - low_wvn) / spec_res + 1), endpoint=True
)  # mirrors the RFM wavenumber grid
wvls = (1 / RFM_wvnm) * 1e4  # convert RFM wavenumber to wavelength in [um]

########################################################################################
# set up calculation of optical properties at grid same to the spectral grid (full res.)
########################################################################################

# calculate optical properties either standard or parallelized, where
#   op_dict["beta_ext"] - extinction coefficient
#   op_dict["ss_alb"] - single scatter albedo
#   op_dict["phase_function"] - phase function
#   op_dict["legendre_coefficient"] - coefficients of the Legendre expansion of the phase function, optional
full_res_op_dict = optical_properties.ewp_hs(
    wavelengths=wvls,
    composition=comp,
    distribution=sd,
    legendre_coefficients_flag=leg_coeffs,
    legendre_coefficients_type=leg_coeffs_type,
    radii=radii,
    eta=eta,
    phase_quad_N=phase_quad_N,
    phase_quad_type=phase_quad_type,
    radii_quad_type=radii_quad_type,
    aria=aria_fldr,
)


########################################################################################
# set up calculation of optical properties at lower resolution (low res.)
########################################################################################

ewp_res = 1  # resolution of the grid to calculate optical properties at
ewp_wvnm = np.linspace(
    low_wvn, upp_wvn, int((upp_wvn - low_wvn) / ewp_res + 1), endpoint=True
)  # grid for the optical properties
ewp_grid = (1 / ewp_wvnm) * 1e4  # convert ewo_wvnm wavenumbers to wavelength in [um]

# calculate optical properties either standard or parallelized, where
#   op_dict["beta_ext"] - extinction coefficient
#   op_dict["ss_alb"] - single scatter albedo
#   op_dict["phase_function"] - phase function
#   op_dict["legendre_coefficient"] - coefficients of the Legendre expansion of the phase function, optional
op_dict = optical_properties.ewp_hs(
    wavelengths=ewp_grid,
    composition=comp,
    distribution=sd,
    legendre_coefficients_flag=leg_coeffs,
    legendre_coefficients_type=leg_coeffs_type,
    radii=radii,
    eta=eta,
    phase_quad_N=phase_quad_N,
    phase_quad_type=phase_quad_type,
    radii_quad_type=radii_quad_type,
    aria=aria_fldr,
)

# regrid op_dict from the ewp_grid to wvls
# input("waiting for input:")
# old_dict = op_dict.copy()
op_dict, d_dict = optical_properties.regrid(
    op_dict, wvls, track_diff=True, diff_type="pct"
)

########################################################################################
# calculate difference between the two ewp calculations
########################################################################################
plt.ion()
diff_dict, fig = optical_properties.calc_op_diff(full_res_op_dict, op_dict, plot=True)
plt.show()
