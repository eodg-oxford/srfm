"""
Name: optical_properties
Author: Don Grainger
Contributors: Antonin Knizek
Date: 24 January 2025
Purpose: Used to calculate particle optical properties, such as extinction coefficient,
         scattering coefficient, scattering phase function and Legendre expansion 
         coefficients for the phase function.
"""

import numpy as np
from . import quadrature as quad
from . import ARIA_module as ARIA  # Import the RI class from ri_module
from . import mie_module  # Assuming mie_module is the compiled Fortran module
from . import utilities as utils


def legendre_polynomial_expansion(inp, qv, qw, phase):
    """
    Expands a function as a Legendre series. The function is assumed to have been
    evaluated at Legendre quadrature points. The evaluation stops when the number
    of terms exceeds the number of quadrature points or when the absolute value of
    the Legendre coefficient is less than 1e-9. Uses the Bonnet's recusion formula
    to generate the Legendre polynomials.

    Parameters:
        inp (int): Number of points at which the phase function is evaluated.
        qv (array-like): Legendre points (points at which the phase function is
                         evaluated, quadrature points).
        qw (array-like): Legendre weights, same shape as qv.
        phase (array-like): Input (phase) function.

    Returns:
        tuple: lc (Legendre coefficients), inlc (Number of coefficients used).
    """
    # Catch insensible inputs
    Imaxnp = 10000  # maximum allowed Legendre points
    if inp > Imaxnp:
        raise ValueError(
            "Error in legendre_polynomial_expansion: Too many quadrature points"
        )

    # Initialize array of Legendre coefficients
    lc = np.zeros(inp, dtype=np.float64)

    # Compute the first two Legendre coefficients
    lc[0] = np.sum(phase * qw) / 2.0
    lc[1] = 3.0 * np.sum(phase * qv * qw) / 2.0

    # Calculate the first two Legendre polynomials, P_0 and P_1
    lpnm2 = np.ones(inp, dtype=np.float64)  # P_0 = 1
    lpnm1 = qv  # P_1 = qv

    # Iterate to compute higher-order coefficients and polynomials
    # Higher-order Legendre polynomials (lpn) are calculated using the
    # Bonnet's recursion formula. n is the order of the respective polynomial
    n = 2
    while n < inp:
        # Compute the nth Legendre polynomial
        lpn = ((2 * n - 1) / n) * qv * lpnm1 - ((n - 1) / n) * lpnm2
        lpnm2 = lpnm1
        lpnm1 = lpn

        # Compute the nth Legendre coefficient
        lc[n] = (2 * n + 1) * np.sum(phase * lpn * qw) / 2.0

        # Stop if the coefficient is below the threshold
        if abs(lc[n]) < 1e-9:
            break

        n += 1

    # Number of coefficients used
    inlc = n - 1
    #    print(f'lc = {lc}')
    return lc, inlc


def normalised_legendre_polynomial_expansion(inp, qv, qw, phase):
    """
    Expands a function as a Legendre series. The function is assumed to have been
    evaluated at Legendre quadrature points. The evaluation stops when the number
    of terms exceeds the number of quadrature points or when the absolute value of
    the Legendre coefficient is less than 1e-9. Uses the Bonnet's recusion formula
    to generate the Legendre polynomials. Uses normalised Legendre polynomial
    coefficients (as in Wiscombe 1977
    https://doi.org/10.1175/1520-0469(1977)034<1408:TDMRYA>2.0.CO;2)

    Parameters:
        inp (int): Number of points at which the phase function is evaluated.
        qv (array-like): Legendre points (points at which the phase function is
                         evaluated, quadrature points).
        qw (array-like): Legendre weights, same shape as qv.
        phase (array-like): Input (phase) function.

    Returns:
        lc - normalised Legendre coefficients
        inlc - Number of coefficients used
    """
    # Catch insensible inputs
    Imaxnp = 10000  # maximum allowed Legendre points
    if inp > Imaxnp:
        raise ValueError(
            """"Error in normalised legendre_polynomial_expansion: 
            Too many quadrature points"""
        )

    # Initialize array of Legendre coefficients
    lc = np.zeros(inp, dtype=np.float64)

    # Compute the first two Legendre coefficients
    lc[0] = np.sum(phase * qw) / 2.0
    lc[1] = np.sum(phase * qv * qw) / 2.0

    # Calculate the first two Legendre polynomials, P_0 and P_1
    lpnm2 = np.ones(inp, dtype=np.float64)  # P_0 = 1
    lpnm1 = qv  # P_1 = qv

    # Iterate to compute higher-order coefficients and polynomials
    # Higher-order Legendre polynomials (lpn) are calculated using the
    # Bonnet's recursion formula. n is the order of the respective polynomial
    n = 2
    while n < inp:
        # Compute the nth Legendre polynomial
        lpn = ((2 * n - 1) / n) * qv * lpnm1 - ((n - 1) / n) * lpnm2
        lpnm2 = lpnm1
        lpnm1 = lpn

        # Compute the nth Legendre coefficient
        lc[n] = np.sum(phase * lpn * qw) / 2.0

        # Stop if the coefficient is below the threshold
        if abs(lc[n]) < 1e-9:
            break

        n += 1

    # Number of coefficients used
    inlc = n - 1
    #    print(f'lc = {lc}')
    return lc, inlc


# Determine the optical properties of a distribution of particles as a function of
# composition and structure
def ewp_hs(
    wavelength,
    composition,
    distribution,
    refractive_index=None,
    angle=None,
    legendre_coefficients_flag=False,
    legendre_coefficients_type="normalised",
    radii=200,
    eta=1e-6,
    phase_quad_N=181,
    phase_quad_type="L",
    radii_quad_type="T",
    aria=None
):
    """
    return the extinction, singlescatter albedo and phase function for a particle
    size of a given 'distribution' and 'composition'
    the phase function is given at the specified angles if given,
    else at 1 degree resolution
    if the legendre_coefficients_flag is True angle is ignored and Lobatto quadrature
    angles are used
    radii - number of radii of the particles in the distribution
    eta - cutoff number for the size distribution
    legendre_coefficients_type - normalised or regular, different functions are called
    for the expansion
    phase_quad_N - number of phase function moments (scattering angles), used when
    legendre expansion requested.
    phase_quad_type - when legendre expansion requested, this type of quadrature is used
    valid input are "G" (Gaussian), "R" (Radau), "L" (Lobatto)
    refractive_index - if compositon is "ri", then refractive indices are required
    from the user.
    szd_quad_type = quadrature type for the particle size distribution, currently "T"
    """
    # get size of the wavelength array (amount of elements)
    wavelengths = wavelength.size

    # calculates cosines of the angles to return the phase function at.
    # if Legendre expansion of the phase function is also requested
    # (i.e. if legendre_coefficients_flag == True)
    # the user angles are ignored. The quadrature points are calculated
    # by the Quadrature module.
    # There are several possible quadratures, here Lobatto is used by default and
    # the angles are from 0 to 180 in steps of 1 degree.
    # elif the Legendre expansion of the phase function is not required and angles are
    # given, simply calculates the cosines,
    # lastly, elif Legendre expansion is not needed and angles are not given, one angle
    # is used.
    # If angles are given, calculates the phase function at those angles otherwise use
    # default 0 to 180 in steps of 1 degree (and flip)
    if legendre_coefficients_flag:
        cos_angle_value, cos_angle_weight = quad.quadrature101(
            phase_quad_type, phase_quad_N
        )

        # raise warning that user angles are being ignored
        if angle is not None:
            print(
                """Warning in ewp_hs: Specified angles ignored as Legendre
                 coefficients requested"""
            )
    else:
        if angle is not None:
            cos_angle_value = np.cos(angle * np.pi / 180)
        #            cos_angle_weight = np.zeros(cos_angle_value.size)
        else:
            cos_angle_value, cos_angle_weight = (0, 1)  # 1 angle

    cos_angle_value = np.array(
        cos_angle_value, dtype=np.float64
    )  # ensure correct type for Mie call
    angles = cos_angle_value.size

    # calculates the lower and upper particle radius from the distribution at points
    # where n(r) = eta*r_mode
    radius_lower_bound = np.exp(
        np.log(distribution.r)
        - np.log(distribution.s) ** 2
        - np.sqrt(-2 * np.log(eta) * np.log(distribution.s) ** 2)
    )
    radius_upper_bound = np.exp(
        np.log(distribution.r)
        - np.log(distribution.s) ** 2
        + np.sqrt(-2 * np.log(eta) * np.log(distribution.s) ** 2)
    )

    # determine quadrature points for the distribution, 
    # radii = number of required quadrature points
    radius, weight_temp = quad.quadrature(
        radii_quad_type, radii, radius_lower_bound, radius_upper_bound
    )
    radius_weight = weight_temp * distribution.value(radius)
    rad_wt_sum = np.sum(radius_weight)

    # load refractive indices (expected in the form (n - ik)
    if composition == "ri":
        if refractive_index is None:
            print("Error: refractive index not defned")
        else:
            ri_n = np.array(np.full(wavelengths, np.real(refractive_index)))
            ri_k = np.array(np.full(wavelengths, np.imag(refractive_index)))
    else:
        ri_object = ARIA.RI()
        ri_n, ri_k = ri_object.load_refractive_indices(
            composition, aria=aria, wave=wavelength, mode="wavelength"
        )
        ri_k = -ri_k

    # check refractive index before proceeding
    if np.any(ri_k > 0):
        raise ValueError(
            """Refractive index (k) contains invalid values. 
               Ensure all values are zero or negative."""
        )

    # Initialize arrays with the same size as 'wavelengths'
    beta_ext = np.zeros(wavelengths)  # 1D array for extinction
    beta_sca = np.zeros(wavelengths)  # 1D array for scattering
    phase_function = np.zeros((wavelengths, angles))  # 2D array for phase function

    # Initialize arrays with the same size as 'radii'
    Q_ext_value = np.zeros(radii)  # 1D array for extinction per radius
    Q_sca_value = np.zeros(radii)  # 1D array for scattering per radius
    phase_function_value = np.zeros(
        (radii, angles)
    )  # 2D array for phase function per radius

    # initialize single element arrays (values) for Q_sca and Q_ext
    Q_ext = np.zeros(1, dtype=np.float64)  # Single output scalar wrapped in an array
    Q_sca = np.zeros(1, dtype=np.float64)  # Single output scalar wrapped in an array

    # initialize arrays for the phase function and error
    phase_function_particle = np.zeros(angles, dtype=np.float64)  # Complex output array
    Error = np.array([0], dtype=np.int32)  # Error code as an integer array

    # initialize array for legendre coefficients (take care not to overlow memory)
    if legendre_coefficients_flag:
        legendre_coefficient = utils.memory_safe_np_zeros_2d(
            constraints = [wavelengths], max_sec_dim=20000
        )
        
    
    max_lc = 0 # tracks length of legendre expansion in the main computational loop
    # (number of Legendre polynomial coefficients used. Used to truncate the final 
    # array.
    
    # main computational loop, first iterate over wavelengths
    for i in range(wavelengths):
        print(f"Calculating optical properties for {wavelength[i]:.4f} um.")
        # make refractive index a complex number
        refractive_index = complex(ri_n[i], ri_k[i])

        # loop over radii
        for j in range(radii):
            size_parameter = 2 * np.pi * radius[j] / wavelength[i]

            # Perform Mie calculations
            Q_ext, Q_sca, phase_function_particle, error = mie_module.mie_ewp(
                size_parameter, refractive_index, cos_angle_value
            )

            # Store results for this radius
            Q_ext_value[j] = Q_ext
            Q_sca_value[j] = Q_sca
            phase_function_value[j, :] = (
                phase_function_particle  # Ensure correct slicing
            )

        # Weighted sums for extinction and scattering
        beta_ext[i] = (
            np.sum(radius_weight * np.pi * radius**2 * Q_ext_value) * 1e-6
        )  # 1/m
        beta_sca[i] = (
            np.sum(radius_weight * np.pi * radius**2 * Q_sca_value) * 1e-6
        )  # 1/m

        # Weighted sums for phase function
        for k in range(angles):
            phase_function[i, k] = (
                np.sum(radius_weight * phase_function_value[:, k]) / rad_wt_sum
            )

        # If requested expand phase functions as Legendre coefficients
        if legendre_coefficients_flag:
            if legendre_coefficients_type == "regular":
                legendre_coefficient_temp, legendre_coefficient_number = (
                    legendre_polynomial_expansion(
                        angles, cos_angle_value, cos_angle_weight, phase_function[i, :]
                    )
                )
            elif legendre_coefficients_type == "normalised":
                legendre_coefficient_temp, legendre_coefficient_number = (
                    normalised_legendre_polynomial_expansion(
                        angles, cos_angle_value, cos_angle_weight, phase_function[i, :]
                    )
                )
            else:
                raise ValueError(
                    """Invalid Legendre coefficient type. Accepted values 
                are 'normalised' and 'regular'."""
                )

            legendre_coefficient[i, 0:legendre_coefficient_number] = (
                legendre_coefficient_temp[0:legendre_coefficient_number]
            )
            
            if legendre_coefficient_number > max_lc:
                max_lc = legendre_coefficient_number
    
    # truncate zeros from legendre coefficient
    legendre_coefficient = legendre_coefficient[:,:max_lc]
    
    # Return the computed results
    if legendre_coefficients_flag:
        return (
            beta_ext,
            beta_sca / beta_ext,
            phase_function,
            legendre_coefficient,
            cos_angle_value,
            cos_angle_weight,
        )
    else:
        return beta_ext, beta_sca / beta_ext, phase_function


def phase_from_legendre(inlc, lc, inp, qv):
    """Recomputes the phase function from Legendre polynomial expansion.
    Inputs:
        inlc - number of Legendre coefficients
        lc - Legendre expansion coefficients
        inp - number of quadrature value, == number of angles
        qv - cos(angles)
    Returns:
        phase - recomputed phase function
    """

    # initialize arrays
    imaxnp = 20000
    phase = np.zeros(inp)
    lpn = np.zeros(imaxnp)
    lpnm1 = np.zeros(imaxnp)
    lpnm2 = np.zeros(imaxnp)

    # loop over angles
    for i in range(0, inp):
        # add zeroth polynomial, P_0(x)
        lpnm2[i] = 1
        phase[i] = lc[0]

        # add first polynomial, P_1(x)
        lpnm1[i] = qv[i]
        phase[i] = phase[i] + lc[1] * lpnm1[i]

        for n in range(2, inlc):
            lpn[i] = ((2 * n - 1) * qv[i] * lpnm1[i] - (n - 1) * lpnm2[i]) / n
            phase[i] = phase[i] + lc[n] * lpn[i]
            lpnm2[i] = lpnm1[i]
            lpnm1[i] = lpn[i]

    return phase


def phase_from_normalised_legendre(inlc, lc, inp, qv):
    """Recomputes the phase function from Legendre polynomial expansion.
    Inputs:
        inlc - number of Legendre coefficients
        lc - normalised Legendre expansion coefficients

        inp - number of quadrature value, == number of angles
        qv - cos(angles)
    Returns:
        phase - recomputed phase function
    """

    # initialize arrays
    imaxnp = 20000
    phase = np.zeros(inp)
    lpn = np.zeros(imaxnp)
    lpnm1 = np.zeros(imaxnp)
    lpnm2 = np.zeros(imaxnp)

    # loop over angles
    for i in range(0, inp):
        # add zeroth polynomial, P_0(x)
        lpnm2[i] = 1
        phase[i] = lc[0]

        # add first polynomial, P_1(x)
        lpnm1[i] = qv[i]
        phase[i] = phase[i] + 3 * lc[1] * lpnm1[i]

        for n in range(2, inlc):
            lpn[i] = ((2 * n - 1) * qv[i] * lpnm1[i] - (n - 1) * lpnm2[i]) / n
            phase[i] = phase[i] + (2 * n + 1) * lc[n] * lpn[i]
            lpnm2[i] = lpnm1[i]
            lpnm1[i] = lpn[i]

    return phase

