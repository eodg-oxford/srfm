"""Module for calculation of particle optical properties.

These are the extinction coefficient, scattering coefficient, scattering phase function
and Legendre expansion coefficients for the phase function.

- Name: optical_properties
- Parent package: srfm
- Author: Don Grainger
- Contributors: Antonin Knizek
- Date: 24 January 2025
"""

import numpy as np
from . import quadrature as quad
from . import ARIA_module as ARIA  # Import the RI class from ri_module
from . import mie_module  # Assuming mie_module is the compiled Fortran module
from . import utilities as utils
import multiprocessing
import warnings
from numba import complex128, jit
import time
import matplotlib.pyplot as plt


@jit(nopython=True, error_model="numpy", parallel=False, fastmath=True)
def legendre_polynomial_expansion(inp, qv, qw, phase):
    """Expands a function as a Legendre series.

    The function is assumed to have been
    evaluated at Legendre quadrature points. The evaluation stops when the number
    of terms exceeds the number of quadrature points or when the absolute value of
    the Legendre coefficient is less than 1e-9. Uses the Bonnet's recusion formula
    to generate the Legendre polynomials.

    Args:
        inp (int): Number of points at which the phase function is evaluated.
        qv (array-like): Legendre points (points at which the phase function is
                         evaluated, quadrature points).
        qw (array-like): Legendre weights, same shape as qv.
        phase (array-like): Input (phase) function.

    Returns:
        lc (array): Legendre coefficients.
        inlc (int): Number of coefficients used.

    Raises:
        ValueError: Raised when >10,000 quadrature points requested.

    """
    # Catch insensible inputs
    Imaxnp = 10000  # maximum allowed Legendre points
    if inp > Imaxnp:
        raise ValueError("""Error in legendre_polynomial_expansion: Too many quadrature
        points. Limit is 10,000.""")

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


@jit(nopython=True, error_model="numpy", parallel=False, fastmath=True)
def normalised_legendre_polynomial_expansion(inp, qv, qw, phase):
    """Expands a function as a Legendre series.

    The function is assumed to have been
    evaluated at Legendre quadrature points. The evaluation stops when the number
    of terms exceeds the number of quadrature points or when the absolute value of
    the Legendre coefficient is less than 1e-9. Uses the Bonnet's recusion formula
    to generate the Legendre polynomials. Uses normalised Legendre polynomial
    coefficients (as in `Wiscombe 1977`_).

    .. _Wiscombe 1977: https://doi.org/10.1175/1520-0469(1977)034<1408:TDMRYA>2.0.CO;2

    Args:
        inp (int): Number of points at which the phase function is evaluated.
        qv (array-like): Legendre points (points at which the phase function is
                         evaluated, quadrature points).
        qw (array-like): Legendre weights, same shape as qv.
        phase (array-like): Input (phase) function.

    Returns:
        lc (array): Legendre coefficients.
        inlc (int): Number of coefficients used.

    Raises:
        ValueError: Raised when >10,000 quadrature points requested.

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
@utils.show_runtime
def ewp_hs(
    wavelengths,
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
    return_dict=None,
    multiprocess=False,
):
    """Main function to calculate extinction, single scatter albedo and phase function.

    The calculation is performed for a set of particles of given size distribution and
    compositon.

    Args:
        wavelengths (array-like): wavelengths in :math:`\\mu`\ m
        composition (str): any of ash, sulphuric acid, ice or ARIA .ri filename
        distribution (obj): size distribution object from srfm.size_distribution module
        refractive_index (array-like): If compositon is "ri", then refractive indices
            are required from the user. Default is None.
        angle (array-like, float, int): Scattering angles, if None, then phase function
            returned at (0,180) degrees with one degree spacing. Can't be set when
            Legendre polynomial expansion of the phase function is requested (in the
            case quadrature angles are computed and used). Default is None.
        legendre_coefficients_flag (bool): If True, Legendre polynomial expansion of the
            phase function is performed and expansion coefficients are returned.
        radii (int, float): Number of radii of the particles in the distribution.
            Default is 200.
        eta (float): Cut-off number for the size distribution. Default is 1e-6.
        legendre_coefficients_type (str): Specifies the type of requested Legendre
            polynomial expansion coefficients. Can be *normalised* or *regular*. Default
            is *normalised*.
        phase_quad_N (int): Number of phase function moments (scattering angles), used
            when Legendre polynomial expansion requested. Default is 181.
        phase_quad_type (str): when Legendre polynomial expansion requested, this type
            of quadrature is used. Can be "G" (Gaussian), "R" (Radau) or "L" (Lobatto).
            Default is L.
        radii_quad_type (str): This type of quadrature is used to calculate the size
            distribution. Default is "T" (trapezium). Can be any quadrature accepted by
            `srfm.quadrature` module.
        return_dict (multiprocess.Manager().dict()): (optional) object to return values
            in when doing multiprocessing (parallel computations). Default is None.
        multiprocess (bool): If True, optical properties are calculated as a shell
        subprocess and the output is different (see Returns below).

    Returns:
        native_dict: Dictionary with calculated optical properties. Returned only if
            multiprocessing == False. Elif multiprocessing == True, returns an updated
            return_dict.

    Raises:
        ValueError: Raised when multiprocess has an invalid value.

    """

    # get size of the wavelengths array (amount of elements)
    wavelengths_size = np.int64(wavelengths.size)

    # calculate scattering angles quadrature
    cos_angle_value, cos_angle_weight, angles = get_quad(
        legendre_coefficients_flag=legendre_coefficients_flag,
        phase_quad_type=phase_quad_type,
        phase_quad_N=phase_quad_N,
        angle=angle,
    )

    # calculate size distribution quadrature
    radius, radius_weight, rad_wt_sum = get_radii(
        distribution=distribution, eta=eta, radii_quad_type=radii_quad_type, radii=radii
    )
    # load refractive indices
    refractive_index_arr = get_ri(
        composition=composition,
        refractive_index=refractive_index,
        wave=wavelengths,
        wave_size=wavelengths_size,
    )

    # calculate size_parameters
    wavelengths_arr = np.resize(wavelengths, (radius.shape[0], wavelengths_size))
    radius_arr = np.resize(radius, (wavelengths_size, radius.shape[0])).T
    size_parameters = 2 * np.pi * radius_arr / wavelengths_arr

    # Initialize arrays with the same size as 'wavelengths_size'
    beta_ext = np.zeros(wavelengths_size, dtype=np.float64)  # 1D array for extinction
    beta_sca = np.zeros(wavelengths_size, dtype=np.float64)  # 1D array for scattering
    phase_function = np.zeros(
        (wavelengths_size, angles), dtype=np.float64
    )  # 2D array for phase function

    # initialize array for legendre coefficients (take care not to overlow memory)
    if legendre_coefficients_flag:
        legendre_coefficient = utils.memory_safe_np_zeros_2d(
            constraints=[wavelengths_size], max_sec_dim=20000
        )

    beta_ext, beta_sca, p_f, l_c, max_lc = loop_mie_over_wavelengths(
        wavelengths_size,
        radii,
        size_parameters,
        refractive_index_arr,
        angles,
        cos_angle_value,
        cos_angle_weight,
        radius_weight,
        radius,
        rad_wt_sum,
        legendre_coefficients_flag,
        legendre_coefficients_type,
        beta_ext,
        beta_sca,
        phase_function,
        legendre_coefficient,
    )

    # truncate trailing zeros from legendre coefficient
    if legendre_coefficients_flag:
        l_c = l_c[:, :max_lc]

    # Return the computed results
    if multiprocess == True:
        native_return_dict = {}
        if legendre_coefficients_flag:
            native_return_dict["wavelengths"] = wavelengths
            native_return_dict["beta_ext"] = beta_ext
            native_return_dict["ssalb"] = beta_sca / beta_ext
            native_return_dict["phase_function"] = p_f
            native_return_dict["legendre_coefficient"] = l_c
        else:
            native_return_dict["wavelengths"] = wavelengths
            native_return_dict["beta_ext"] = beta_ext
            native_return_dict["ssalb"] = beta_sca / beta_ext
            native_return_dict["phase_function"] = p_f
        return_dict.update(native_return_dict)
        return

    elif multiprocess == False:
        native_return_dict = {}
        if legendre_coefficients_flag:
            native_return_dict["wavelengths"] = wavelengths
            native_return_dict["beta_ext"] = beta_ext
            native_return_dict["ssalb"] = beta_sca / beta_ext
            native_return_dict["phase_function"] = p_f
            native_return_dict["legendre_coefficient"] = l_c
        else:
            native_return_dict["wavelengths"] = wavelengths
            native_return_dict["beta_ext"] = beta_ext
            native_return_dict["ssalb"] = beta_sca / beta_ext
            native_return_dict["phase_function"] = p_f
        return native_return_dict

    else:
        raise ValueError("multiprocess must be True or False.")


def phase_from_legendre(inlc, lc, inp, qv):
    """Recomputes the phase function from Legendre polynomial expansion.

    Args:
        inlc (int): Number of Legendre coefficients.
        lc (array-like): Legendre expansion coefficients.
        inp (int): Number of quadrature value, == number of angles.
        qv (array-like): Cos(angles).

    Returns:
        phase (array): Recomputed phase function.

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

    Args:
        inlc (int): Number of Legendre coefficients.
        lc (array-like): Normalised Legendre expansion coefficients.
        inp (int): Number of quadrature value, == number of angles.
        qv (array-like): Cos(angles).

    Returns:
        phase (array): Recomputed phase function.

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


# @utils.show_runtime
#@jit(nopython=True, error_model="numpy", parallel=False, fastmath=True)
def mie_ewp(Dx, SCm, Dqv):
    """Mie calculation (python version).

    Args:
        Dx (float): Particle size parameter.
        SCm (complex): Particle refractive index.
        Dqv (array-like): Cosines of the scattering angles.

    Returns:
        Dqxt (float): Extinction coefficient.
        Dqsc (float): Scattering coefficient.
        pf (array): Phase function.

    Raises:
        ValueError: Raised when size parameter greater than 105000 or when
            Nmx > Itermax.

    """
    if isinstance(Dqv, list):
        Inp = len(Dqv)
    elif isinstance(Dqv, np.ndarray):
        if Dqv.ndim == 1:
            Inp = Dqv.shape[0]
        else:
            raise ValueError(f"""Dqv must be a 1D Array or list, but currently is
            a {Dqv.ndim}-D array.""")
    else:
        raise TypeError(f"""Dqv must be a list or 1-D np.ndarray, but currently is 
            {type(Dqv)}.""")
    
    Imaxx = 105000
    Itermax = 106000

    if Dx > Imaxx:
        raise ValueError("Dx is greater than 105000.")

    Ir = 1 / SCm  # np.complex128
    Y = Dx * SCm  # np.complex128

    NStop = (
        np.int64(max(2, (Dx + 4.05 * Dx ** (1 / 3) + 2)))
        if Dx < 4200
        else np.int64(Dx + 4 * Dx ** (1 / 3) + 2)
    )
    NmX = np.int64(max(NStop, abs(Y).real) + 15)

    if NmX > Itermax:
        raise ValueError(f"Nmx > Itermax.")

    D = np.zeros(NmX + 1, dtype=np.complex128)  # array, np.complex128

    for N in range(NmX - 1, 0, -1):
        A1 = (np.float64(N) + 1) / Y
        D[N] = A1 - 1 / (A1 + D[N + 1])

    Sp = np.zeros(Inp, dtype=np.complex128)
    Sm = np.zeros(Inp, dtype=np.complex128)
    Pi0 = np.zeros(Inp, dtype=np.float64)
    Pi1 = np.ones(Inp, dtype=np.float64)

    Psi0, Psi1 = np.cos(Dx), np.sin(Dx)  # np.float64
    Chi0, Chi1 = -Psi1, Psi0  # np.float64
    Xi0, Xi1 = np.complex128(Psi0 + 1j * Chi0), np.complex128(
        Psi1 + 1j * Chi1
    )  # np.complex128

    Dqxt = np.float64(0.0)  # np.float64
    Dqsc = np.float64(0.0)  # np.float64
    Tnp1 = np.int64(1)  # np.int64

    for N in range(1, NStop + 1):
        DN = np.float64(N)
        Tnp1 = Tnp1 + 2  # np.int64
        A2 = Tnp1 / (DN * (DN + 1.0))  # np.float64
        Turbo = (DN + 1.0) / DN  # np.float64
        Rnx = DN / Dx  # np.float64

        Psi = (Tnp1 - 2) * Psi1 / Dx - Psi0  # np.float64
        Chi = (Tnp1 - 2) * Chi1 / Dx - Chi0  # np.float64
        Xi = np.complex128(Psi + 1j * Chi)  # np.complex64

        A = ((D[N] * Ir + Rnx) * Psi - Psi1) / (
            (D[N] * Ir + Rnx) * Xi - Xi1
        )  # np.complex128
        B = ((D[N] * SCm + Rnx) * Psi - Psi1) / (
            (D[N] * SCm + Rnx) * Xi - Xi1
        )  # np.complex128

        Dqxt = Dqxt + Tnp1 * (A + B).real  # np.float64
        Dqsc = Dqsc + Tnp1 * (A * np.conj(A) + B * np.conj(B)).real  # np.float64

        APB, AMB = A2 * (A + B), A2 * (A - B)  # np.complex128

        S = Dqv * Pi1  # array, np.float64
        T = S - Pi0  # array np.float64
        Taun = N * T - Pi0  # array, np.float64
        Sp = Sp + APB * (Pi1 + Taun)  # array, np.complex128
        Sm = Sm + AMB * (Pi1 - Taun)  # array, np.complex128
        Pi0, Pi1 = Pi1, S + T * Turbo  # array, np.float64

        Psi0, Psi1 = Psi1, Psi  # np.float64
        Chi0, Chi1 = Chi1, Chi  # np.float64
        Xi1 = np.complex128(Psi1 + 1j * Chi1)  # np.complex128

    Dx2 = Dx**2  # np.float64
    Dqsc = 2 * Dqsc / Dx2  # np.float64
    Dqxt = 2 * Dqxt / Dx2  # np.float64
    AA = 2 / (Dx2 * Dqsc)  # np.float64

    Xs1, Xs2 = 0.5 * (Sp + Sm), 0.5 * (Sp - Sm)  # array, np.complex128
    pf = AA * (np.abs(Xs1) ** 2 + np.abs(Xs2) ** 2)  # array, np.float64

    return Dqxt, Dqsc, pf


# @utils.show_runtime
def get_quad(
    legendre_coefficients_flag=False, phase_quad_type="L", phase_quad_N=181, angle=None
):
    """Calculates cosines of the angles to return the phase function at.

    If Legendre expansion of the phase function is also requested
    (i.e. if legendre_coefficients_flag == True)
    the user angles are ignored. The quadrature points are calculated
    by the Quadrature module.
    There are several possible quadratures, here Lobatto is used by default and
    the angles are from 0 to 180 in steps of 1 degree.
    elif the Legendre expansion of the phase function is not required and angles are
    given, simply calculates the cosines,
    lastly, elif Legendre expansion is not needed and angles are not given, one angle
    is used.
    If angles are given, calculates the phase function at those angles otherwise use
    default 0 to 180 in steps of 1 degree (and flip).

    Args:
        legendre_coefficients_flag (bool): If True, signifies that Legendre polynomial
            expansion is requested. Default is False.
        angle (array-like, float, int): Scattering angles. If None, then phase function
            returned at (0,180) degrees with one degree spacing. Can't be set when
            Legendre polynomial expansion of the phase function is requested (in the
            case quadrature angles are computed and used). Default is None.
        phase_quad_N (int): Number of phase function moments (scattering angles). Used
            when Legendre polynomial expansion requested. Default is 181.
        phase_quad_type (str): when Legendre polynomial expansion requested, this type
            of quadrature is used. Can be "G" (Gaussian), "R" (Radau) or "L" (Lobatto).
            Default is L.

    Returns:
        cos_angle_value (array-like): Cosines of scattering angle value of the
            quadrature.
        cos_angle_weight (array-like): Weight of the quadrature points.
        angles (int): Number of scattering angles (quarature points).

    """
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
            cos_angle_weight = np.zeros(cos_angle_value.size, dtype=np.float64)
        else:
            cos_angle_value = 0
            cos_angle_weight = 1  # 1 angle

    # ensure correct type for Mie call
    cos_angle_value = np.array(cos_angle_value, dtype=np.float64)

    cos_angle_weight = np.array(cos_angle_weight, dtype=np.float64)

    angles = np.int64(cos_angle_value.size)

    return cos_angle_value, cos_angle_weight, angles


def get_radii(distribution, eta=1e-6, radii_quad_type="T", radii=200):
    """Calculates the lower and upper particle radii from the distribution.

    The upper and lower particle radii are calculated at points where
    :math:`n(r) = eta r_mode`.
    where *eta* is the cut-off and *r_mode* is the mode of the distribution.

    Args:
        distribution (obj): An instance of `srfm.size_distribution.SizeDistribution`.
        eta (int, float): Size distribution cut-off value. The size distribution is
            a function that technically spans the (-inf,+inf) size interval. The eta
            value is a value of the size distribution beyond whose corresponding
            size (radius) the distribution is truncated. Default is 1e-6.
        radii_quad_type (str): Quadrature type for the radii calculation function.
            Accetped values are values implemented in the ``srfm.quadrature``
            module. Currently implemented (Apr 2025) are "L" (Lobatto quadrature
            rule), "G" (Gauss), "R" (Radau) and "T" (Trapezium). Default is T.
        radii (int): Number of radii in size distribution. Default is 200.

    Returns:
        radius (array-like): size distribution quadrature points
        radius_weight (array-like): quadrature weights
        red_wt_sum (float): sum of radius_weight, used to normalize relevant sums in
            the code.

    """

    eta = np.float64(eta)

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
    radii = np.int64(radii)

    radius, weight_temp = quad.quadrature(
        radii_quad_type, radii, radius_lower_bound, radius_upper_bound
    )
    radius_weight = weight_temp * distribution.value(radius)
    rad_wt_sum = np.sum(radius_weight)

    return radius, radius_weight, rad_wt_sum


def get_ri(composition, refractive_index=None, wave=None, wave_size=None):
    """Load refractive indices.

    Indices expected in the form (n - ik).

    Args:
        composition (str): One of accepted values by
            ``ARIA_module.get_ri_filepathname()``.
        refractive_index (array-like): If compositon is "ri", then refractive indices
            are required from the user. Default is None.
        wave (array-like): Wavelength grid. Default is None.
        wave_size (int): size of wave (equal to len(wave) is wave is list or 1D array
            and wave.size if wave is an array).

    Returns:
        refractive_index_arr (array-like): Array of refractive indices. dtype *complex*.

    Raises:
        RuntimeError: Raised when composition is "ri", but refractive index is not
            given.
        ValueError: Raised when and imaginary part of a refractive index is > 0, because
            indices are expected in the form (n - ik).

    """

    if composition == "ri":
        if refractive_index is None:
            raise RuntimeError("Error: refractive index not defined.")
        else:
            ri_n = np.array(np.full(wave_size, np.real(refractive_index)))
            ri_k = np.array(np.full(wae_size, np.imag(refractive_index)))
    else:
        ri_object = ARIA.RI()
        ri_n, ri_k = ri_object.load_refractive_indices(
            composition, wave=wave, mode="wavelength"
        )
        ri_k = -ri_k

    # check refractive index before proceeding
    if np.any(ri_k > 0):
        raise ValueError(
            """Refractive index (k) contains invalid values. 
               Ensure all values are zero or negative."""
        )

    refractive_index_arr = np.array(ri_n + 1j * ri_k, dtype=np.complex128)
    return refractive_index_arr


# @jit(nopython=True, error_model="numpy", parallel=False, fastmath=True)
def loop_mie_over_radii(radii, sp, ri, angles, cos_angle_value):
    """Calculate the Mie scattering for a given wavelength over a range of radii.

    Args:
        radii (int): Number of size distribution quadrature radii.
        sp (array-like): Array of size parameters.
        ri (array-like): Array of refractive indices.
        angles (int): Number of scattering angles (quadrature points).
        cos_angle_value (array-like): Cosines of scattering angles at which to calculate
            the scattering properties, comes from quadrature.

    Returns:
        Q_ext_value (array): Extinction coefficients array, shape (radii,).
        Q_sca_value (array): Scattering coefficients array, shape (radii,).
        phase_function_value (array): Phase function array, shape (radii, angles).

    """

    # Initialize arrays with the same size as 'radii'
    Q_ext_value = np.zeros(
        radii, dtype=np.float64
    )  # 1D array for extinction per radius
    Q_sca_value = np.zeros(
        radii, dtype=np.float64
    )  # 1D array for scattering per radius
    phase_function_value = np.zeros(
        (radii, angles), dtype=np.float64
    )  # 2D array for phase function per radius

    for j in range(radii):
        # Perform Mie calculations
        Q_ext_value[j], Q_sca_value[j], phase_function_value[j, :] = (
            mie_module.mie_ewp(sp[j], ri, cos_angle_value)
        )

    return Q_ext_value, Q_sca_value, phase_function_value


# @jit(nopython=True, error_model="numpy", parallel=False, fastmath=True)
@utils.show_runtime
def regrid(op_dict, wvls, track_diff=False, diff_type="pct"):
    """Linearly interpolates calculated values from ewp_hs to a new grid.

    Args:
        op_dict (dict): output from ``srfm.optical_properties.ewp_hs``. Should have
            keys:

            - beta_ext
            - ssalb
            - phase_function
            - legendre_coefficient
            - wavelengths

            If ewp_hs changes in the future, this function needs to change as well.

        wvls (array-like): new grid, [\ :math:`\\mu`\ m]
        track_diff (bool): If True, calculates differences arising from interpolation

    Returns:
        op_dict (dict): New version of op_dict on a new grid.
        diff_dict (dict): Difference between the interpolated and non-interpolated
            dictionaries. Same keys as op_dict.

    Raises:
        RuntimeError: Raised when invalid key encountered in op_dict.
        ValueError: Raised when wavelength grid is not monotonic.

    """
    # check for invalid keys (basically safeguard if ewp_hs changes in the future
    valid_keys = [
        "beta_ext",
        "ssalb",
        "phase_function",
        "legendre_coefficient",
        "wavelengths",
    ]

    for key in op_dict.keys():
        if key not in valid_keys:
            raise RuntimeError(
                """Invalid key encountered in op_dict. Has ewp_hs changed
            since regriding was implemented?"""
            )

    # check if wvls is monotonic and increasing or decreasing
    wvls_mono = utils.monotonic(wvls)
    if wvls_mono == 0:
        raise ValueError("Wvls is not monotonic.")
    elif wvls_mono == 2:
        wvls = np.flip(wvls)

    dict_mono = utils.monotonic(op_dict["wavelengths"])
    if dict_mono == 0:
        raise ValueError("op_dict is not monotonic.???")
    elif dict_mono == 2:
        op_dict["wavelengths"] = np.flip(op_dict["wavelengths"])
        op_dict["beta_ext"] = np.flip(op_dict["beta_ext"])
        op_dict["ssalb"] = np.flip(op_dict["ssalb"])
        op_dict["phase_function"] = np.flipud(op_dict["phase_function"])
        op_dict["legendre_coefficient"] = np.flipud(op_dict["legendre_coefficient"])

    diff_dict = {}
    if track_diff == True:
        old_op_dict = op_dict.copy()

    op_dict["beta_ext"] = np.interp(wvls, op_dict["wavelengths"], op_dict["beta_ext"])
    op_dict["ssalb"] = np.interp(wvls, op_dict["wavelengths"], op_dict["ssalb"])
    op_dict["phase_function"] = np.array(
        [
            np.interp(wvls, op_dict["wavelengths"], op_dict["phase_function"][:, i])
            for i in range(op_dict["phase_function"].shape[1])
        ]
    ).T
    op_dict["legendre_coefficient"] = np.array(
        [
            np.interp(
                wvls, op_dict["wavelengths"], op_dict["legendre_coefficient"][:, i]
            )
            for i in range(op_dict["legendre_coefficient"].shape[1])
        ]
    ).T
    op_dict["wavelengths"] = wvls

    if track_diff == True:
        diff_dict = track_regrid_diff(old_op_dict, op_dict)

    if wvls_mono == 1:
        return op_dict, diff_dict
    elif wvls_mono == 2:
        op_dict["wavelengths"] = np.flip(op_dict["wavelengths"])
        op_dict["beta_ext"] = np.flip(op_dict["beta_ext"])
        op_dict["ssalb"] = np.flip(op_dict["ssalb"])
        op_dict["phase_function"] = np.flipud(op_dict["phase_function"])
        op_dict["legendre_coefficient"] = np.flipud(op_dict["legendre_coefficient"])
        diff_dict["wavelengths"] = np.flip(op_dict["wavelengths"])
        diff_dict["beta_ext"] = np.flip(op_dict["beta_ext"])
        diff_dict["ssalb"] = np.flip(op_dict["ssalb"])
        diff_dict["phase_function"] = np.flipud(op_dict["phase_function"])
        diff_dict["legendre_coefficient"] = np.flipud(op_dict["legendre_coefficient"])
        return op_dict, diff_dict


# @jit(nopython=True, error_model="numpy", parallel=False, fastmath=True)
def loop_mie_over_wavelengths(
    wavelengths_size,
    radii,
    size_parameters,
    refractive_index_arr,
    angles,
    cos_angle_value,
    cos_angle_weight,
    radius_weight,
    radius,
    rad_wt_sum,
    legendre_coefficients_flag,
    legendre_coefficients_type,
    beta_ext,
    beta_sca,
    phase_function,
    legendre_coefficient,
):
    """Main computational loop for ``srfm.optical_properties.ewp_hs()``.

    Loops over wavelengths and calcualtes optical properties.

    Args:
        wavelengths_size (int): Number of wavelengths to loop over.
        radii (int, float): Number of radii of the particles in the distribution.
            Default is 200.
        size_parameters (array-like): Array of size parameters.
        refractive_index_arr (array-like): Array of refractive indices. dtype *complex*.
        angles (int): Number of scattering angles (quarature points).
        cos_angle_value (array-like): Cosines of scattering angle value of the
            quadrature.
        cos_angle_weight (array-like): Weight of the quadrature points.
        radius (array-like): size distribution quadrature points
        radius_weight (array-like): quadrature weights
        red_wt_sum (float): sum of radiu_weight, used to normalize relevant sums in
            the code.
        legendre_coefficients_flag (bool): If True, Legendre polynomial expansion of the
            phase function is performed and expansion coefficients are returned.
        legendre_coefficients_type (str): Specifies the type of requested Legendre
            polynomial expansion coefficients. Can be *normalised* or *regular*. Default
            is *normalised*.
        beta_ext (array): Empty array to store extinction coefficient in, shape
            (wavelengths_size,).
        beta_sca (array): Empty array to store extinction coefficient in, shape
            (wavelengths_size,).
        phase_function (array): Empty array to store the phase function in, shape
            (wavelengths_size, angles).
        legendre_coefficient (array): Empty array to store Legendre coefficients in,
            shape (wavelengths_size, _), where the second dimension is big enough to
            contain the Legendre coefficients. This needs to be guessed and eventually
            adapted by the user, because if the wavelengths_size is too big, then
            the total array size can cause overflows in memory. It is possible to use
            ``srfm.utilities.memory_safe_np_zeros_2d()`` to generate a 2D array
            with dimensions adapted to not overflow the available memory, but be advised
            that this may cause the array to have the second dimension too small,
            resulting in truncation of Legendre coefficients and loss of precision.
            Will not cause problems on high-RAM stations. Mind that in a np.zeros array
            the default data type is np.float64, which has a size of 8 bytes. The array
            size is number of elements*8 in bytes.

    Raises:
        ValueError: Raised when invalid input into legendre_coeffs_type is detected.

    """

    max_lc = 0  # tracks length of legendre expansion in the main computational loop
    # (number of Legendre polynomial coefficients used. Used to truncate the final
    # array.

    for i in range(wavelengths_size):

        # Report progress as percentage completed
        pct = [
            0,
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
        pct_val = [np.int64(i * wavelengths_size / 100) for i in pct]

        if i in pct_val:
            val_idx = pct_val.index(i)
            print(f"Calculating particle optical properties. {pct[val_idx]}% done...")

#        Q_ext_value, Q_sca_value, phase_function_value = loop_mie_over_radii(
#                                            radii=radii,
#                                            sp=size_parameters[:,i],
#                                            ri=refractive_index_arr[i],
#                                            angles=angles,
#                                            cos_angle_value=cos_angle_value
#                                        )

        Q_ext_value, Q_sca_value, phase_function_value, error = mie_module.mie_ewp(
            size_parameters[:, i], refractive_index_arr[i], cos_angle_value
        )

        # calculate weighted sums of some variables:

        # Weighted sums for extinction and scattering
        beta_ext[i] = np.pi * np.sum(radius_weight * radius**2 * Q_ext_value) * 1e-6
        beta_sca[i] = np.pi * np.sum(radius_weight * radius**2 * Q_sca_value) * 1e-6

        # Weighted sums for phase function
        for k in range(angles):
            phase_function[i, k] = (
                np.pi
                * np.sum(
                    radius_weight * radius**2 * Q_sca_value * phase_function_value[:, k]
                )
                / (beta_sca[i] * 1e6)
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
    return beta_ext, beta_sca, phase_function, legendre_coefficient, max_lc


def track_regrid_diff(old_dict, new_dict, diff_type="pct"):
    """Calculates the difference before and after interpolation.

    Called by optical_properties.regrid().

    Args:
        old_dict (dict): op_dict before interpolation.
        new_dict (dict): op_dict after interpolation.
        diff_type (str): Can be "pct" (differences in percent from the old value) or
            "abs" - absolute difference. Default is "pct".

    Raises:
        ValueError: Raised when invalid diff_type encountered.

    """
    diff_dict = {}
    diff_dict["wavelengths"] = new_dict["wavelengths"]

    temp_dict = {}
    temp_dict["wavelengths"] = new_dict["wavelengths"]
    temp_dict["beta_ext"] = []
    temp_dict["ssalb"] = []
    temp_dict["phase_function"] = []
    temp_dict["legendre_coefficient"] = []

    for i in new_dict["wavelengths"]:
        x = np.argmin(
            np.abs(old_dict["wavelengths"] - np.full(old_dict["wavelengths"].size, i))
        )
        temp_dict["beta_ext"].append(old_dict["beta_ext"][x])
        temp_dict["ssalb"].append(old_dict["ssalb"][x])
        temp_dict["phase_function"].append(old_dict["phase_function"][x])
        temp_dict["legendre_coefficient"].append(old_dict["legendre_coefficient"][x])

    for key in temp_dict.keys():
        if key == "wavelengths":
            pass
        else:
            temp_dict[key] = np.asarray(temp_dict[key])

            if diff_type == "abs":
                diff_dict[key] = np.abs(temp_dict[key] - new_dict[key])
            elif diff_type == "pct":
                numerator = np.abs(temp_dict[key] - new_dict[key])
                diff_dict[key] = (
                    np.nan_to_num(np.divide(numerator, temp_dict[key])) * 100
                )
            else:
                raise ValueError("diff_type must be pct or abs.")

    return diff_dict


def calc_op_diff(first_dict, sec_dict, diff_type="pct", plot=True):
    """Calculates the difference between optical properties.

    Used to get differences coming from different sources, i.e. one that's calculated
    from ewp_hs at higher resolution and one that's caluclated at lower resolution and
    interpolated.
    Assumes identical wavelength grids.

    Args:
        first_dict (dict): dictionary with optical properties, op_dict from ewp_hs
            or equivalent
        sec_dict (dict): dictionary with optical properties, op_dict from ewp_hs
            or equivalent
        diff_type (str): Can be "pct" (differences in percent from the old value) or
            "abs" - absolute difference. Default is "pct".
        plot (bool): If True, makes a plot.

    Raises:
        ValueError: Raised when a key is found in one dictionary that is missing from
            the other one (i.e. the two input dictionaries must have the same keys).
        ValueError: Raised when incorrect value enountered in diff_type.
    """

    expected_keys = [
        "wavelengths",
        "beta_ext",
        "ssalb",
        "phase_function",
        "legendre_coefficient",
    ]

    for k in expected_keys:
        if k not in first_dict.keys() and k not in sec_dict.keys():
            raise ValueError(f"{k} not found in both dictionaries.")

    assert (first_dict["wavelengths"] == sec_dict["wavelengths"]).all()

    diff_dict = {}
    diff_dict["wavelengths"] = first_dict["wavelengths"]

    for key in first_dict.keys():
        if key == "wavelengths":
            pass
        else:
            if diff_type == "abs":
                diff_dict[key] = np.abs(first_dict[key] - sec_dict[key])
            elif diff_type == "pct":
                numerator = np.abs(sec_dict[key] - first_dict[key])
                diff_dict[key] = (
                    np.nan_to_num(np.divide(numerator, first_dict[key])) * 100
                )
            else:
                raise ValueError("diff_type must be pct or abs.")
    if plot == True:
        fig = plot_diff_dict(diff_dict)
        return diff_dict, fig
    else:
        return diff_dict


def plot_diff_dict(diff_dict, **kwargs):
    """Plots diff_dict - difference in calcualted optical properties.

    Args:
        diff_dict (dict): Dictionary with differences in optical properties, comes e.g.
            from calc_op_diff or track_regrid_diff.

    Returns:
        fig (object): Figure object.

    """
    fig = plt.figure(figsize=(11.7, 8.3))

    plt.subplot(3, 1, 1)
    plt.plot(diff_dict["wavelengths"], diff_dict["beta_ext"], label="beta_ext")
    plt.xlabel("wavelength (um)", fontsize="12")
    plt.ylabel(r"$\Delta \beta ^{ext} (\lambda)$ (%)", fontsize="12")
    # plt.yscale("log")
    plt.title("Extinction coefficient", fontsize="15")

    plt.subplot(3, 1, 2)
    plt.plot(diff_dict["wavelengths"], diff_dict["ssalb"], label="ssalb")
    plt.xlabel("wavelength (um)", fontsize="12")
    plt.ylabel(r"$\Delta \omega (\lambda)$ (%)", fontsize="12")
    # plt.yscale("log")
    plt.title("Single scattering albedo", fontsize="15")

    plt.subplot(3, 1, 3)
    [
        plt.plot(
            diff_dict["wavelengths"],
            diff_dict["phase_function"][:, i],
            label=f"phase function, {i} deg",
        )
        for i in range(diff_dict["phase_function"].shape[1])
    ]
    plt.xlabel("wavelength (um)", fontsize="12")
    plt.ylabel(r"$\Delta p(\mu, \lambda)$ (%)", fontsize="12")
    # plt.yscale("log")
    # plt.legend()
    plt.title("Phase function - each line represents one angle", fontsize="15")

    # plt.subplot(4,1,4)
    # [plt.plot(diff_dict["wavelengths"],diff_dict["legendre_coefficient"][:,i],label=f"legendre coefficient no. {i}") for i in range(diff_dict["legendre_coefficient"].shape[1])]
    # plt.xlabel("wavelength (um)")
    # plt.ylabel("legendre coefficients")
    # plt.yscale("log")
    ##plt.legend()

    plt.suptitle(
        "Difference in optical properties arising from interpolation", fontsize="18"
    )
    plt.tight_layout()
    return fig
