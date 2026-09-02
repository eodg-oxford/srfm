"""Calculates quadrature points.

Can do Gaussian, Lobatto, Radau, and Trapezium quadratures.

- Name: quadrature
- Parent package: srfm
- Author: Don Grainger
- Contributors: Antonin Knizek
- Date: 24 January 2025
"""

import numpy as np
from numba import jit_module


# Function to calculate BesselZero
def bessel_zero(s):
    """Calculates BesselZero.

    Args:
        s (float): ?

    Returns:
        bslz (float): ?

    """
    A1 = -0.60792710185402662866327677925836
    A2 = 6.159589352810601113491669960271217e-2
    A3 = -0.981080301612647885671934857251079
    A4 = 11.7507817547326698965409187559573
    A5 = -2731.24978203593727776707267899459

    b = 4 * s + 1
    bslz = (
        0.25
        * np.pi
        * b
        * (1 + A1 / b**2 + A2 / b**4 + A3 / b**6 + A4 / b**8 + A5 / b**10)
    )
    return bslz


# Function to calculate FirstGuess
def first_guess(quad_type, n, term):
    """Function to calculate first guess of the abscissa points.

    Args:
        quad_type (str): Quadrature type. Accepted values are "G" (Gaussian), "R"
            (Radau), and "L" (Lobatto).
        n (float): ?
        term (float): ?

    Returns:
        fg (float): First guess of the abscissa points.

    """

    quad_type = quad_type.upper()
    if quad_type == "G":
        return np.cos(np.pi * (term - 0.25) / (n + 0.5))
    elif quad_type == "R":
        return (
            np.cos(np.pi * (term - 0.25) / (n + 0.5))
            + np.cos(np.pi * (term - 0.25) / (n - 0.5))
        ) / 2
    elif quad_type == "L":
        return np.cos(
            bessel_zero(term) / np.sqrt((n - 0.5) ** 2 + (0.25 - 1 / np.pi**2))
        )
    else:
        raise ValueError(f"Invalid quadrature type: {quad_type}")


# Function to calculate Newton correction
def newton_g(quad_type, n, x):
    """Performs Newton correction on the abscissa.

    Args:
        quad_type (str): Quadrature type. Accepted values are "G" (Gaussian), "R"
            (Radau), and "L" (Lobatto).
        n (float): ?
        x (float): ?

    Returns:

    """

    quad_type = quad_type.upper()
    if quad_type == "G":
        pl, pm, pn = legendre(n, x)
        return (1 - x**2) * pn / (n * (pm - x * pn))
    elif quad_type == "R":
        pl, pm, pn = legendre(n, x)
        return (
            (1 + x)
            * (pm + pn)
            / (
                ((n - 1) * pl - (n - 1) * x * pm + n * pm - n * x * pn) / (1 - x)
                - (pm + pn)
            )
        )
    elif quad_type == "L":
        pl, pm, pn = legendre(n - 1, x)
        return (
            (1 - x**2)
            * (pm - x * pn)
            / (n * pl + 2 * (1 - n) * x * pm + (x**2 * (n - 1) - 1) * pn)
        )
    else:
        raise ValueError(f"Invalid quadrature type: {quad_type}")


# Function to calculate Legendre polynomials
def legendre(n, x):
    """Calculate the Legenedre polynomials using Bonnet's recursion formula.

    Args:
        n (int): number of points in the quadrature
        x (float): specific point in the quadrature

    Returns:
        pl (float): Legendre polynomial :math:`P_{n-2}(x)`.
        pm (float): Legendre polynomial :math:`P_{n-1}(x)`.
        pn (float): Legendre polynomial :math:`P_{n}(x)`.

    """

    if n == 0:
        return None, None, 1.0
    elif n == 1:
        return None, 1.0, x
    else:
        pl = 1.0
        pm = 1.0
        pn = x
        for i in range(2, n + 1):
            coeff = 1 / float(i)
            pl, pm = pm, pn
            pn = (2 - coeff) * x * pm - (1 - coeff) * pl
        return pl, pm, pn


# translated from Don Grainger's IDL by Antonin Knizek
def quadrature101(quad_type, npts):
    """Asign npts abscissae and weights for integration on the interval [-1,1] for
    a variety of quadrature types.

    Args:
        quad_type (str): Quadrature type. Accepted values are "G" (Gaussian), "R"
            (Radau), "L" (Lobatto), "T" (Trapezium), "S" (Simpson).
        npts (int): Number of abscissa points.

    Returns:
        abscissa (array): Calculated abscissa.
        weight (array): Corresponding weights.

    """

    if not isinstance(npts, (int, np.integer)):
        raise TypeError("npts must be an integer.")
    if npts < 1:
        raise ValueError("npts must be positive.")
    if npts > 20000:
        raise ValueError("Too many quadrature points.")

    sigfig = 14  # precision limit (idl legacy, but value kept)

    quad_type = quad_type.upper()

    if quad_type == "T" and npts < 2:
        raise ValueError("Trapezium quadrature requires at least two points.")
    if quad_type == "L" and npts < 2:
        raise ValueError("Lobatto quadrature requires at least two points.")

    if quad_type == "T":
        abscissa = -1 + 2 * np.arange(npts, dtype=np.float64) / (npts - 1)
        weight = np.full(npts, (2 / (npts - 1)), dtype=np.float64)
        weight[0] = 1 / (npts - 1)
        weight[-1] = 1 / (npts - 1)

    elif quad_type == "S":
        raise ValueError("Simpson method not implemented yet.")

    elif quad_type in ["G", "R", "L"]:
        abscissa = np.zeros(npts, dtype=np.float64)  # quadrature points
        weight = np.zeros(npts, dtype=np.float64)  # quadrature weigths

        # set-up boundary conditions
        if quad_type == "G":
            term = 0
            zeros = npts
        elif quad_type == "R":
            abscissa[-1] = -1
            weight[-1] = 2 / npts**2
            term = 0
            zeros = npts - 1
        elif quad_type == "L":
            abscissa[0] = 1
            weight[0] = 2 / (npts * (npts - 1))
            abscissa[-1] = -1
            weight[-1] = 2 / (npts * (npts - 1))
            term = 1
            zeros = npts - 2

        # loop over points in the quadrature
        for zero in range(1, zeros + 1):
            lx = float("inf")
            x = first_guess(quad_type, npts, zero)

            while abs(x - lx) > 10 ** (-sigfig):
                lx = x
                x = x - newton_g(quad_type, npts, x)

            abscissa[term] = x
            _, pm, _ = legendre(npts, x)
            if quad_type == "G":
                weight[term] = 2 * (1 - x**2) / (npts * pm) ** 2
            elif quad_type == "R":
                weight[term] = (1 - x) / (npts * pm) ** 2
            elif quad_type == "L":
                weight[term] = 2 / (npts * (npts - 1) * pm**2)

            term += 1

        abscissa = abscissa[::-1]
        weight = weight[::-1]

    else:
        raise ValueError("Error in quadrature: Invalid quadrature type.")

    return abscissa, weight


def shift_quadrature(abscissa, weight, lower_bound, upper_bound):
    """Shifts quadrature abscissa and weight.

    Shifts from the interval [-1, 1] to [lower_bound, upper_bound].

    Args:
        abscissa (array): Quadrature points on the interval [-1, 1].
        weight (array): Quadrature weight on the interval [-1, 1].
        lower_bound (float): Lower bound of the new interval.
        upper_bound (float): Upper bound of the new interval.

    Returns:
        new_abscissa (array): New, shifted abscissa
        new_weights (array): New, shifted weights.

    """
    new_abscissa = (
        (lower_bound + upper_bound) + (upper_bound - lower_bound) * abscissa
    ) / 2
    new_weight = (upper_bound - lower_bound) * weight / 2
    return new_abscissa, new_weight


def quadrature(quad_type, n_pts, lower_bound, upper_bound):
    """Calls the appropriate quadrature function.

    First calls a function to calculate the quadrature on a [-1,1] interval, then shifts
    the calculated quadrature and weigths.

    Args:
        quad_type (str): Quadrature type. Accepted values are "G" (Gaussian), "R"
            (Radau), "L" (Lobatto), "T" (Trapezium), "S" (Simpson).
        npts (int): Number of abscissa points.
        lower_bound (float): Lower bound of the new interval.
        upper_bound (float): Upper bound of the new interval.

    Returns:
        abscissa_new (array): New, shifted abscissa
        weight_new (array): New, shifted weights.

    """
    abscissa, weight = quadrature101(quad_type, n_pts)
    abscissa_new, weight_new = shift_quadrature(
        abscissa, weight, lower_bound, upper_bound
    )
    return abscissa_new, weight_new


# jit_module(nopython=True, error_model="numpy", parallel=False, fastmath=True)
