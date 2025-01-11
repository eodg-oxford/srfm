# quadrature.py

import numpy as np

# Constants
PI = np.pi

# Function to calculate BesselZero
def bessel_zero(s):
    A1 = -0.60792710185402662866327677925836
    A2 = 6.159589352810601113491669960271217e-4
    A3 = -0.981080301612647885671934857251079
    A4 = 11.7507817547326698965409187559573
    A5 = -2731.24978203593727776707267899459

    b = 4 * s + 1
    bslz = 0.25 * PI * b * (1 + A1 / b**2 + A2 / b**4 + A3 / b**6 + A4 / b**8 + A5 / b**10)
    return bslz

# Function to calculate FirstGuess
def first_guess(quad_type, n, term):
    quad_type = quad_type.upper()
    if quad_type == 'G':
        return np.cos(PI * (term - 0.25) / (n + 0.5))
    elif quad_type == 'R':
        return (np.cos(PI * (term - 0.25) / (n + 0.5)) +
                np.cos(PI * (term - 0.25) / (n - 0.5))) / 2
    elif quad_type == 'L':
        return np.cos(bessel_zero(term) / np.sqrt((n - 0.5)**2 + (0.25 - 1 / PI**2)))
    else:
        raise ValueError(f"Invalid quadrature type: {quad_type}")

# Function to calculate Newton correction
def newton_g(quad_type, n, x):
    quad_type = quad_type.upper()
    if quad_type == 'G':
        pl, pm, pn = legendre(n, x)
        return (1 - x**2) * pn / (n * (pm - x * pn))
    elif quad_type == 'R':
        pl, pm, pn = legendre(n, x)
        return ((1 + x) * (pm + pn) /
                (((n - 1) * pl - (n - 1) * x * pm + n * pm - n * x * pn) / (1 - x) - (pm + pn)))
    elif quad_type == 'L':
        pl, pm, pn = legendre(n - 1, x)
        return ((1 - x**2) * (pm - x * pn) /
                (n * pl + 2 * (1 - n) * x * pm + (x**2 * (n - 1) - 1) * pn))
    else:
        raise ValueError(f"Invalid quadrature type: {quad_type}")

# Function to calculate Legendre polynomials
def legendre(n, x):
    if n == 0:
        return None, None, 1.0
    elif n == 1:
        return None, 1.0, x
    else:
        pl = 1.0
        pm = 1.0
        pn = x
        for i in range(2, n + 1):
            coeff = 1 / i
            pl, pm = pm, pn
            pn = (2 - coeff) * x * pm - (1 - coeff) * pl
        return pl, pm, pn

# Main Quadrature function on interval [-1,1]
def quadrature101(quad_type, n_pts):
    n = float(n_pts)
    abscissa = np.zeros(n_pts)
    weight = np.zeros(n_pts)
    sigfig = 14  # Precision for convergence

    quad_type = quad_type.upper()
    if quad_type == 'T':
        abscissa = -1 + 2 * np.arange(n_pts) / (n - 1)
        weight[:] = 2 / (n - 1)
        weight[0] = weight[-1] = 1 / (n - 1)
        return abscissa, weight

    elif quad_type == 'S':
        raise NotImplementedError("Simpson method is not implemented.")

    elif quad_type in ['G', 'R', 'L']:
        term = 0
        if quad_type == 'R':
            abscissa[-1] = -1
            weight[-1] = 2 / (n * n)
        elif quad_type == 'L':
            abscissa[0] = 1
            weight[0] = 2 / (n * (n - 1))
            abscissa[-1] = -1
            weight[-1] = 2 / (n * (n - 1))
            term = 1

        for zero in range(term, n_pts - (1 if quad_type == 'R' else 0)):
            x = first_guess(quad_type, n_pts, zero + 1)
            lx = float('inf')
            while abs(x - lx) > 10**-sigfig:
                lx = x
                x = x - newton_g(quad_type, n_pts, x)

            abscissa[term] = x
            _, pm, _ = legendre(n_pts, x)
            if quad_type == 'G':
                weight[term] = 2 * (1 - x**2) / (n * pm)**2
            elif quad_type == 'R':
                weight[term] = (1 - x) / (n * pm)**2
            elif quad_type == 'L':
                weight[term] = 2 / (n * (n - 1) * pm**2)
            term += 1

        abscissa = abscissa[::-1]
        weight = weight[::-1]
        return abscissa, weight

    else:
        raise ValueError(f"Invalid quadrature type: {quad_type}")
        
def shift_quadrature(abscissa, weight, lower_bound, upper_bound):
    """
    Shifts quadrature abscissa and weight from the interval [-1, 1] to [lower_bound, upper_bound].

    Parameters:
        abscissa (numpy.ndarray): Quadrature points on the interval [-1, 1].
        weight (numpy.ndarray): Quadrature weight on the interval [-1, 1].
        lower_bound (float): Lower bound of the new interval.
        upper_bound (float): Upper bound of the new interval.

    Returns:
        tuple: (new_abscissa, new_weights)
    """
    new_abscissa = ((lower_bound + upper_bound) + (upper_bound - lower_bound) * abscissa) / 2
    new_weight = (upper_bound - lower_bound) * weight / 2
    return new_abscissa, new_weight
    
def quadrature(quad_type, n_pts, lower_bound, upper_bound):
    abscissa, weight = quadrature101(quad_type, n_pts)
    abscissa_new, weight_new = shift_quadrature(abscissa, weight, lower_bound, upper_bound)
    return abscissa_new, weight_new
