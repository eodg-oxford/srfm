"""
Name: quadrature
Parent package: srfm
Author: Don Grainger
Contributors: 
Date: 24 January 2025
Purpose: Calculates quadrature points, can do Gaussian, Lobatto and Radau quadratures
""" 

import numpy as np
from numba import jit_module

# Function to calculate BesselZero
def bessel_zero(s):
    A1 = -0.60792710185402662866327677925836
    A2 = 6.159589352810601113491669960271217e-2
    A3 = -0.981080301612647885671934857251079
    A4 = 11.7507817547326698965409187559573
    A5 = -2731.24978203593727776707267899459

    b = 4 * s + 1
    bslz = (
        0.25 * np.pi * b * (1 + A1 / b**2 + A2 / b**4 + A3 / b**6 + A4 / b**8 + A5 / b**10)
    )
    return bslz

def test_bessel_zero():
    assert bessel_zero(3) == 10.173467949597212, "bessel_zero returns incorrect value."
    return

# Function to calculate FirstGuess
def first_guess(quad_type, n, term):
    quad_type = quad_type.upper()
    if quad_type == "G":
        return np.cos(np.pi * (term - 0.25) / (n + 0.5))
    elif quad_type == "R":
        return (
            np.cos(np.pi * (term - 0.25) / (n + 0.5))
            + np.cos(np.pi * (term - 0.25) / (n - 0.5))
        ) / 2
    elif quad_type == "L":
        return np.cos(bessel_zero(term) / np.sqrt((n - 0.5) ** 2 + (0.25 - 1 / np.pi**2)))
    else:
        raise ValueError(f"Invalid quadrature type: {quad_type}")

def test_first_guess():
    assert float(first_guess("G", 181, 50)) == 0.6515842860871356, """Gaussian
        first_guess returns incorrect value."""
    assert float(first_guess("R", 181, 50)) == 0.6497710863213735, """Radau
        first_guess returns incorrect value."""
    assert float(first_guess("L", 181, 50)) == 0.6413165972728538, """Lobatto
        first_guess returns incorrect value."""
    return
    
# Function to calculate Newton correction
def newton_g(quad_type, n, x):
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

def test_newton_g():
    assert newton_g("G", 181, 0.1) == -0.004336062494010145, """Gaussian newton_g
     correction returns incorrect value."""
    assert newton_g("R", 181, 0.1) == 0.0003714827832841103, """Radau newton_g
     correction returns incorrect value."""
    assert newton_g("L", 181, 0.1) == -0.005288777969312314, """Lobatto newton_g
     correction returns incorrect value."""

# Function to calculate Legendre polynomials
def legendre(n, x):
    """Calculate the Legenedre polynomials using Bonnet's recurison formula.
        Inputs:
        n - number of points in the quadrature
        x - specific point in the quadrature
        
        Returns:
        pl == P_(n-2)(x)
        pm == P_(n-1)(x)
        pn == P_(n)(x)        
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

def test_legendre():
    assert legendre(181,0.1) == (0.0456111402767717, 
                                 0.0427731037482122, 
                                 -0.03682815582601351
                                 ), """Function legendre returns incorrect values."""

# translated from Don Grainger's IDL by Antonin Knizek
def quadrature101(quad_type, npts):
    """Asign n[ts abscissae and weights for integration on the interval [-1,1] for
    a variety of quadrature types:
    G: Gauss
    S: Simpson
    T: Trapezium
    R: Radau
    L: Lobatto
    """
    
    if npts > 20000:
        raise ValueError("Too many quadrature points.") 
    
    sigfig = 14 # precision limit (idl legacy, but value kept)
    
    quad_type = quad_type.upper()
    
    if quad_type == "T":
        abscissa = -1 + 2 * np.arange(npts,dtype=np.float64) / (npts-1)
        weight = np.full(npts, (2/(npts-1)),dtype=np.float64)
        weight[0] = 1/(npts-1)
        weight[-1] = 1/(npts-1)
    
    elif quad_type == "S":
        raise ValueError("Simpson method not implemented yet.")
    
    elif quad_type in ["G", "R", "L"]:        
        abscissa = np.zeros(npts,dtype=np.float64) # quadrature points
        weight = np.zeros(npts,dtype=np.float64) # quadrature weigths       
        
        # set-up boundary conditions
        if quad_type == "G":
            term = 0
            zeros = npts
        elif quad_type == "R":
            abscissa[-1] = -1
            weight[-1] = 2 / npts**2
            term = 0
            zeros = npts-1
        elif quad_type == "L":
            abscissa[0] = 1
            weight[0] = 2 / (npts * (npts-1))
            abscissa[-1] = -1
            weight[-1] = 2 / (npts * (npts-1))
            term = 1
            zeros = npts-2
        
        # loop over points in the quadrature
        for zero in range(1,zeros+1):
            lx = float("inf")
            x = first_guess(quad_type, npts, zero)
            
            while abs(x-lx) > 10**(-sigfig):
                lx = x
                x = x - newton_g(quad_type, npts, x)
            
            abscissa[term] = x
            _, pm, _ = legendre(npts, x)
            if quad_type == "G":
                weight[term] = 2 * (1-x**2) / (npts*pm)**2
            elif quad_type == "R":
                weight[term] = (1-x) / (npts*pm)**2
            elif quad_type == "L":
                weight[term] = 2 / (npts * (npts-1) * pm**2)
            
            term += 1
            
        abscissa = abscissa[::-1]
        weight = weight[::-1]
        
    else:
        raise ValueError("Error in quadrature: Invalid quadrature type.")

    return abscissa, weight

def test_quadrature101():
    assert float(quadrature101("G",181)[0][50]) == -0.6383546791355884, """Gaussian
        abscissa in quadrature101 returns a wrong value."""

    assert float(quadrature101("G",181)[1][50]) == 0.013323424039823995, """Gaussian
        weight in quadrature101 returns a wrong value."""

    assert float(quadrature101("R",181)[0][50]) == -0.6431669568316702, """Radau
        abscissa in quadrature101 returns a wrong value."""

    assert float(quadrature101("R",181)[1][50]) == 0.013290800444751914, """Radau
        weight in quadrature101 returns a wrong valunew_wavelengthse."""

    assert float(quadrature101("L",181)[0][50]) == -0.641312349855755, """Lobatto
        abscissa in quadrature101 returns a wrong value."""

    assert float(quadrature101("L",181)[1][50]) == 0.01335472617564372, """Lobatto
        weight in quadrature101 returns a wrong value."""

        
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
    new_abscissa = (
        (lower_bound + upper_bound) + (upper_bound - lower_bound) * abscissa
    ) / 2
    new_weight = (upper_bound - lower_bound) * weight / 2
    return new_abscissa, new_weight

def test_shift_quadrature():
    assert float(
        shift_quadrature(quadrature101("L",181)[0],
                         quadrature101("L",181)[1],
                         0,
                         180
        )[0][50]
    ) == 32.28188851298205, "shift_quadrature returns incorrect abscissa."
    
    assert float(
        shift_quadrature(quadrature101("L",181)[0],
                         quadrature101("L",181)[1],
                         0,
                         180
        )[1][50]
    ) == 1.2019253558079348, "shift_quadrature returns incorrect weight."


def quadrature(quad_type, n_pts, lower_bound, upper_bound):
    abscissa, weight = quadrature101(quad_type, n_pts)
    abscissa_new, weight_new = shift_quadrature(
        abscissa, weight, lower_bound, upper_bound
    )
    return abscissa_new, weight_new

def test_quadrature():
    assert float(quadrature("L",181, 0, 180)[0][50]) == 32.28188851298205, """Function
        quadrature returns incorrect abscissa."""
    
    assert float(quadrature("L",181, 0, 180)[1][50]) == 1.2019253558079348, """Function
        quadrature returns incorrect weight."""
    
if __name__ == "__main__":
    test_bessel_zero()
    test_first_guess()
    test_legendre()
    test_newton_g()
    test_quadrature101()
    test_shift_quadrature()
    test_quadrature()
    
    print("Module quadrature has passed all unit tests.")

#jit_module(nopython=True, error_model="numpy", parallel=False, fastmath=True)    
