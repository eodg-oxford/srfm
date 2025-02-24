"""
Name: units
Parent package: srfm
Author: Antonin Knizek
Contributors: 
Date: 18 February 2025
Purpose: Provides functions for unit conversion.
""" 
def decimal_degree_to_DMS(dec):
    negative = dec < 0
    dec = abs(dec)
    M, S = divmod(dec * 3600, 60)
    D, M = divmod(M, 60)
    if negative:
        if D < 0:
            D = -D
        elif M < 0:
            M = -M
        elif S < 0:
            S = -S
    return (int(D), int(M), int(S))


def inv_cm_to_micron(wavenumber):
    """convert wavenumber in cm-1 to wavelength in micrometers"""
    return (1 / wavenumber) * 1e4


def inv_cm_to_nm(wavenumber):
    """convert wavenumber in cm-1 to wavelength in nm"""
    return (1 / wavenumber) * 1e7


def micron_to_inv_cm(wavelength):
    """convert wavelength in micrometers to wavenumbers in cm-1"""
    return 1 / wavelength * 1e4


def nm_to_inv_cm(wavelength):
    """convert wavelength in nm to wavenumbers in cm-1"""
    return 1 / wavelength * 1e7
