"""Provides functions for unit conversion.

- Name: units
- Parent package: srfm
- Author: Antonin Knizek
- Contributors: 
- Date: 18 February 2025 
""" 
def decimal_degree_to_DMS(dec):
    """Convert a value in decimal degrees to degree-minute-seconds.
    
    Args:
        dec (int, float): Value in decimal degrees.
    
    Returns:
        dms (tuple of ints): Tuple containing the result as (degrees, minutes, seconds).
    
    """
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
    
    dms = (int(D), int(M), int(S))
    
    return dms


def inv_cm_to_micron(wavenumber):
    """Convert wavenumber in cm\ :sup:`-1` to wavelength in micrometers.
    
    Args:
        wavenumber (int, float): Wavenumber value, units [cm\ :sup:`-1`].
    
    Returns:
        um (float): Wavelength value, units [\ :math:`\\mu`\ m].
    
    """
    um = (1 / wavenumber) * 1e4
    return um


def inv_cm_to_nm(wavenumber):
    """Convert wavenumber in cm\ :sup:`-1` to wavelength in nanometers.
    
    Args:
        wavenumber (int, float): Wavenumber value, units [cm\ :sup:`-1`].
    
    Returns:
        nm (float): Wavelength value, units [nm].
    
    """
    nm = (1 / wavenumber) * 1e7
    return nm


def micron_to_inv_cm(wavelength):
    """Convert wavelength in :math:`\\mu`\ m to wavenumber in cm\ :sup:`-1`.
    
    Args:
        wavelength (float): Wavelength value, units [\ :math:`\\mu`\ m].
        
    
    Returns:
        wvnm (int, float): Wavenumber value, units [cm\ :sup:`-1`].
    
    """
    wvnm = 1 / wavelength * 1e4
    return wvnm


def nm_to_inv_cm(wavelength):
    """Convert wavelength in nm to wavenumber in cm\ :sup:`-1`.
    
    Args:
        wavelength (float): Wavelength value, units [nm].
        
    
    Returns:
        wvnm (int, float): Wavenumber value, units [cm\ :sup:`-1`].
    
    """
    wvnm = 1 / wavelength * 1e7
    return wvnm
