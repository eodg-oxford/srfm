"""
Name: utilities
Parent package: srfm
Author: Antonin Knizek
Contributors: 
Date: 18 February 2025
Purpose: Provides functions for srfm that do not fall in any other category, incl.
decorator functions, memory-safe declarations, some physical formulas, etc.
""" 
import numpy as np
from . import units
import warnings
import psutil
import time
from numba import njit


def closest(lst_lon, lon, lst_lat, lat):
    """find index of closest pixel to lon and lat from IASI L1c data"""
    lonlat = np.column_stack((lst_lon, lst_lat))
    Tree = KDTree(lonlat)
    dd, ii = Tree.query([lon, lat])
    return ii


def convert_spectral_radiance_to_bbt(B, v):
    """Converts intensity (spectral radiance) to brightness temperature.
    B(T,v) (intens, spectral radiance) = [ W m-2 sr-1 cm ] == [ kg s-3 sr-1 cm ]
    note: the 'cm' comes from the intensity being per cm-1, i.e. 1/cm-1
    T - temperature, [K]
    h - Planck constant, 6.62607015e-34 kg m2 s-1
    c - speed of light
    kb - boltzmann constant
    v - wavenumber, [cm-1]

    T = h * c / kb * v * 1 / ( ln( 1 + 2 * h * c**2 * v**3 / B(T,v) ) )
    this expression was derived from:
    Dudhia 2017 (10.1016/j.jqsrt.2016.06.018), eq. 25 (and equally eq. 23)
    """
    h = 6.62607015e-30  # kg cm2 s-1
    c = 2.99792458e10  # cm s-1
    kb = 1.380649e-19  # kg cm2 s-2 K-1

    return h * c / kb * v / (np.log(1 + 2 * h * c**2 * v**3 / B))


def calc_tot_Rayleigh_opt_depth(ps, l):
    """Calculates total atmospheric Rayleigh optical depth
    Formula from: (according to Don, put ref here)
    ps - surface pressure in hPa/mbar
    l - wavelength in micrometers
    """
    tau = (ps / 1013.0) / (117.03 * l**4 - 1.316 * l**2)
    return tau


def calc_layer_opt_thick_Rayleigh(pl, pu, ps, tau0):
    """Calculates one layer optical thickness (optical depth) from Rayleigh 
    scattering.
    Formula from: (according to Don, put ref here)
    Layer defined by upper pressure pu and lower pressure pl.
    ps - surface pressure
    tau0 - total atmospheric Rayleigh optical depth
    formula:
    delta tau_Rayleigh(wavelength,pu,pl) = tau0*(pl-pu)/ps
    """
    tau = tau0 * (pl - pu) / ps
    return tau
    
def calc_Rayleigh_opt_depths(ps, pl, pu, l):
    """Calculates optical depths from Rayleigh scattering.
    Uses formulas:
    note: once Don provides reference, rename this function.
    ps - surface pressure
    pl - layer lower boundary pressure
    pu - layer upper boundary pressure
    l - wavenumber, [cm-1]
    ##
    tau_r - rayleigh optical depth
    tau0 - total atmospheric Rayleigh optical depth
    """
    # calculate total atmospheric Rayleigh optical depth at this wavenumber
    tau0 = calc_tot_Rayleigh_opt_depth(ps=ps, l=units.inv_cm_to_micron(l))

    # calculate layer Rayleigh optical thickness
    tau_r = calc_layer_opt_thick_Rayleigh(ps=ps, pl=pl, pu=pu, tau0=tau0)
    # convert from pandas series to numpy nd.array
    tau_r = tau_r.to_numpy()
    
    return tau_r

def line_break_str(txt, chars, delim, indent=0):
    r"""This function adds line ends ('\n') at desired places in a string.
    The original intention of this function is to ensure that no string
    is longer than {chars} characters for the RFM driver table, whose line length
    is limited to 200 characters.
    The function takes in a string txt, and finds the nearest lower occurence of 
    the required delimiter.
    It then adds {indent} spaces and a line end character, so that if written to 
    a txt file, the lines do not exceed {chars} characters.  
    """
    # check input value format
    if not isinstance(txt, str):
        raise TypeError("txt must be a string.")
    if not isinstance(chars, int):
        raise TypeError("chars must be type int.")
    if not isinstance(indent, int):
        raise TypeError("indent must be type int")
    if not isinstance(delim, str):
        raise TypeError("delimiter must be a string.")
    
    orig_txt = txt
    lines = []
    while len(txt) > chars:
        occur = (txt.rfind(delim))
        if occur > -1:
           line = txt[0:(occur+len(delim)-1)]
           lines.append(line)
           txt = txt[(occur+len(delim)):]
        elif occur == -1:
            msg = """Could not split the string - strings longer than the
            required chars limit without the occurence of the required delimiter
            detected. RFM line character limit is 200 characters."""
            warnings.warn(msg)
            return orig_txt
    lines.append(txt)
    ind = " "*indent
    
    return f"{ind}\n".join(lines)

def memory_safe_np_zeros_2d(constraints=None, pct=99, max_sec_dim = None):
    """Initialzes a numpy.zeros 2D array of maximum allowed size so as not to overflow 
    system RAM.
    constraints = list/tuple/1D-np.array with constraints on mimimum required shape
    pct = maximum allowed percentage of available RAM to be used, defualt 99%
    max_sec_dim = maximum allowed size of the second dimension
    """
    
    # check inputs:
    if not isinstance(pct,(int,float,type(None))):
        raise TypeError("Parameter pct must be int, float or None")
    if not isinstance(constraints, (list, tuple, np.ndarray)):
        raiseTypeError("Parameter constraints must be list, tuple or np.ndarray")
    if len(constraints) > 2:
        raise ValueError("""Number of constraints exceeds the number of required 
                        dimensions.""")
    if not isinstance(max_sec_dim, int):
        raise TypeError("Parameter max_sec_dim must be an array.")
    
    
    #get machine ram data
    ram = psutil.virtual_memory()
    
    # determine the size in bytes of one array item, should be 8, but in case this 
    # changes in the future, it's determined here from runtime
    isize = np.zeros(1).itemsize
    
    # calcualte the maximum number of elements the array can have
    max_array_items = ram.available / isize
    
    # determine shape of array to be declared and declare
    if constraints == None:
        prime_factors = find_prime_factors(max_array_ram)
        return np.zeros(
                (int(max(prime_factors)), int(max_array_items / max(prime_factors))
            ),
            dtype=np.float64
        )
    
    elif len(constraints) == 1:
        # floor division to determine second dimension
        second_dim = int(max_array_items // constraints[0])
        
        # reduce second dimension size if max value given
        if second_dim > max_sec_dim:
            second_dim = max_sec_dim
        return np.zeros(
                   (constraints[0], second_dim),
                   dtype=np.float64
               )

    elif len(constraints) == 2:
        if (constraints[0]*constraints[1] > max_array_items):
            raise ValueError("Requested array too large, not enough memory available.")
        else:
            return np.zeros(
                        (constraints[0],constraints[1]),
                        dtype=np.float64
                )
    

def find_prime_factors(num):
    """Find prime factors of a number
    inputs:
        num - number
    outputs:
        factors - list of prime factors
    """
    factors = []
    factor = 2
    while num >= 2:
        if (num % factor == 0):
            factors.append(factor)
            num = num / factor
        else:
            factor += 1
    return factors
    
def calc_layer_extent(a,t):
    """Calculates the upper and lower level of a layer (e.g. aerosol layer)
    inputs:
        a - altitude of the center of the layer
        t - layer thickness
    outputs:
        u - upper boundary altitude
        l - lower boundary altitude
    Input and output units match.
    """
    
    if not isinstance(a, (int, float)):
        raise TypeError("Altitude must be int or float.")
    if not isinstance(t, (int, float)):
        raise TypeError("Thickness must be int or float.")
    if a < 0:
        warnings.warn("""Negative altitude doesn't make sense on rocky planets, 
        are you sure?""")
    
    u = a + (t / 2)
    l = a - (t / 2)
    
    if l < 0:
        warnings.warn("""Layer is so thick that lower boundary is negative. This doesn't
        make sense on rocky planets, are the inputs correct?""")
    
    return u, l

def calc_layer_bounds(u,l):
    """Calculate atmospheric layer center altitude and vertical extent.
    inputs:
        u - layer upper boundary altitude
        l  - layer lower boundary altitude
    outputs:
        a - layer center altitude
        t - layer vertical extent (thickness)
    
    """
    if not isinstance(u, (int,float)):
        raise TypeError("Layer upper boundary altitude must be int or float.")
    if not isinstance(l, (int,float)):
        raise TypeError("Layer lower boundary altitude must be int or float.")
    if u <= l:
        raise ValueError("""Layer upper boundary altitude must
                          be > than lower boundary altitude.""")
    t = u - l
    a = l + t/2
    
    return a, t
    
def add_lyr(old_lev,u,l):
    """Adds a layer to an existing atmoshperic level structure.
    Deletes any levels from the existing structure that would fall within the new layer.
    inputs:
        old_lev - list with current user-desired layers
        u, l - new layer upper and lower boundary
    outputs:
         lev - new atmospheric level structure
         ilyr - particle layer index, useful when atmosphere specified such that layer 0
         is at the surface/lowerst altitude (e.g. in RFM)
         inv_ilyr - particle layer index from an inverted profile, useful when 
         atmosphere specified such that layer 0 is at the TOA (e.g. in DISORT).
    """
    if not isinstance(old_lev, (list, np.ndarray)):
        raise TypeError("old_lev must be a list or a 1D np.ndarray")
    if isinstance(old_lev, np.ndarray):
        if old_lev.ndim != 1:
            raise ValueError("if old_lev is a np.ndarray, it must be 1D.")
        elif old_lev.ndim == 1:
            try:
                old_lev = old_lev.tolist()
            except:
                raise RuntimeError("Could not convert lev to list.")
    if not isinstance(u, (int,float)):
        raise TypeError("Parameter 'u' must be int or float.")
    if not isinstance(l, (int,float)):
        raise TypeError("Parameter 'l' must be int or float.")

    to_remove = [i for i in old_lev if i <= u and i >= l]
    lev = [float(i) for i in old_lev if i not in to_remove]
    lev = lev + [float(u),float(l)]
    lev.sort()
    ilyr = lev.index(l)
    inv_ilyr = lev[::-1].index(u)
    return lev, ilyr, inv_ilyr

def calc_tot_dtauc(tau_g,tau_R,tau_p):
    """Calculate total optical depth of model layers, delta tau (dtau).
    dtau = tau_g + tau_R + tau_p
    dtau = layer optical depth
    tau_g = layer optical depth from gas absorption
    tau_R = layer optical depth from Rayleigh scattering
    tau_p = layer optical depth from particle scattering
    """
    def check_convert_dtype(obj):
        """Inner function that checks the input data types and tries to
        convert them the a suitable dtype.
        """
        if isinstance(obj, np.ndarray):
            pass
        elif isinstance(obj, pd.Series):
            obj = obj.to_numpy()
        elif isinstance(obj, list):
            obj = np.asarray(obj)
        elif isinstance(obj, (int,float)):
            pass
        else:
            raise TypeError(f'inputs must be np.ndarrays, pd.Series, lists, ints or floats.')
        
        return obj
    
    tau_g = check_convert_dtype(tau_g)
    tau_R = check_convert_dtype(tau_R)
    tau_p = check_convert_dtype(tau_p)

    return tau_g + tau_R + tau_p
    
def show_runtime(func):
    """Wrapper function to time other functions."""
    def wrapper(*args,**kwargs):
        t_start = time.perf_counter()
        result = func(*args,**kwargs)
        t_end = time.perf_counter()
        elapsed = (t_end - t_start)
        print(f"Time taken to execute {func.__name__}: {elapsed:.6f} seconds.")
        return result
    return wrapper
    
def number_conc_from_mass_loading(l, rho, thick, r=None, d=None):
    """Calculates particle number concentration from particle column loading.
    Inputs:
        l - particle column loading, [g m-2]
        rho - particle density, [kg m-3]
            - uniform density assumed
            - can be either a value in [kg m-3] or one of the following strings:
                "pumice" - 950 kg m-3
                "glass" - 2400 kg m-3
                "mineral" - 3000 kg m-3
                "rock" - 2900 kg m-3
            - note that the densities are average densities from https://volcanoes.usgs.
            gov/volcanic_ash/density_hardness.html#:~:text=Volcanic%20Ash,-Ash%20Particl
            e%20Size&text=For%20example%2C%20700%2D1200%20kilograms,material%20if%20depo
            sited%20on%20water.
            who in their turn take their data from Shipley and Sarna-Wojcicki, 1982
            (and don't give details on this reference at all).
            using the strings should serve only an illustrative purpose and should not
            be relied on as the data in a given eruption may vary!
        thick - particle layer thickness, [km]
        r - particle radius, [um], either r or d has to be given
        d - particle diameter, [um], either r or d has to be given
    Outputs:
        n - particle number concentration, [particles cm-3]
    """
    
    rho_dict = {"pumice" : 950,
                "glass" : 2400,
                "mineral" : 3000,
                "rock" : 2900
            }
    
    for i in [l, thick]:
        if not isinstance(i,(float, int)):
            raise TypeError(f"l and thick must be ints or floats.")
    
    if not isinstance(rho, (int, float, str)):
        raise TypeError("rho must be an int, float or str")
    elif isinstance(rho, str) and rho not in rho_dict.keys():
        raise ValueError(f"If rho is a string, it must be one of {rho_dict.keys()}.")
    elif isinstance(rho, str) and rho in rho_dict.keys():
        rho = rho_dict[rho]
    
    if r == None and d == None:
        raise ValueError(f"Both r and d are None. One has to be given.")
    elif r != None and d != None:
        if d/2 == r:
            pass
        else:
            raise ValueError(f"""Conflicting r and d are given. 
            Providing one is sufficient.""")
    elif r == None and d != None:
        r = d/2
    elif r != None and d == None:
        pass
    
    #factor 1e6 comes from unit conversions
    return 3 / 4 * l / (rho * thick * np.pi * r**3) * 1e6

def mass_loading_from_number_conc(n, thick, rho, r=None, d=None):
    """Calculate particle mass loading from number concentration.
     Inputs:
        n - particle number concentration, [particles cm-3]
        rho - particle density, [kg m-3]
            - uniform density assumed
            - can be either a value in [kg m-3] or one of the following strings:
                "pumice" - 950 kg m-3
                "glass" - 2400 kg m-3
                "mineral" - 3000 kg m-3
                "rock" - 2900 kg m-3
            - note that the densities are average densities from https://volcanoes.usgs.
            gov/volcanic_ash/density_hardness.html#:~:text=Volcanic%20Ash,-Ash%20Particl
            e%20Size&text=For%20example%2C%20700%2D1200%20kilograms,material%20if%20depo
            sited%20on%20water.
            who in their turn take their data from Shipley and Sarna-Wojcicki, 1982
            (and don't give details on this reference at all).
            using the strings should serve only an illustrative purpose and should not
            be relied on as the data in a given eruption may vary!
        thick - particle layer thickness, [km]
        r - particle radius, [um], either r or d has to be given
        d - particle diameter, [um], either r or d has to be given
    Outputs:
        l - particle column loading, [g m-2]
    """
    
    rho_dict = {"pumice" : 950,
                "glass" : 2400,
                "mineral" : 3000,
                "rock" : 2900
            }
    
    for i in [l, thick]:
        if not isinstance(i,(float, int)):
            raise TypeError(f"l and thick must be ints or floats.")
    
    if not isinstance(rho, (int, float, str)):
        raise TypeError("rho must be an int, float or str")
    elif isinstance(rho, str) and rho not in rho_dict.keys():
        raise ValueError(f"If rho is a string, it must be one of {rho_dict.keys()}.")
    elif isinstance(rho, str) and rho in rho_dict.keys():
        rho = rho_dict[rho]
    
    if r == None and d == None:
        raise ValueError(f"Both r and d are None. One has to be given.")
    elif r != None and d != None:
        if d/2 == r:
            pass
        else:
            raise ValueError(f"""Conflicting r and d are given. 
            Providing one is sufficient.""")
    elif r == None and d != None:
        r = d/2
    elif r != None and d == None:
        pass

    return 4 / 3 * np.pi * r**3 * thick * rho * n * 1e-6

#@njit
def monotonic(x):
    """Check is list/1D array is monotonic and increasing or decreasing.
    input - x, input array/list, must be 1D
    returns values:
        0 - not monotonic
        1 - strictly increasing
        2 - strictly decreasing
    
    """
    if isinstance(x,list):
        x = np.array(x, dtype=np.float64)
        
    dx = np.diff(x)
    
    if np.all(dx<=0): # is decreasing
        return 2
    elif np.all(dx >= 0): # is increasing
        return 1
    else:
        return 0
    
def calc_grids(lo, hi, res, units):
    """Calculates spectral grids from a given lower and upper limit and resolution
    in given units.
    inputs:
        lo - lower limit
        hi - upper limit
        res - resolution
        units - units of input parameters, can be "cm-1", "nm", "um"
    outputs:
        wvnm - wavenumber grid, [cm-1]
        wvls - wavelength grid, [um]
    The output grids are regular in the input units.
    """    
    # perform checks
    if not isinstance(lo, (int, float)):
        raise TypeError("lo must be int or float.")
    if not isinstance(hi, (int, float)):
        raise TypeError("lo must be int or float.")
    if not isinstance(res, (int, float)):
        raise TypeError("lo must be int or float.")
    if not isinstance(units, str):
        raise TypeError("units must be a string.")
    if units not in ["cm-1", "um", "nm"]:
        raise ValueError("Accepted values for units are 'cm-1', 'nm', and 'um'.")    
    
    if units == "cm-1":
        if not lo < hi and not lo >= 0.001:
            raise ValueError("lo must satisfy 0.001 cm-1 <= lo < hi.")
        if hi > 50000:
            raise ValueError("hi must satisfy lo < hi <= 50,000.")    
    elif units == "um":
        if not lo < hi and not lo >= 0.2:
            raise ValueError("lo must satisfy 0.2 um <= lo < hi.")
        if hi > 1e7:
            raise ValueError("hi must satisfy lo < hi <= 1e7 um.")    
    elif units == "nm":
        if not lo < hi and not lo >= 200:
            raise ValueError("lo must satisfy 200 nm <= lo < hi.")
        if hi > 1e10:
            raise ValueError("hi must satisfy lo < hi <= 1e10 nm.")
    
    # calculate the grids
    if units == "cm-1":
        wvnm = np.linspace(lo, hi, int((hi-lo)/res+1))
        wvls = (1/wvnm)*1e4
    elif units == "um":
        wvls = np.linspace(lo, hi, int((hi-lo)/res+1))[::-1]
        wvnm = (1/wvls)*1e4
    elif units == "nm":
        lo, hi, res = lo*1e-3, hi*1e-3, res*1e-3
        wvls = np.linspace(lo, hi, int((hi-lo)/res+1))[::-1]
        wvnm = (1/wvls)*1e4
    return wvnm, wvls
    
