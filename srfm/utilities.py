import numpy as np
from . import units
import warnings
import psutil


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
    
    lines = []
    while len(txt) > chars:
        line = txt[0:txt.rfind(delim)+len(delim)+1]
        lines.append(line)
        txt = txt[txt.rfind(delim)+len(delim)+1:]
    ind = " "*indent
    
    if len(lines) == 0:
        warnings.warn("""No occurences of the required delimiter
             were found before the line break limit.
             String not split.""")
        return txt
    else:
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
            )
        )
    
    elif len(constraints) == 1:
        # floor division to detemrine second dimension
        second_dim = int(max_array_items // constraints[0])
        
        # reduce second dimension size if max value given
        if second_dim > max_sec_dim:
            second_dim = max_sec_dim
        return np.zeros(
                   (constraints[0],second_dim)
               )

    elif len(constraints) == 2:
        if (constraints[0]*constraints[1] > max_array_items):
            raise ValueError("Requested array too large, not enough memory available.")
        else:
            return np.zeros(
                        (constraints[0],constraints[1]
                    )
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
