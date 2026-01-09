"""Provides functions for srfm that do not fall in any other category.

Contains e.g. decorator functions, memory-safe declarations,
some physical formulas, etc.

- Name: utilities
- Parent package: srfm
- Author: Antonin Knizek
- Contributors:
- Date: 18 February 2025
"""

import numpy as np
from . import units
import warnings
import psutil
import time
from numba import njit
from bisect import bisect
from functools import wraps
from importlib.resources import files, as_file
from scipy.signal import convolve


def closest(lst_lon, lon, lst_lat, lat):
    """Find index of closest pixel to lon and lat from IASI L1c data.

    Originally written to work with IASI data.
    IASI data are spectra which are taken at given longitude and latitude. This function
    takes in longitude and latitude from the user and finds an IASI spectrum closest to
    the query.

    Can be used in general to find a pair of values closest to an existing pair in
    two list/2D array.

    Note that this function exists because the longitude and latitude lists are
    attributes of a pixel and tue purpose is to find a closest pixel, i.e. it is not
    equal to finding a minimum in each of the two lists, but rather is about finding
    a pixel whose distance is minimum to the required (lon, lat).

    Args:
        lst_lon (array-like): List of longitudes.
        lon (int, float): User-specified longitude.
        lst_lat (array-like): List of latitudes.
        lat (int, float): User-specified latitudes.

    Returns:
        ii (int): Index of closest pixel. The longitude and latitude value of this pixel
            are simply lst_lon[ii] and lst_lat[ii], respectively.

    """
    lonlat = np.column_stack((lst_lon, lst_lat))
    Tree = KDTree(lonlat)
    dd, ii = Tree.query([lon, lat])

    return ii


def convert_spectral_radiance_to_bbt(B, v):
    """Converts intensity (spectral radiance) to brightness temperature.

    Calculated according to:

    .. math::

        T_{BB} = \\frac{h  c  \\nu}{k_b  \\ln{(1 +
        \\frac{2  h  c^2  \\nu^3}{B(T,\\nu)})}}

    where:
        - :math:`B(T,\\nu)`: Intensity, spectral radiance, units
          [W m\ :sup:`-2` sr\ :sup:`-1` cm ], equiv. to
          [ kg s\ :sup:`-3` sr\ :sup:`-1` cm ].
        - :math:`T_{B}`: Brightness temperature (equiv. blackbody temperature),
          units [K].
        - *h*: Planck constant,
          6.62607015 :math:`\\times 10^{-34}` kg m\ :sup:`2` s\ :sup:`-1`.
        - *c*: speed of light, 2.99792458 :math:`\\times 10^{10}` cm s\ :sup:`-1`.
        - :math:`k_b`: Boltzmann constant, 1.380649 :math:`\\times 10^{-19}` kg
          cm\ :sup:`2` s\ :sup:`-2` K\ :sup:`-1`.
        - :math:`\\nu`: Wavenumber, units [cm\ :sup:`-1`].

    This expression was derived from `Dudhia 2017`_, eq. 25 (and equally eq. 23).

    .. _Dudhia 2017: https://doi.org/10.1016/j.jqsrt.2016.06.018

    Args:
        B (int, float): Intensity, spectral radiance, units
            [ W m\ :sup:`-2` sr\ :sup:`-1` cm ]. Note that the 'cm' comes from the
            intensity being per cm\ :sup:`-1`, i.e. :math:`\\frac{1}{cm^{-1}}`.
        v (int, float): Wavenumber, units [cm\ :sup:`-1`].

    Returns:
        T (float): Brightness temperature (equivalent black body temperature),
            units [K].

    """
    h = 6.62607015e-30  # kg cm2 s-1
    c = 2.99792458e10  # cm s-1
    kb = 1.380649e-19  # kg cm2 s-2 K-1

    T = h * c / kb * v / (np.log(1 + 2 * h * c**2 * v**3 / B))

    return T


def calc_tot_Rayleigh_opt_depth(ps, l):
    """Calculates total atmospheric Rayleigh optical depth.

    That is the optical Rayleight optical depth of the whole atmosphere.
    Formula from: (according to Don, put ref here)

    Args:
        ps (int, float): Surface pressure in hPa/mbar.
        l (int, float): Wavelength in micrometers.

    Returns:
        tau (float): Total Rayleigh optical depth.

    """
    tau = (ps / 1013.0) / (117.03 * l**4 - 1.316 * l**2)

    return tau


def calc_layer_opt_thick_Rayleigh(pl, pu, ps, tau0):
    """Calculates one layer optical thickness from Rayleigh scattering.

    Formula from to Don, put ref here.
    The calculation formula is

    .. math::

        \\begin{eqnarray}
            \\Delta \\tau_R(\\lambda,p_u,p_l) = \\tau_0 \\frac{(p_l-p_u)}{p_s}
        \\end{eqnarray}

    where:
        - :math:`\\Delta \\tau_R(\\lambda,p_u,p_l)` is the layer optical depth as
          a function of wavelength and layer upper and lower boundary pressures.
        - :math:`\\tau_0` is the total atmospheric Rayleight optical depth.
        - :math:`p_u, p_l` are the layer upper and lower boundary pressures.
        - :math:`p_s` is the surface pressure.

    Args:
        pl (int, float, array-like): Layer bottom boundary pressure.
        pu (int, float, array-like): Layer upper boundary pressure.
        ps (int, float): Surface pressure.
        tau0 (int, float): total atmospheric Rayleigh optical depth.

    Returns:
        tau (float): Layer Rayleigh optical depth.

    """
    tau = tau0 * (pl - pu) / ps

    return tau


def calc_Rayleigh_opt_depths(ps, pl, pu, l):
    """Calculates optical depths from Rayleigh scattering.

    Args:
        ps (int, float): Surface pressure.
        pl (int, float): Layer lower boundary pressure.
        pu (int, float): Layer upper boundary pressure.
        l (int, float): Wavenumber, [cm\ :sup:`-1`].

    Returns:
        tau_r (float): Rayleigh optical depth.
    """
    # calculate total atmospheric Rayleigh optical depth at this wavenumber
    tau0 = calc_tot_Rayleigh_opt_depth(ps=ps, l=units.inv_cm_to_micron(l))

    # calculate layer Rayleigh optical thickness
    tau_r = calc_layer_opt_thick_Rayleigh(ps=ps, pl=pl, pu=pu, tau0=tau0)
    # convert from pandas series to numpy nd.array
    tau_r = tau_r.to_numpy()

    return tau_r


def line_break_str(txt, chars, delim, indent=0):
    r"""This function adds line breaks at desired places in a string.

    The original intention of this function is to ensure that no string
    is longer than {chars} characters for the RFM driver table, whose line length
    is limited to 200 characters.
    The function takes in a string txt, and finds the nearest lower occurence of
    the required delimiter.
    It then adds *indent* spaces and a line end character, so that if written to
    a txt file, the lines do not exceed *chars* characters.

    Args:
        txt (str): String to be split
        chars (int): Number of character allowed per one line.
        delim (str): Character to split the string at.
        indent (int): The amount of spaces placed at the beginning of each split section.
            Default is 0.

    Returns:
        fin (str): New string with line breaks inserted in appropriate places.

    Raises:
        TypeError: Raised when inputs are incorrect datatypes.

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
        occur = txt.rfind(delim)
        if occur > -1:
            line = txt[0 : (occur + len(delim) - 1)]
            lines.append(line)
            txt = txt[(occur + len(delim)) :]
        elif occur == -1:
            msg = """Could not split the string - strings longer than the
            required chars limit without the occurence of the required delimiter
            detected. RFM line character limit is 200 characters."""
            warnings.warn(msg)
            return orig_txt
    lines.append(txt)
    ind = " " * indent

    fin = f"{ind}\n".join(lines)

    return fin


def memory_safe_np_zeros_2d(constraints=None, pct=99, max_sec_dim=10000):
    """Initialzes a numpy.zeros 2D array os size that does not overflow RAM.

    Creates a numpy.zeros array of maximum allowed size so as not to overflow
    system RAM.
    Note that by default, a np.zeros array has dtype np.float64. Therefore each
    element of the array takes up 8 bytes. The size of the array in the memory is
    simply the number of elements times the 8 bytes.

    Args:
        constraints (array-like): Constraints on mimimum required shape. Default is
            None.
        pct (int, float): Maximum allowed percentage of available RAM to be used,
            Defualt is 99%.
        max_sec_dim (int) Maximum allowed size of the second dimension. Default is
            10,000.

    Returns:
        arr (np.ndarray): Numpy.zeros array of appropriate shape.

    Raises:
        TypeError: Raised when inputs are incorrect dtype.
        ValueError: Raised when the number of constraints exceeds two (the amount of
            output array dimensions.
        RuntimeError: Raised when the size of array that would satisfy the required
            constraints would exceed the amount of available memory.

    Todo:
        Generalize this function, so that it takes the required number of dimensions as
        a parameter and can return ND array rather than just 2D array.

    """

    # check inputs:
    if not isinstance(pct, (int, float, type(None))):
        raise TypeError("Parameter pct must be int, float or None")
    if not isinstance(constraints, (list, tuple, np.ndarray)):
        raise TypeError("Parameter constraints must be list, tuple or np.ndarray")
    if len(constraints) > 2:
        raise ValueError(
            """Number of constraints exceeds the number of required 
                        dimensions."""
        )
    if not isinstance(max_sec_dim, int):
        raise TypeError("Parameter max_sec_dim must be an array.")

    # get machine ram data
    ram = psutil.virtual_memory()

    # determine the size in bytes of one array item, should be 8, but in case this
    # changes in the future, it's determined here from runtime
    isize = np.zeros(1).itemsize

    # calcualte the maximum number of elements the array can have
    max_array_items = ram.available / isize

    # determine shape of array to be declared and declare
    if constraints is None:
        prime_factors = find_prime_factors(max_array_ram)
        return np.zeros(
            (int(max(prime_factors)), int(max_array_items / max(prime_factors))),
            dtype=np.float64,
        )

    elif len(constraints) == 1:
        # floor division to determine second dimension
        second_dim = int(max_array_items // constraints[0])

        # reduce second dimension size if max value given
        if second_dim > max_sec_dim:
            second_dim = max_sec_dim
        return np.zeros((constraints[0], second_dim), dtype=np.float64)

    elif len(constraints) == 2:
        if constraints[0] * constraints[1] > max_array_items:
            raise RuntimeError(
                """Requested array too large, not enough memory
                available."""
            )
        else:
            return np.zeros((constraints[0], constraints[1]), dtype=np.float64)


def find_prime_factors(num):
    """Find prime factors of a number.

    Args:
        num (int): Input number.

    Returns:
        factors (list): List of prime factors of num.

    """
    factors = []
    factor = 2
    while num >= 2:
        if num % factor == 0:
            factors.append(factor)
            num = num / factor
        else:
            factor += 1

    return factors


def calc_layer_extent(a, t):
    """Calculates the upper and lower boundary of a layer (e.g. aerosol layer).

    Input and output units match.

    Args:
        a (int, float): Altitude of the center of the layer.
        t (int, float): Layer thickness.

    Return:
        u (int, float): Upper boundary altitude.
        l (int, float): Lower boundary altitude.

    Raises:
        TypeError: Raised if inputs are incorrect dtype.

    """

    if not isinstance(a, (int, float)):
        raise TypeError("Altitude must be int or float.")
    if not isinstance(t, (int, float)):
        raise TypeError("Thickness must be int or float.")
    if a < 0:
        warnings.warn(
            """Negative altitude doesn't make sense on rocky planets, 
        are you sure?"""
        )

    u = round(a + (t / 2), 3)
    l = round(a - (t / 2), 3)

    if l < 0:
        warnings.warn(
            """Layer is so thick that lower boundary is negative. This doesn't
        make sense on rocky planets, are the inputs correct?"""
        )

    return u, l


def calc_layer_bounds(u, l):
    """Calculate atmospheric layer center altitude and vertical extent.

    Args:
        u (int, float): Upper boundary altitude.
        l (int, float): Lower boundary altitude.

    Returns:
        a (int, float): Altitude of the center of the layer.
        t (int, float): Layer thickness.

    Raises:
        TypeError: Raised if inputs are incorrect dtype.
        ValueError: Raised when input lower boundary is greated than the upper boundary.

    """
    if not isinstance(u, (int, float)):
        raise TypeError("Layer upper boundary altitude must be int or float.")
    if not isinstance(l, (int, float)):
        raise TypeError("Layer lower boundary altitude must be int or float.")
    if u <= l:
        raise ValueError(
            """Layer upper boundary altitude must
                          be > than lower boundary altitude."""
        )
    t = round(u - l, 3)
    a = round(l + t / 2, 3)

    return a, t


def add_lyr_from_Layer(lev, track_lev, new_lyr):
    """Adds a layer to an existing atmoshperic level structure.

    Used to add layer properties from the *srfm.layer.Layer()* object
    Deletes any levels from the existing structure that would fall within the new layer.
    Updates a tracking array.

    Args:
        lev (list): List of current layers.
        track_lev (list): Helper list to track inserted levels, None for all levels
            except those inserted from a scattering layer.
        new_lyr (obj): ``srfm.Layer.MieLayer`` object, which contains upper and lower
            layer boundary altitudes, alt_upp and alt_low.

    Returns:
         lev (list): New atmospheric level structure.
         track_lev (list): Modified tracker list for level structure. None for all user
            specified levels, layer upper and lower boundary specifiers
            (``MieLayer.name``) for the new inserted levels.

    """

    to_remove_idcs = [
        lev.index(i) for i in lev if i <= new_lyr.alt_upp and i >= new_lyr.alt_low
    ]
    for ii in to_remove_idcs[::-1]:
        del track_lev[ii]
        del lev[ii]

    idx_low = bisect(lev, new_lyr.alt_low)
    lev.insert(idx_low, new_lyr.alt_low)
    track_lev.insert(idx_low, new_lyr.name)

    idx_upp = bisect(lev, new_lyr.alt_upp)
    lev.insert(idx_upp, new_lyr.alt_upp)
    track_lev.insert(idx_upp, new_lyr.name)

    return lev, track_lev


def calc_tot_dtauc(tau_g, tau_R, tau_p):
    """Calculate total optical depth of model layers.

    Formula is:

    .. math::

        \\Delta \\tau = \\tau_g + \\tau_R + \\tau_p

    where:
        - :math:`\\Delta \\tau`: Layer optical depth.
        - :math:`\\tau_g`: Layer optical depth from gas absorption.
        - :math:`\\tau_R`: Layer optical depth from Rayleigh scattering.
        - :math:`\\tau_p`: Layer optical depth from particle scattering.

    Args:
        tau_g (int, float, array-like): Layer optical depth from gas absorption.
        tau_R (int, float, array-like): Layer optical depth from Rayleigh scattering.
        tau_p (int, float, array-like): Layer optical depth from particle scattering.

    Returns:
        dtau (float, array-like): Layer optical depths.

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
        elif isinstance(obj, (int, float)):
            pass
        else:
            raise TypeError(
                """Inputs must be np.ndarrays, pd.Series, lists, ints 
                or floats."""
            )

        return obj

    tau_g = check_convert_dtype(tau_g)
    tau_R = check_convert_dtype(tau_R)
    tau_p = check_convert_dtype(tau_p)

    dtau = tau_g + tau_R + tau_p

    return dtau


def show_runtime(func):
    """Decorator function to time other functions."""

    @wraps(func)
    def wrapper(*args, **kwargs):
        t_start = time.perf_counter()
        result = func(*args, **kwargs)
        t_end = time.perf_counter()
        elapsed = t_end - t_start
        print(f"Time taken to execute {func.__name__}: {elapsed:.6f} seconds.")

        return result

    return wrapper


def number_conc_from_mass_loading(l, rho, thick, s, dist_type, r=None, d=None):
    """Calculates particle number concentration from particle column loading.

    Either *r* or *d* has to be given.

    Args:
        l (int, float): Particle column mass loading, units [g m\ :sup:`-2`].
        rho (int, float, str): Particle density, units [kg m\ :sup:`-3`]. Uniform density
            is assumed. Can be either a value in [kg m\ :sup:`-3`] or one of the
            following strings:

                - "pumice" - 950 kg m\ :sup:`-3`
                - "glass" - 2400 kg m\ :sup:`-3`
                - "mineral" - 3000 kg m\ :sup:`-3`
                - "rock" - 2900 kg m\ :sup:`-3`

            Note that the densities are average densities from the `USGS VAI`_.

            .. _USGS VAI: https://volcanoes.usgs.gov/volcanic_ash/density_hardness.html

            They in turn cite `Shipley and Sarna_Wojcicki (1982)`_.

            .. _Shipley and Sarna_Wojcicki (1982): https://pubs.usgs.gov/mf/1983/1435/report.pdf

            The same information is also repeated in `Wilson et al. (2012)`_.

            .. _Wilson et al. (2012): http://dx.doi.org/10.1016/j.pce.2011.06.006

            Using the strings should serve only an illustrative purpose and should not
            be relied on as the data in a given eruption may vary!
        thick (int, float): Particle layer thickness, units [km].
        s (int, float): Lognormal distribution spread. Default is None.
        dist_type (str): Size distribution type. If None, then a monodisperse distribution
            is implied and parameter s is not used. Other permitted value is "lognormal".
            Default is None.
        r (int, float): Particle radius, units [\ :math:`\\mu`\ m]. Default is None.
        d (int, float): Particle diameter, units [\ :math:`\\mu`\ m], Default is None.

    Returns:
        n (float): Particle number concentration, units [particles cm\ :sup:`-3`].

    Raises:
        TypeError: Raised if inputs are incorrect format.
        ValueError: Raised if neither r or d are given.
        ValueError: Raised if both r and d are given, but d != 2*r.
        ValueError: Raise if distribution type isn't recognized or implemented.

    Todo:
        Add other size distributions as an option.

    """

    rho_dict = {"pumice": 950, "glass": 2400, "mineral": 3000, "rock": 2900}

    for i in [l, thick]:
        if not isinstance(i, (float, int)):
            raise TypeError(f"l and thick must be ints or floats.")

    if not isinstance(rho, (int, float, str)):
        raise TypeError("rho must be an int, float or str")
    elif isinstance(rho, str) and rho not in rho_dict.keys():
        raise ValueError(f"If rho is a string, it must be one of {rho_dict.keys()}.")
    elif isinstance(rho, str) and rho in rho_dict.keys():
        rho = rho_dict[rho]

    if dist_type is not None and not isinstance(dist_type, str):
        raise TypeError("dist_type must be None or str.")

    if s is not None and not isinstance(s, (int, float)):
        raise TypeError("s must be None, int or float.")

    if r is None and d is None:
        raise ValueError(f"Both r and d are None. One has to be given.")
    elif r is not None and d is not None:
        if d / 2 == r:
            pass
        else:
            raise ValueError(
                f"""Conflicting r and d are given. 
            Providing one is sufficient."""
            )
    elif r is None and d is not None:
        r = d / 2
    elif r is not None and d is None:
        pass

    # factor 1e6 comes from unit conversions
    if dist_type == None:
        n = 3 / 4 * l / (rho * thick * np.pi * r**3) * 1e6
    elif dist_type == "log_normal":
        n = (
            3
            / 4
            * l
            / (rho * thick * np.pi * r**3 * np.exp(9 / 2 * (np.log(s)) ** 2))
            * 1e6
        )
    elif dist_type == "normal":
        raise ValueError("Normal distribution hasn't been implemented yet. Sorry. :(")
    else:
        raise ValueError("Unrecognized or unsupported distribution type.")

    return n


def mass_loading_from_number_conc(n, thick, rho, s, dist_type, r=None, d=None):
    """Calculate particle mass loading from number concentration.

    Either *r* or *d* has to be given.

    Args:
        n (int, float): Particle number concentration, units [particles cm\ :sup:`-3`].
        rho (int, float, str): Particle density, units [kg m\ :sup:`-3`]. Uniform density
            is assumed. Can be either a value in [kg m\ :sup:`-3`] or one of the
            following strings:

                - "pumice" - 950 kg m\ :sup:`-3`
                - "glass" - 2400 kg m\ :sup:`-3`
                - "mineral" - 3000 kg m\ :sup:`-3`
                - "rock" - 2900 kg m\ :sup:`-3`

            Note that the densities are average densities from the `USGS VAI`_.

            .. _USGS VAI: https://volcanoes.usgs.gov/volcanic_ash/density_hardness.html

            They in turn cite `Shipley and Sarna_Wojcicki (1982)`_.

            .. _Shipley and Sarna_Wojcicki (1982): https://pubs.usgs.gov/mf/1983/1435/report.pdf

            The same information is also repeated in `Wilson et al. (2012)`_.

            .. _Wilson et al. (2012): http://dx.doi.org/10.1016/j.pce.2011.06.006

            Using the strings should serve only an illustrative purpose and should not
            be relied on as the data in a given eruption may vary!
        thick (int, float): Particle layer thickness, units [km].
        s (int, float): Lognormal distribution spread. Default is None.
        dist_type (str): Size distribution type. If None, then a monodisperse distribution
            is implied and parameter s is not used. Other permitted value is "lognormal".
            Default is None.
        r (int, float): Particle radius, units [\ :math:`\\mu`\ m]. Default is None.
        d (int, float): Particle diameter, units [\ :math:`\\mu`\ m], Default is None.

    Returns:
        l (int, float): Particle column mass loading, units [g m\ :sup:`-2`].

    Raises:
        TypeError: Raised if inputs are incorrect format.
        ValueError: Raised if neither r or d are given.
        ValueError: Raised if both r and d are given, but d != 2*r.
        ValueError: Raise if distribution type isn't recognized or implemented.

    Todo:
        Add other size distributions as an option.

    """

    rho_dict = {"pumice": 950, "glass": 2400, "mineral": 3000, "rock": 2900}

    for i in [n, thick]:
        if not isinstance(i, (float, int)):
            raise TypeError(f"l and thick must be ints or floats.")

    if not isinstance(rho, (int, float, str)):
        raise TypeError("rho must be an int, float or str")
    elif isinstance(rho, str) and rho not in rho_dict.keys():
        raise ValueError(f"If rho is a string, it must be one of {rho_dict.keys()}.")
    elif isinstance(rho, str) and rho in rho_dict.keys():
        rho = rho_dict[rho]

    if dist_type is not None and not isinstance(dist_type, str):
        raise TypeError("dist_type must be None or str.")

    if s is not None and not isinstance(s, (int, float)):
        raise TypeError("s must be None, int or float.")

    if r is None and d is None:
        raise ValueError(f"Both r and d are None. One has to be given.")
    elif r is not None and d is not None:
        if d / 2 == r:
            pass
        else:
            raise ValueError(
                f"""Conflicting r and d are given. 
            Providing one is sufficient."""
            )
    elif r is None and d is not None:
        r = d / 2
    elif r is not None and d is None:
        pass

    if dist_type == None:
        l = 4 / 3 * np.pi * r**3 * thick * rho * n * 1e-6
    elif dist_type == "log_normal":
        l = (
            4
            / 3
            * np.pi
            * r**3
            * np.exp(9 / 2 * (np.log(s)) ** 2)
            * thick
            * rho
            * n
            * 1e-6
        )
    elif dist_type == "normal":
        raise ValueError("Normal distribution hasn't been implemented yet. Sorry. :(")
    else:
        raise ValueError("Unrecognized or unsupported distribution type.")

    return l


# @njit
def monotonic(x):
    """Check if list/1D array is monotonic and increasing or decreasing.

    Args:
        x (array-like): input array/list, if array, then must be 1D.

    Returns:
        val (int): Int value which signifies the result. Can be:
            - 0 - Input not monotonic.
            - 1 - Input strictly increasing.
            - 2 - Input strictly decreasing.

    """
    if isinstance(x, list):
        x = np.array(x, dtype=np.float64)

    dx = np.diff(x)

    if np.all(dx <= 0):  # is decreasing
        return 2
    elif np.all(dx >= 0):  # is increasing
        return 1
    else:
        return 0


def calc_grids(lo, hi, res, units):
    """Calculates spectral grids.

    Grids are calculated from a given lower and upper limit and resolution in given
    units.
    The output grids are regular in the input units.

    Args:
        lo (int, float): Grid lower limit.
        hi (int, float): Grid upper limit.
        res (int, float): Grid resolution.
        units (str): Units of input parameters, can be "cm\ :sup:`-1`", "nm", or
            "\ :math:`\\mu`\ m".

    Returns:
        wvnm (array): Wavenumber grid, units [cm\ :sup:`-1`].
        wvls (array): Wavelength grid, units [\ :math:`\\mu`\ m].

    Raises:
        TypeError: Raised if inputs are incorrect dtypes.
        ValueError: Raised if units is an unknown value.
        ValueError: Raised if lo and hi exceed prescribed limits. The limits are set
            to match the limits of RFM.

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
        npts = int(np.floor((hi - lo) / res)) + 1 # expected number of points in the grid
        wvnm = lo + np.arange(npts) * res
        wvls = (1 / wvnm) * 1e4
    elif units == "um":
        npts = int(np.floor((hi - lo) / res)) + 1 # expected number of points in the grid
        wvls = lo + np.arange(npts) * res
        wvnm = (1 / wvls) * 1e4
    elif units == "nm":
        lo, hi, res = lo * 1e-3, hi * 1e-3, res * 1e-3
        npts = int(np.floor((hi - lo) / res)) + 1 # expected number of points in the grid
        wvls = lo + np.arange(npts) * res
        wvnm = (1 / wvls) * 1e4
    return wvnm, wvls


def track_lev_to_track_lyr(lev):
    """Converts a levels tracking list to layer tracking list.

    In the levels tracking list, each element represents an atmospheric level and
    evaluates to *None* for all levels except the tracked level.
    The tracking list for layer follows the same logic.

    Example logic:
    (*L* is the layer name in place of the lower and upper layer boundary
    altitudes.)

    track_lev = [*None*, *None*, *L*, *L*, *None*, *None*]

    is converted to

    track_lyr = [*None*, *None*, *L*, *None*, *None*]

    Args:
        lev (list): Levels tracking list.

    Returns:
        lyr (list): Layers tracking array.

    """

    lyr = []
    for i in range(1, len(lev)):
        if lev[i] is not None and lev[i - 1] is not None and lev[i] == lev[i - 1]:
            lyr.append(lev[i])
        else:
            lyr.append(None)
    return lyr


def read_ils(filename):
    """Read instrument line shape file in RFM format.

    For description of the format, see `here`_.

    .. _here: https://eodg.atm.ox.ac.uk/RFM/sum/ilsfil.html

    Args:
        filename (str): Filename of the instrument line shape file in RFM format.

    Returns:
        x (array): Wavenumber/wavelength values.
        y (array): Intensity (line shape).
        lo (int): Lower wavenumber/wavelength (where the line shape starts).
        hi (int): Upper wavenumber/wavelength (where the line shape ends).

    """

    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    head = lines[3].split()
    npts = int(head[0])
    lo = float(head[1])
    res = float(head[2])
    hi = lo + (npts - 1) * res
    if len(head) > 3:
        wnorange = head[3]
    x = np.linspace(lo, hi, npts)
    y = np.array([float(x) for xs in [i.split() for i in lines[4:]] for x in xs])
    return x, y, lo, hi


def scale_solar_spectrum(spc, yday):
    """Scales solar spectrum to the correct year day number.

    Input solar spectrum is assumed for to be for 1 AU Sun-Earth distance. Since the
    Sun-Earth distance changes throughout the year, the spectrum needs to be scaled
    accordingly.
    Calculation formula from Don Grainger.

    Calculation formula:

    .. math::

        F(d) = \\left\\{1 - 0.0167086 \\cos\\left[\\frac{2\pi(d - 4)}
               {365.256363}\\right]\\right\\}^{-2} F_0

    where:
        - :math:`F(d)`: scaled flux
        - :math:`d`: year day number
        - :math:`F_0`: flux at 1 AU


    Args:
        spc (array): Solar spectrum, assumed valid for 1 AU Sun-Earth distance.
        yday (int, float): Year day number, Jan 1 is day 1.

    Returns:
        sc_spc (array): Input spectrum scaled to the input year day number.

    Raises:
        ValueError: Raised when 0 < yday < 366 is not True.

    """
    if yday < 0 or yday > 366:
        raise ValueError("Yday must satisfy 0 < yday < 366.")

    sc_spc = (1 - 0.0167086 * np.cos((2 * np.pi * (yday - 4)) / 365.256363)) ** -2 * spc

    return sc_spc


def get_altitude_prf(p, T, z0=0, M=28.97, g=9.81):
    """Calculates atmospheric altitude profile given pressure and temperature profiles.

    Manually integrates/sums the hydrostatic equation (using ideal gas law) to
    calculate atmospheric altitude profile given pressure and temperature profiles.

    Calculation derived from a combination of:

    ..math::

        p = z \\rho g

    and

    ..math::

        p V = \\frac{m R T}{M}

    where:
        - :math:`p`: atmospheric pressure (Pa)
        - :math:`z`: altitude (m)
        - :math:`\\rho`: air density
        - :math:`g`: gravity (:math:`m s^{-2}`)
        - :math:`M`: air molar mass (:math:`g mol^{-1}`)
        - :math:`R`: universal gas constant, :math:`8.314 J K mol^{-1}`
        - :math:`V`: gas volume
        - :math:`T`: tempeature (K)
        - :math:`m`: mass (g)

    The final calculation formula is then:

    ..math::

        z_{i+1} = - \\frac{\\ln(\\frac{p_{i+1}}{i}) R \\frac{T_i + T_{i+1}}{2}}{M g} + z_i

    Args:
        p (list): Atmospheric pressure profile, units [Pa].
            Expected sorted from lower to higher altitudes.
        T (list): Atmospheric temperature profile, units [K].
            Expected sorted from lower to higher altitudes.
        z0 (float): Surface altitude. Units [m]. Default is 0.
        M (float): Air molar mass. Units [:math:`g mol^{-1}`]. Default
            is :math:`28.97 g mol^{-1}`. Note that air molar mass is considered
            independent of altitude (is constant). THis formula therefore fails above
            the homopause.
        g (float): Gravity. Units :math:`m s^{-2}`. Default is :math:`9.81 m s^{-2}`.
            Note that gravity is independent of altitude in this calculation.

    Returns:
        alt (list): Atmospheric altitude profile. Units [km].

    """

    R = 8.314  # universal gas constant

    if isinstance(p, np.ndarray):
        print(
            """Supplied pressure profile is an array. Conversion to list will be 
        attempted."""
        )
        p = list(p)

    if isinstance(T, np.ndarray):
        print(
            """Supplied pressure profile is an array. Conversion to list will be 
        attempted."""
        )
        T = list(T)

    if len(p) != len(T):
        raise ValueError("Pressure and temperature are not of the same length.")

    if not all(x > y for x, y in zip(p, p[1:])):
        warnings.warn(
            """Pressure not strictly decreasing with altitude. Is this 
            intentional?"""
        )

    # set first altitude
    alt = [z0]

    for i in range(len(p) - 1):
        alt.append(
            (-(np.log(p[i + 1] / p[i]) * R * ((T[i] + T[i + 1]) / 2)) / (M * g))
            + alt[i]
        )

    return alt


def load_solar_spectrum_Gueymard20018():
    """Loads and converts solar spectrum from Gueymard 2018."""
    with as_file(files("srfm.data") / "Gueymard2018.sssi") as path:
        solar_spc_fl = np.loadtxt(path, skiprows=3)
        solar_spc = (
            solar_spc_fl[:, 1] * 10
        )  # conversion from mW cm-2 um-1 to W m-2 um-1
        solar_spc = (
            solar_spc / 1e4 * (solar_spc_fl[:, 0] ** 2)
        )  # conversion from W m-2 um-1 to W m-2 cm-1, wavenumbers are decreasing
        solar_spc_wvnm = 1 / (solar_spc_fl[:, 0] * 1e-4)  # wavenumbers are decreasing
    return solar_spc, solar_spc_wvnm

def convolve_spectrum(spc, x, ils):
    """Convolves a calculated spectrum with an instrument line shape.

    Note that the convolved spectrum suffers from boundary effects (scipy
    interpolate does zero padding of the data at the boundaries). Best avoided 
    by calculating your original spectra at a wider interval and then
    interpolating/truncating. For interpolation (best used to get spectra at 
    your simulated satellite grid).

    Args:
        spc (array-like): Spectrum to be convolved. Could be as short as single
            value.
        x (array-like): Spectral grid (wavenumbers/wavelengths/..). Must be
            regularly spaced.
        ils (str, ): Path to the instrument line shape. Make sure the units of 
            the ils and spectrum grids are identical. If the file extension is
            .atm, then the file is assumed to be and RFM-formatted .atm file
            (see `here`_).
            
            .. _here: https://eodg.atm.ox.ac.uk/RFM/sum/atmfil.html
            
            If the file is in any other format, then the file is read through
            numpy.loadtxt, so must conform to that format with the first column
            being the grid and the second column the intensity/radiance.
        
    Returns:
        spc_conv (array-like): Convolved spectrum.
    
    """
    
    # read instrument line shape
    if ils.endswith(".ils"):
        ils_x, ils_y, ils_lo, ils_hi = read_ils(ils)
    else:
        ils_all = np.loadtxt(ils)
        ils_x = ils_all[:,0]
        ils_y = ils_all[:,1]
        ils_lo = ils_x[0]
        ils_hi = ils_x[-1]

    # check if spectral grid is regular
    a = np.diff(x, n=2)  # calculate 2nd discrete difference
    a[a < 1e12] = 0  # remove small numbers (arising from computer precision limits)
    assert not np.all(a), """Wavenumber grid is not regular."""  # check is all values in a are 0 (0 evaluates to False)

    # determine resolution from model wavenumber grid
    num = len(x)
    lo = x[0]
    hi = x[-1]
    res = np.round((hi - lo) / num, decimals=8)  
    # this is inadvertedly introduces a limit
    # on the minimum resolution used in the code as 1e-8 cm-1, which should be
    # enough though, and also this may not be the numerically most stable way to go

    # generate new grid for ils
    npts = int(np.floor((ils_hi - ils_lo) / res)) + 1 # expected number of points in the grid
    new_x = ils_lo + np.arange(npts) * res

    # interpolate ils to new grid
    new_y = np.interp(new_x, ils_x, ils_y)

    # calculate sum of instrument line shape for normalization later
    norm = np.sum(new_y)
        
    # convolve
    spc_c = convolve(spc, new_y, mode="same") / norm

    return spc_c

