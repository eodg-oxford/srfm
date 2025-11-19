"""Defines the Layer class and all its subclasses.

Used to contain scattering information for a single layer.
- Name: layer
- Parent package: srfm
- Author: Antonin Knizek
- Contributors:
- Date: 26 Mar 2025
"""

import numpy as np
from . import optical_properties as op
from . import utilities as utils
import warnings
from . import size_distribution as sz
from multiprocessing import Process, Manager
import matplotlib.pyplot as plt


class Layer:
    """Base superclass for an atmospheric layer."""

    def __init__(self, name=None, **parameters):
        self.name = name
        for key, val in parameters.items():
            self.key = val


class MieLayer(Layer):
    """Class that contains Mie scattering parameters, calculations and outputs.

    Is subclass of Layer.

    """

    def __init__(self, name=None, **parameters):
        super().__init__(name)
        for key, val in parameters.items():
            self.key = val

    def set_name(self, name):
        """Assign a name to the layer.

        Args:
            name (str): Layer name.

        """
        self.name = name

    def set_spc_lim(self, lo, hi):
        """Set spectral calculation grid lower and upper limits.

        Lo is the lower value and hi is the upper value, whatever the units. Units are
        then specified in set_spec_units, so this remains consistent with both
        wavenumbers and wavelengths.

        Args:
            lo (int, float): Lower wavenumber or wavelength.
            hi (int, float): Upper wavenumber or wavelength.

        """
        self.low_spc = lo
        self.upp_spc = hi

    def set_spec_units(self, units):
        """Set spectral calculation grid units.

        Args:
            units (str): Units. should be one of [\ :math:`\\mu`\ m, cm\ :sup:`-1`, nm].
                This is not enforced but currently the rest of this package is unable to
                handle any other options.

        """
        self.spec_units = units

    def set_res(self, res):
        """Set spectral calculation grid resolution.

        Args:
            res (int, float): Resolution of the spectral grid. Should be same units as
                values in set_spc_lim. This is not enforced, but is the only logically
                consistent option.

        """
        self.res = res

    def set_mass_loading(self, load):
        """Set particle mass loading.

        Args:
            load (int, float): Particle mass loading, units [g m\ :sup:`-2`].

        """
        self.mass_loading = load

    def set_n(self, n):
        """Set particle number concentration.

        Args:
            n (int, float): Particle number concentration, units [cm\ :sup:`-1`].

        """
        self.n = n

    def set_r(self, r):
        """Set mean particle radius.

        Args:
            r (int, float): Mean particle radius, units [\ :math:`\\mu`\ m]

        """
        self.r = r

    def set_s(self, s):
        """Set particle size distribution spread.

        Args:
            s (int, float): Particle size distribution spread.

        """
        self.s = s

    def set_rho(self, rho):
        """Set particle mean density.

        Args:
            rho (str, int, float): Particle mean density units [kg m\ :sup:`-3`\ ] or one
                of strings accepted by
                ``srfm.utilities.number_conc_from_mass_loading()``.

        """
        self.rho = rho

    def set_s_a_den(self, s_a_den):
        """Set surface area density of the particles.

        Args:
            s_a_den (int, float): Surface area density of the particles, units
                [\ :math:`\\mu`\ m\ :sup:`2` cm\ :sup:`-3`\ ].

        """
        self.s_a_den = s_a_den

    def set_v_den(self, v_den):
        """Set volume density of particles.

        Args:
            v_den (int, float): Particle volume density, units
                [\ :math:`\\mu`\ m\ :sup:`3` cm\ :sup:`-3`\ ].

        """
        self.v_den = v_den

    def set_dist_type(self, dist_type):
        """Set particle size distribution type.

        Args:
            dist_type (str): Particle size distribution type. Accepted values are
                *lognormal* and *normal*.

        """
        self.dist_type = dist_type

    def set_comp(self, comp):
        """Set particle composition type.

        Args:
            comp (str): One of accepted values by ``ARIA_module.get_ri_filepathname()``.

        """
        self.comp = comp

    def set_center_alt(self, center_alt):
        """Set average (center) particle layer altitude.

        Args:
            center_alt (int, float): Center particle layer altitude, units [km].

        """
        self.center_alt = center_alt

    def set_thick(self, thick):
        """Set layer thickness (vertical extent).

        Args:
            thick (int, float): Layer thickness (vertical extent), units [km].

        """
        self.thick = thick

    def set_alt_lim(self, alt_low, alt_upp):
        """Set lower and upper layer boundary (altitudes).

        Args:
            alt_low (int, float): Layer lower boundary altitude, units [km].
            alt_upp (int, float): Layer upper boundary altitude, units [km].

        """
        self.alt_upp = alt_upp
        self.alt_low = alt_low

    def set_radii(self, radii=200):
        """Set number of radii for size disitribution calculation.

        Args:
            radii (int): Number of radii in size distribution. Default is 200.

        """
        self.radii = radii

    def set_eta(self, eta=1e-6):
        """Set size distribution cut-off (eta value).

        Args:
            eta (int, float): Size distribution cut-off value. The size distribution is
                a function that technically spans the (-inf,+inf) size interval. The eta
                value is a value of the size distribution beyond whose corresponding
                size (radius) the distribution is truncated. Default is 1e-6.

        """
        self.eta = eta

    def set_phase_quad_N(self, phase_quad_N=181):
        """Set number of points for the scattering quadrature.

        The number of quadrature points is the number of scattering angles

        Args:
            phase_quad_N (int): Number of quadrature points. Default is 181.

        """
        self.phase_quad_N = phase_quad_N

    def set_phase_quad_type(self, phase_quad_type="L"):
        """Set quadrature type that will be used to calculate the phase fucntion.

        Args:
            phase_quad_type (str): Quadrature type for the phase function. Accetped
                values are values implemented in the ``srfm.quadrature`` module.
                Currently implemented (Apr 2025) are "L" (Lobatto quadrature rule),
                "G" (Gauss), "R" (Radau) and "T" (Trapezium). Default is L.

        """
        self.phase_quad_type = phase_quad_type

    def set_radii_quad_type(self, radii_quad_type="T"):
        """Set quadrature type for the calculation of the size distribution.

        Args:
            radii_quad_type (str): Quadrature type for the radii calculation function.
                Accetped values are values implemented in the ``srfm.quadrature``
                module. Currently implemented (Apr 2025) are "L" (Lobatto quadrature
                rule), "G" (Gauss), "R" (Radau) and "T" (Trapezium). Default is T.


        """
        self.radii_quad_type = radii_quad_type

    def set_leg_coeffs(self, leg_coeffs=True):
        """Toggle expansion of the phase fucntion into Legendre polynomials.

        Args:
            leg_coeffs (bool): If True, Legendre polynomial expansion of the phase
                function is calculated by ``srfm.optical_properties.ewp_hs()``. Default is
                True.

        """
        self.leg_coeffs = leg_coeffs

    def set_leg_coeffs_type(self, leg_coeffs_type="normalised"):
        """Set type of requested Legendre polynomial expansion coefficients.

        Args:
            leg_coeffs_type (str): Requested type of Legendre polynomial expansion
                coefficients. Accepted values are *regular* and *normalised*. Default
                is normalised.

        """
        self.leg_coeff_type = leg_coeffs_type

    def set_multiproccess(self, multiprocess):
        """Toggle multiprocessing.

        Args:
            multiprocess (bool): if True, then ``srfm.optical_properties.ewp_hs()`` is
                calculated in a separate subprocess. The output type of ewp_hs changes!!

        """
        self.multiprocess = multiprocess

    def set_input_from_dict(self, inp_dict):
        """Set all necessary input_parameters from an input dictionary.

        Args:
            inp_dict (dict): Input dictonary with layer properties. Must contain
                all necessary keys (cannot be incomplete). Required keys are:

                    - name
                    - low_spc
                    - upp_spc
                    - spec_units
                    - res
                    - mass_loading
                    - n
                    - r
                    - s
                    - rho
                    - s_a_den
                    - v_den
                    - dist_type
                    - comp
                    - center_alt
                    - thick
                    - alt_upp
                    - alt_low
                    - radii
                    - eta
                    - phase_quad_N
                    - phase_quad_type
                    - radii_quad_type
                    - leg_coeffs
                    - leg_coeffs_type
                    - multiprocess

                For explanation of each of those parameters please refer to the
                respective function which set them explicitly (set_{parameter name}).

        """
        self.name = inp_dict["name"]
        self.low_spc = inp_dict["low_spc"]
        self.upp_spc = inp_dict["upp_spc"]
        self.spec_units = inp_dict["spec_units"]
        self.res = inp_dict["res"]
        self.mass_loading = inp_dict["mass_loading"]
        self.n = inp_dict["n"]
        self.r = inp_dict["r"]
        self.s = inp_dict["s"]
        self.rho = inp_dict["rho"]
        self.s_a_den = inp_dict["s_a_den"]
        self.v_den = inp_dict["v_den"]
        self.dist_type = inp_dict["dist_type"]
        self.comp = inp_dict["comp"]
        self.center_alt = inp_dict["center_alt"]
        self.thick = inp_dict["thick"]
        self.alt_upp = inp_dict["alt_upp"]
        self.alt_low = inp_dict["alt_low"]
        self.radii = inp_dict["radii"]
        self.eta = inp_dict["eta"]
        self.phase_quad_N = inp_dict["phase_quad_N"]
        self.phase_quad_type = inp_dict["phase_quad_type"]
        self.radii_quad_type = inp_dict["radii_quad_type"]
        self.leg_coeffs = inp_dict["leg_coeffs"]
        self.leg_coeffs_type = inp_dict["leg_coeffs_type"]
        self.multiprocess = inp_dict["multiprocess"]

    def test_complete_input_format(self):
        """Checks the input for correct format before running any calculation.

        Tests all input variables. Creates a boolean value called passmark, which is
        initially False and if the test is passed, is changed to True and returned.

        Returns:
            passmark (bool): If True, input has passed the format test.

        Raises:
            TypeError: Raised when test fails for each parameter separately.

        """
        passmark = False

        if not isinstance(self.name, str):
            raise TypeError("Name must be str.")

        if not isinstance(self.low_spc, (int, float)):
            raise TypeError("Low_spc must be int or float.")

        if not isinstance(self.upp_spc, (int, float)):
            raise TypeError("Upp_spc must be int or float.")

        if not isinstance(self.res, (int, float)):
            raise TypeError("res must be int or float.")

        if not isinstance(self.spec_units, str):
            raise TypeError("Spec_units must be str.")

        if not isinstance(self.mass_loading, (int, float, type(None))):
            raise TypeError("mass_loading must be int, float or type(None).")

        if not isinstance(self.n, (int, float)):
            raise TypeError("n must be int or float.")

        if not isinstance(self.r, (int, float)):
            raise TypeError("r must be int or float.")

        if not isinstance(self.s, (int, float)):
            raise TypeError("s must be int or float.")

        if not isinstance(self.rho, (int, float, str, type(None))):
            raise TypeError("rho must be int, float or type(None).")

        if not isinstance(self.s_a_den, (int, float, type(None))):
            raise TypeError("s_a_den must be int, float or type(None).")

        if not isinstance(self.v_den, (int, float, type(None))):
            raise TypeError("v_den must be int, float or type(None).")

        if not isinstance(self.dist_type, str):
            raise TypeError("Dist_type must be str.")

        if not isinstance(self.comp, str):
            raise TypeError("Comp must be str.")

        if not isinstance(self.center_alt, (int, float)):
            raise TypeError("Center_alt must be int or float.")

        if not isinstance(self.thick, (int, float)):
            raise TypeError("thick must be int or float.")

        if not isinstance(self.alt_upp, (int, float)):
            raise TypeError("Alt_upp must be int or float.")

        if not isinstance(self.alt_low, (int, float)):
            raise TypeError("Alt_low must be int or float.")

        if not isinstance(self.radii, int):
            raise TypeError("Radii must be int.")

        if not isinstance(self.eta, (int, float)):
            raise TypeError("Eta must be int or float.")

        if not isinstance(self.phase_quad_N, int):
            raise TypeError("Phase_quad_N must be int.")

        if not isinstance(self.phase_quad_type, str):
            raise TypeError("Phase_quad_type must be str.")

        if not isinstance(self.radii_quad_type, str):
            raise TypeError("radii_quad_type must be int or float.")

        if not isinstance(self.leg_coeffs, bool):
            raise TypeError("Leg_coeffs must be bool.")

        if not isinstance(self.leg_coeffs_type, str):
            raise TypeError("Leg_coeffs_type must be str.")

        if not isinstance(self.multiprocess, bool):
            raise TypeError("multiprocess must be bool.")

        passmark = True

        return passmark

    def test_input_values(self):
        """This function tests input values of the MieLayer class.

        The tests written in this function are added as bugs are encountered. The test
        is far from exhaustive and the user is encouraged to sanity check their values
        independently. Creates a boolean value called passmark, which is
        initially False and if the test is passed, is changed to True and returned.

        Returns:
            passmark (bool): If True, input has passed the format test.

        Raises:
            TypeError: Raised when test fails for each parameter separately.

        """

        passmark = False

        if self.low_spc < 0:
            raise ValueError("Attribute low_spc must be non-negative.")

        if self.upp_spc < 0:
            raise ValueError("Attribute upp_spc must be non-negative.")

        if hasattr(self, "mass_loading") and self.mass_loading < 0:
            raise ValueError("Attribute mass_loading must be non-negative.")

        if self.r < 0:
            raise ValueError("Particle mean radius (r) must be non-negative.")

        if self.s < 1:
            raise ValueError("Distribution spread must be <= 1.")

        if self.radii < 1:
            raise ValueError("Number of requested radii must be at least 1.")

        if self.phase_quad_N < 1:
            raise ValueError(
                """Number of scattering angles (phase function quadrature 
            points must be at least 1."""
            )

        if self.thick < 0.002:
            raise ValueError(
                """Thickness must be >= 0.002. This is because due to 
            various roundings in the code, RFM may insert two levels of equal altitude
            and then fail."""
            )

        passmark = True

        return passmark

    def calculate_op(self):
        """Used to calculate the optical properties for the MieLayer class.

        The function calls specialized functions to perform the actual calculations.

        The structure is the following:
            1. Checks the inputs consistency and format and calculates what is missing.
            2. Checks the inputs format (this is a little logically inconsistent
               because some checks must be performed in step one before missing
               properties are calculated).
            3. Calculate size distribution.
            4. Calculate spectral calculation grid.
            5. Calculate the optical properties, which are returned as class attributes.

        Raises:
            RuntimeError: Raised when input fails the format or input tests.

        """

        # calcualtes layer vertical extent related properties
        self.calc_layer_extent()
        self.nsv_or_ml()
        if self.test_complete_input_format():
            pass
        else:
            raise RuntimeError("Input did not pass format test.")

        if self.test_input_values():
            pass
        else:
            raise RuntimeError("Input did not pass value test.")

        self.calc_size_distribution()
        self.calc_grids()
        self.calc_optical_properties()

    def nsv_or_ml(self):
        """Check input for size distribution.

        Checks if either number concentration (n), surface area density (s_a_den) or
        volume density (v_den) are set or is mass loading (mass_loading) is set.
        The logic that s_a_den, v_den and n can be calculated from each other.
        mass_loading can be calculated from n and vice versa.

        Raises:
            RuntimeError: Raised if r and s are not set or if all values are missing.

        """
        if not (hasattr(self, "r") and self.r is not None):
            raise RuntimeError("Mean particle radius (r) must be set.")
        if not (hasattr(self, "s") and self.s is not None):
            raise RuntimeError("Size distribution spread (s) must be set.")

        if (
            (hasattr(self, "n") and self.n is not None)
            or (hasattr(self, "s_a_den") and self.s_a_den is not None)
            or (hasattr(self, "v_den") and self.v_den is not None)
        ):

            self.n_s_v()

            try:
                self.mass_loading = utils.mass_loading_from_number_conc(
                    self.n, self.thick, self.rho, self.s, self.dist_type, self.r
                )
            except AttributeError:
                warnings.warn(
                    """Could not calculate particle loading, because some of 
                thick, rho and r are not set, calculation 
                continues without it."""
                )
            except:
                warnings.warn(
                    """Could not calculate mass_loading, calculation 
                continues without it."""
                )
        elif hasattr(self, "mass_loading") and self.mass_loading is not None:
            self.n = utils.number_conc_from_mass_loading(
                self.mass_loading, self.rho, self.thick, self.s, self.dist_type, self.r
            )
            self.n_s_v()

        else:
            raise RuntimeError(
                """One of number concentration (n), surface area density (s_a_den), 
            volume density (v_den) or mass loading (mass_loading) must be set."""
            )

    def n_s_v(self):
        """Checks inputs for n, s or v.

        Checks if particle either particle number concentration (n), surface area
        density (s_a_den) or v_den are set and calculates the missing of the three.
        Requires meand particle radius (r) and distribution spread (s) to be set.
        If none of the three are set, but particle mass loading (mass_loading) is set
        instead, n will be calculated from particle loading.

        Raises:
            RuntimeError: Raised when calculation fails due to insensible values.
        """

        if hasattr(self, "n") and self.n is not None:
            self.s_a_den = (
                self.n * 4 * np.pi * self.r**2 * np.exp(2 * np.log(self.s) ** 2)
            )
            self.v_den = (
                self.n * 4 * np.pi * self.r**3 * np.exp(9 * np.log(self.s) ** 2 / 2) / 3
            )
        elif hasattr(self, "s_a_den") and self.s_a_den is not None:
            self.n = self.s_a_den / (
                4 * np.pi * self.r**2 * np.exp(2 * np.log(self.s) ** 2)
            )
            self.v_den = (
                self.n * 4 * np.pi * self.r**3 * np.exp(9 * np.log(self.s) ** 2 / 2) / 3
            )
        elif hasattr(self, "v_den") and self.v_den is not None:
            self.n = (
                3
                * self.v_den
                / (4 * np.pi * self.r**3 * np.exp(9 * np.log(self.s) ** 2 / 2))
            )
            self.s_a_den = (
                self.n * 4 * np.pi * self.r**2 * np.exp(2 * np.log(self.s) ** 2)
            )
        else:
            raise RuntimeError(
                """Something is wrong with number concentration, surface 
            area density and volume density. Are input values sensible?"""
            )

        return

    def calc_layer_extent(self):
        """This function attempts to calculate the layer vertical extent.

        If layer center altitude and vertical extent are given, the lower and upper
        boundary altitudes are calculated.
        Else if lower and upper boundaries are given, the center altitude and
        vertical extent are calculated.
        All variables are in units [km].

        Raises:
            RuntimeError: Raised when inputs are insufficient.
        """

        if (hasattr(self, "center_alt") and hasattr(self, "thick")) and (
            self.center_alt is not None and self.thick is not None
        ):
            if (hasattr(self, "alt_upp") and hasattr(self, "alt_low")) and (
                self.alt_upp is not None and self.alt_low is not None
            ):
                assert (self.alt_upp, self.alt_low) == utils.calc_layer_extent(
                    self.center_alt, self.thick
                ), """Both layer (thickness + center altitude) and (upper + lower 
                        boundary altitude) were given, but do not match."""

                return

            else:
                self.alt_upp, self.alt_low = utils.calc_layer_extent(
                    self.center_alt, self.thick
                )
                return

        elif (hasattr(self, "alt_upp") and hasattr(self, "alt_low")) and (
            self.alt_upp is not None and self.alt_low is not None
        ):
            self.center_alt, self.thick = utils.calc_layer_bounds(
                self.alt_upp, self.alt_low
            )
            return

        else:
            raise RuntimeError(
                """To calulate layer altitude characteristics either 
            upper and lower boundary altitude (alt_upp and alt_low) or layer center 
            altitude and vertical extent (center_alt and thick) must be given."""
            )
        return

    def calc_size_distribution(self):
        """Calls functions to calculate the size distribution."""
        self.size_distribution = sz.create_distribution(
            dist_type=self.dist_type,
            n=self.n,
            r=self.r,
            s=self.s,
            surface_area_density=self.s_a_den,
            volume_density=self.v_den,
        )
        return

    def calc_grids(self):
        """Calculates spectral grids for the calculation of optical properties.

        Calculates both wavelengths and wavenumbers.

        """
        self.wvnm, self.wvls = utils.calc_grids(
            lo=self.low_spc, hi=self.upp_spc, res=self.res, units=self.spec_units
        )
        return

    def calc_optical_properties(self):
        """Calculate optical layer properties.

        Calls ``srfm.optical_properties_ewp_hs()`` to calculate extinction coefficient,
        single scatter albedo, phase function and (optional) Legendre polynomial
        coefficients.

        Raises:
            ValueError: Raised when multiprocess isn't specified as attribute.

        """

        if not hasattr(self, "refractive_index"):
            self.refractive_index = None
        if not hasattr(self, "angle"):
            self.angle = None

        if self.multiprocess == False:
            self.op_dict = op.ewp_hs(
                wavelengths=self.wvls,
                composition=self.comp,
                distribution=self.size_distribution,
                refractive_index=self.refractive_index,
                angle=self.angle,
                legendre_coefficients_flag=self.leg_coeffs,
                legendre_coefficients_type=self.leg_coeffs_type,
                radii=self.radii,
                radii_quad_type=self.radii_quad_type,
                eta=self.eta,
                phase_quad_N=self.phase_quad_N,
                phase_quad_type=self.phase_quad_type,
                return_dict=None,
                multiprocess=False,
            )

        elif self.multiprocess == True:
            self.op_dictproxy = Manager().dict()
            self.op_process = Process(
                target=op.ewp_hs,
                kwargs={
                    "wavelengths": self.wvls,
                    "composition": self.comp,
                    "distribution": self.size_distribution,
                    "refractive_index": self.refractive_index,
                    "angle": self.angle,
                    "legendre_coefficients_flag": self.leg_coeffs,
                    "legendre_coefficients_type": self.leg_coeffs_type,
                    "radii": self.radii,
                    "eta": self.eta,
                    "phase_quad_N": self.phase_quad_N,
                    "phase_quad_type": self.phase_quad_type,
                    "radii_quad_type": self.radii_quad_type,
                    "return_dict": self.op_dictproxy,
                    "multiprocess": True,
                },
            )
            self.op_process.start()

        else:
            raise ValueError("multiprocess must be True or False.")
        return

    def add_op_calc_output(self):
        """Assigns results from optical properties calculation as class attributes.

        Takes results from optical properties calcualtion stored in a dictionary and
        assignes them as class attributes. This is a separate function so that when
        multiprocess is True, this step (adding results to class can be done after
        threads join in the main script in one command.
        This function deletes the original dictionary.

        """

        if hasattr(self, "op_process"):
            self.op_process.join()
            self.op_dict = {}
            self.op_dict.update(self.op_dictproxy)
            delattr(self, "op_dictproxy")

        if hasattr(self, "op_dict"):
            self.beta_ext = self.op_dict["beta_ext"]
            self.ssalb = self.op_dict["ssalb"]
            self.phase_function = self.op_dict["phase_function"]
            if "legendre_coefficient" in self.op_dict.keys():
                self.legendre_coefficient = self.op_dict["legendre_coefficient"]
            delattr(self, "op_dict")
        else:
            raise RuntimeError(
                """Optical properties (op_dict) not found. Were they 
            calcualted?"""
            )

    @utils.show_runtime
    def regrid(self, wvls, track_diff=False, diff_type="pct"):
        """Linearly interpolates calculated values from ewp_hs to a new grid.

        Instance must have  "beta_ext", "ssalb", "phase_function", and
        "legendre_coefficient".

        Args:

            wvls (aray-like): new grid, units [\ :math:`\\mu`\ m]
            track_diff (bool): If True, calculates differences arising from
                interpolation.
        Returns:
            old attributes now as "_old"
            new attributes "beta_ext", "ssalb", "phase_function" and
            optional ("legendre_coefficient")
            difference between old and new attributes: "_diff" version of the new
              attributes

        Raises:
            ValueError: Raised when any wavelength array is not monotonic.

        """
        # check if legendre coefficients are present:
        lc_flag = hasattr(self, "legendre_coefficient")

        # preserve original results
        self.wvls_old = self.wvls
        self.beta_ext_old = self.beta_ext
        self.ssalb_old = self.ssalb
        self.phase_function_old = self.phase_function
        if lc_flag == True:
            self.legendre_coefficient_old = self.legendre_coefficient

        # check if new input wavelengths are monotonic and increasing or decreasing
        wvls_mono = utils.monotonic(wvls)
        if wvls_mono == 0:
            raise ValueError("Wvls is not monotonic.")
        elif wvls_mono == 2:
            wvls = np.flip(wvls)

        # check if existing wavelengths are monotoning and increasing or decreasing
        cls_mono = utils.monotonic(self.wvls)
        if cls_mono == 0:
            raise ValueError("class wavelengths are not monotonic.???")
        elif cls_mono == 2:
            self.wvls = np.flip(self.wvls)
            self.beta_ext = np.flip(self.beta_ext)
            self.ssalb = np.flip(self.ssalb)
            self.phase_function = np.flipud(self.phase_function)
            if lc_flag == True:
                self.legendre_coefficient = np.flipud(self.legendre_coefficient)

        # interpolate
        self.beta_ext = np.interp(wvls, self.wvls, self.beta_ext)
        self.ssalb = np.interp(wvls, self.wvls, self.ssalb)
        self.phase_function = np.array(
            [
                np.interp(wvls, self.wvls, self.phase_function[:, i])
                for i in range(self.phase_function.shape[1])
            ]
        ).T
        if lc_flag == True:
            self.legendre_coefficient = np.array(
                [
                    np.interp(wvls, self.wvls, self.legendre_coefficient[:, i])
                    for i in range(self.legendre_coefficient.shape[1])
                ]
            ).T
        self.wvls = wvls

        if wvls_mono == 2:
            self.wvls = np.flip(self.wvls)
            self.beta_ext = np.flip(self.beta_ext)
            self.ssalb = np.flip(self.ssalb)
            self.phase_function = np.flipud(self.phase_function)
            if lc_flag == True:
                self.legendre_coefficient = np.flipud(self.legendre_coefficient)

        if track_diff == True:
            self.track_regrid_diff(diff_type=diff_type)

        return

    def track_regrid_diff(self, diff_type="pct"):
        """Calculates the difference before and after interpolation.

        Called by optical_properties.regrid().

        Args:
            diff_type (str): "pct" (differences in percent from the old value) or
                "abs" (absolute difference)

        Raises:
            ValueError: Raised when diff_type is invalid.

        """
        # check if legendre coefficients present
        lc_flag = hasattr(self, "legendre_coefficient")

        diff_dict = {}
        diff_dict["wavelengths"] = self.wvls

        temp_dict = {}
        temp_dict["wavelengths"] = self.wvls
        temp_dict["beta_ext"] = []
        temp_dict["ssalb"] = []
        temp_dict["phase_function"] = []
        if lc_flag == True:
            temp_dict["legendre_coefficient"] = []

        for i in self.wvls:
            x = np.argmin(np.abs(self.wvls_old - np.full(self.wvls_old.size, i)))
            temp_dict["beta_ext"].append(self.beta_ext_old[x])
            temp_dict["ssalb"].append(self.ssalb_old[x])
            temp_dict["phase_function"].append(self.phase_function_old[x])
            if lc_flag == True:
                temp_dict["legendre_coefficient"].append(
                    self.legendre_coefficient_old[x]
                )

        for key in temp_dict.keys():
            if key == "wavelengths":
                pass
            else:
                temp_dict[key] = np.asarray(temp_dict[key])

        self.beta_ext_diff = np.abs(temp_dict["beta_ext"] - self.beta_ext)
        self.ssalb_diff = np.abs(temp_dict["ssalb"] - self.ssalb)
        self.phase_function_diff = np.abs(
            temp_dict["phase_function"] - self.phase_function
        )
        if lc_flag == True:
            self.legendre_coefficient_diff = np.abs(
                temp_dict["legendre_coefficient"] - self.legendre_coefficient
            )

        if diff_type == "abs":
            return

        elif diff_type == "pct":
            self.beta_ext_diff = (
                np.nan_to_num(np.divide(self.beta_ext_diff, temp_dict["beta_ext"]))
                * 100
            )
            self.ssalb_diff = (
                np.nan_to_num(np.divide(self.ssalb_diff, temp_dict["ssalb"])) * 100
            )
            self.phase_function_diff = (
                np.nan_to_num(
                    np.divide(self.phase_function_diff, temp_dict["phase_function"])
                )
                * 100
            )
            if lc_flag == True:
                self.legendre_coefficient_diff = (
                    np.nan_to_num(
                        np.divide(
                            self.legendre_coefficient_diff,
                            temp_dict["legendre_coefficient"],
                        )
                    )
                    * 100
                )
            return
        else:
            raise ValueError("diff_type must be pct or abs.")

        return

    def plot_diff(self, **kwargs):
        """Plots difference in calcualted optical properties.

        Instance must have differences in optical properties calculated and stored
        as _diff attributes. For further information see
        ``srfm.layer.track_regrid_diff()``.

        Returns:
            fig (obj): fig object
        """
        fig = plt.figure(figsize=(11.7, 8.3))

        plt.subplot(3, 1, 1)
        plt.plot(self.wvls, self.beta_ext_diff, label="beta_ext")
        plt.xlabel("wavelength (um)", fontsize="12")
        plt.ylabel(r"$\Delta \beta ^{ext} (\lambda)$ (%)", fontsize="12")
        # plt.yscale("log")
        plt.title("Extinction coefficient", fontsize="15")

        plt.subplot(3, 1, 2)
        plt.plot(self.wvls, self.ssalb_diff, label="ssalb")
        plt.xlabel("wavelength (um)", fontsize="12")
        plt.ylabel(r"$\Delta \omega (\lambda)$ (%)", fontsize="12")
        # plt.yscale("log")
        plt.title("Single scattering albedo", fontsize="15")

        plt.subplot(3, 1, 3)
        [
            plt.plot(
                self.wvls,
                self.phase_function_diff[:, i],
                label=f"phase function, {i} deg",
            )
            for i in range(self.phase_function_diff.shape[1])
        ]
        plt.xlabel("wavelength (um)", fontsize="12")
        plt.ylabel(r"$\Delta p(\mu, \lambda)$ (%)", fontsize="12")
        # plt.yscale("log")
        # plt.legend()
        plt.title("Phase function - each line represents one angle", fontsize="15")

        #        plt.subplot(4, 1, 4)
        #        [
        #            plt.plot(
        #                diff_dict["wavelengths"],
        #                diff_dict["legendre_coefficient"][:, i],
        #                label=f"legendre coefficient no. {i}",
        #            )
        #            for i in range(diff_dict["legendre_coefficient"].shape[1])
        #        ]
        #        plt.xlabel("wavelength (um)")
        #        plt.ylabel("legendre coefficients")
        #        plt.yscale("log")
        #        # plt.legend()

        plt.suptitle(
            "Difference in optical properties arising from interpolation", fontsize="18"
        )
        plt.tight_layout()
        return fig

    def calc_tau(self):
        """Calculates optical depth from beta_ext and layer vertical extent.

        Both beta_ext and layer boundary altitudes must be specified in the object
        instance.
        """
        self.tau = self.beta_ext * 1e3 * (self.alt_upp - self.alt_low)
        return


class GreyBodyCloud(Layer):
    """Class that represents a grey body cloud.

    Is subclass of Layer.
    The rationale for this class is to emulate the simplest possible case where a cloud
    is represented as a greybody emitter. No scattering is involved.

    Unless changed, emissivity is 1 by default.

    Todo:
        The emissivity is a stub, the attribute isn't implemented anywhere. For the time
        being then, the GreBodyCloud is effectively a black body emitter.

    """

    def __init__(self, name=None, emis=1, **parameters):
        super().__init__(name)
        for key, val in parameters.items():
            self.key = val
        self.emis = emis

    def set_name(self, name):
        """Assign a name to the layer.

        Args:
            name (str): Layer name.

        """
        self.name = name

    def set_spc_lim(self, lo, hi):
        """Set spectral grid lower and upper limits.

        Lo is the lower value and hi is the upper value, whatever the units. Units are
        then specified in set_spec_units, so this remains consistent with both
        wavenumbers and wavelengths.

        Args:
            lo (int, float): Lower wavenumber or wavelength.
            hi (int, float): Upper wavenumber or wavelength.

        """
        self.low_spc = lo
        self.upp_spc = hi

    def set_spec_units(self, units):
        """Set spectral calculation grid units.

        Args:
            units (str): Units. should be one of [\ :math:`\\mu`\ m, cm\ :sup:`-1`, nm].
                This is not enforced but currently the rest of this package is unable to
                handle any other options.

        """
        self.spec_units = units

    def set_res(self, res):
        """Set spectral calculation grid resolution.

        Args:
            res (int, float): Resolution of the spectral grid. Should be same units as
                values in set_spc_lim. This is not enforced, but is the only logically
                consistent option.

        """
        self.res = res

    def set_center_alt(self, center_alt):
        """Set average (center) particle layer altitude.

        Args:
            center_alt (int, float): Center particle layer altitude, units [km].

        """
        self.center_alt = center_alt

    def set_thick(self, thick):
        """Set layer thickness (vertical extent).

        Args:
            thick (int, float): Layer thickness (vertical extent), units [km].

        """
        self.thick = thick

    def set_alt_lim(self, alt_low, alt_upp):
        """Set lower and upper layer boundary (altitudes).

        Args:
            alt_low (int, float): Layer lower boundary altitude, units [km].
            alt_upp (int, float): Layer upper boundary altitude, units [km].

        """
        self.alt_upp = alt_upp
        self.alt_low = alt_low

    def set_emis(self, e):
        """Set emissivity of the clouds.

        Args:
            e (int, float): Grey body emissivity.

        """
        self.emis = e

    def set_tau(self, tau):
        """Set input layer optical depth.

        Args:
            tau (int, float): Layer optical depth. The optical depth is uniform across
                all wavelengths (since this is class represents a grey body cloud).

        """

        self.inp_tau = tau

    def set_input_from_dict(self, inp_dict):
        """Set all necessary input_parameters from an input dictionary.

        Args:
            inp_dict (dict): Input dictonary with layer properties. Must contain
                all necessary keys (cannot be incomplete). Required keys are:

                    - name
                    - low_spc
                    - upp_spc
                    - spec_units
                    - res
                    - center_alt
                    - thick
                    - alt_upp
                    - alt_low
                    - emis
                    - inp_tau

                For explanation of each of those parameters please refer to the
                respective function which set them explicitly (set_{parameter name}).

        """
        self.name = inp_dict["name"]
        self.low_spc = inp_dict["low_spc"]
        self.upp_spc = inp_dict["upp_spc"]
        self.spec_units = inp_dict["spec_units"]
        self.res = inp_dict["res"]
        self.center_alt = inp_dict["center_alt"]
        self.thick = inp_dict["thick"]
        self.alt_upp = inp_dict["alt_upp"]
        self.alt_low = inp_dict["alt_low"]
        self.emis = inp_dict["emis"]
        self.inp_tau = inp_dict["inp_tau"]

    def test_complete_input_format(self):
        """Checks the input for correct format before running any calculation.

        Tests all input variables. Creates a boolean value called passmark, which is
        initially False and if the test is passed, is changed to True and returned.

        Returns:
            passmark (bool): If True, input has passed the format test.

        Raises:
            TypeError: Raised when test fails for each parameter separately.

        """
        passmark = False

        if not isinstance(self.name, str):
            raise TypeError("Name must be str.")

        if not isinstance(self.low_spc, (int, float)):
            raise TypeError("Low_spc must be int or float.")

        if not isinstance(self.upp_spc, (int, float)):
            raise TypeError("Upp_spc must be int or float.")

        if not isinstance(self.res, (int, float)):
            raise TypeError("res must be int or float.")

        if not isinstance(self.spec_units, str):
            raise TypeError("Spec_units must be str.")

        if not isinstance(self.center_alt, (int, float)):
            raise TypeError("Center_alt must be int or float.")

        if not isinstance(self.thick, (int, float)):
            raise TypeError("thick must be int or float.")

        if not isinstance(self.alt_upp, (int, float)):
            raise TypeError("Alt_upp must be int or float.")

        if not isinstance(self.alt_low, (int, float)):
            raise TypeError("Alt_low must be int or float.")

        if not isinstance(self.emis, (int, float)):
            raise TypeError("Emis must be int or float.")

        if not isinstance(self.inp_tau, (int, float)):
            raise TypeError("Inp_tau must be int or float.")

        passmark = True

        return passmark

    def test_input_values(self):
        """This function tests input values of the GreyBodyCloud class.

        The tests written in this function are added as bugs are encountered. The test
        is far from exhaustive and the user is encouraged to sanity check their values
        independently. Creates a boolean value called passmark, which is
        initially False and if the test is passed, is changed to True and returned.

        Returns:
            passmark (bool): If True, input has passed the format test.

        Raises:
            TypeError: Raised when test fails for each parameter separately.

        """

        passmark = False

        if self.low_spc < 0:
            raise ValueError("Attribute low_spc must be non-negative.")

        if self.upp_spc < 0:
            raise ValueError("Attribute upp_spc must be non-negative.")

        if self.emis < 0 or self.emis > 1:
            raise ValueError("Emissivity must be between 0 and 1.")

        if self.thick < 0.002:
            raise ValueError(
                """Thickness must be >= 0.002. This is because due to 
            various roundings in the code, RFM may insert two levels of equal altitude
            and then fail."""
            )

        passmark = True

        return passmark

    def calculate_op(self):
        """Used to calculate the optical properties for the GreyBodyCloud class.

        The function calls specialized functions to perform the actual calculations.

        The structure is the following:
            1. Checks the inputs consistency and format and calculates what is missing.
            2. If emissivity not set, defaults to 1.
            3. Checks the inputs format (this is a little logically inconsistent
               because some checks must be performed in step one before missing
               properties are calculated).
            4. Calculate spectral grid.
            5. Calculate the optical properties, which are returned as class attributes.

        Raises:
            RuntimeError: Raised when input fails the format or input tests.

        """

        # calcualtes layer vertical extent related properties
        self.calc_layer_extent()

        if self.test_complete_input_format():
            pass
        else:
            raise RuntimeError("Input did not pass format test.")

        if self.test_input_values():
            pass
        else:
            raise RuntimeError("Input did not pass value test.")

        self.calc_grids()
        self.calc_optical_properties()

    def calc_layer_extent(self):
        """This function attempts to calculate the layer vertical extent.

        If layer center altitude and vertical extent are given, the lower and upper
        boundary altitudes are calculated.
        Else if lower and upper boundaries are given, the center altitude and
        vertical extent are calculated.
        All variables are in units [km].

        Raises:
            RuntimeError: Raised when inputs are insufficient.
        """

        if (hasattr(self, "center_alt") and hasattr(self, "thick")) and (
            self.center_alt is not None and self.thick is not None
        ):
            if (hasattr(self, "alt_upp") and hasattr(self, "alt_low")) and (
                self.alt_upp is not None and self.alt_low is not None
            ):
                assert (self.alt_upp, self.alt_low) == utils.calc_layer_extent(
                    self.center_alt, self.thick
                ), """Both layer (thickness + center altitude) and (upper + lower 
                        boundary altitude) were given, but do not match."""

                return

            else:
                self.alt_upp, self.alt_low = utils.calc_layer_extent(
                    self.center_alt, self.thick
                )
                return

        elif (hasattr(self, "alt_upp") and hasattr(self, "alt_low")) and (
            self.alt_upp is not None and self.alt_low is not None
        ):
            self.center_alt, self.thick = utils.calc_layer_bounds(
                self.alt_upp, self.alt_low
            )
            return

        else:
            raise RuntimeError(
                """To calulate layer altitude characteristics either 
            upper and lower boundary altitude (alt_upp and alt_low) or layer center 
            altitude and vertical extent (center_alt and thick) must be given."""
            )
        return

    def calc_grids(self):
        """Calculates spectral grids for the calculation of optical properties.

        Calculates both wavelengths and wavenumbers.

        """
        self.wvnm, self.wvls = utils.calc_grids(
            lo=self.low_spc, hi=self.upp_spc, res=self.res, units=self.spec_units
        )
        return

    def calc_optical_properties(self):
        """Calculated optical properties.

        In this case all it does is take inp_tau and stretches it over the spectral
        grid. The fact that this is a separate function is because the implementation
        mirrors implementation scheme in other *srfm.layer.Layer()* objects.

        """

        self.tau = [self.inp_tau] * len(self.wvnm)
        return

    @utils.show_runtime
    def regrid(self, wvls, track_diff=False, diff_type="pct"):
        """Linearly interpolates calculated values from ewp_hs to a new grid.

        Instance must have the "tau" attribute.

        Args:

            wvls (aray-like): new grid, units [\ :math:`\\mu`\ m]
            track_diff (bool): If True, calculates differences arising from
                interpolation.
        Returns:
            old attribute tau now as "tau_old"
            new attribute "tau"
            difference between old and new attribute "tau_diff"

        Raises:
            ValueError: Raised when any wavelength array is not monotonic.

        """

        # preserve original results
        self.wvls_old = self.wvls
        self.tau_old = self.tau

        # check if new input wavelengths are monotonic and increasing or decreasing
        wvls_mono = utils.monotonic(wvls)
        if wvls_mono == 0:
            raise ValueError("Wvls is not monotonic.")
        elif wvls_mono == 2:
            wvls = np.flip(wvls)

        # check if existing wavelengths are monotoning and increasing or decreasing
        cls_mono = utils.monotonic(self.wvls)
        if cls_mono == 0:
            raise ValueError("Instance wavelengths are not monotonic.???")
        elif cls_mono == 2:
            self.wvls = np.flip(self.wvls)
            self.tau = np.flip(self.tau)

        # interpolate
        self.tau = np.interp(wvls, self.wvls, self.tau)
        self.wvls = wvls

        if wvls_mono == 2:
            self.wvls = np.flip(self.wvls)
            self.tau = np.flip(self.tau)

        if track_diff == True:
            self.track_regrid_diff(diff_type=diff_type)

        return

    def track_regrid_diff(self, diff_type="pct"):
        """Calculates the difference before and after interpolation.

        Called by optical_properties.regrid().

        Args:
            diff_type (str): "pct" (differences in percent from the old value) or
                "abs" (absolute difference)

        Raises:
            ValueError: Raised when diff_type is invalid.

        """

        diff_dict = {}
        diff_dict["wavelengths"] = self.wvls

        temp_dict = {}
        temp_dict["wavelengths"] = self.wvls
        temp_dict["tau"] = []

        for i in self.wvls:
            x = np.argmin(np.abs(self.wvls_old - np.full(self.wvls_old.size, i)))
            temp_dict["tau"].append(self.tau_old[x])

        for key in temp_dict.keys():
            if key == "wavelengths":
                pass
            else:
                temp_dict[key] = np.asarray(temp_dict[key])

        self.tau_diff = np.abs(temp_dict["tau"] - self.tau)

        if diff_type == "abs":
            return

        elif diff_type == "pct":
            self.tau_diff = (
                np.nan_to_num(np.divide(self.tau_diff, temp_dict["tau"])) * 100
            )
            return
        else:
            raise ValueError("diff_type must be pct or abs.")

        return

    def plot_diff(self, **kwargs):
        """Plots difference in calcualted optical properties.

        Instance must have differences in optical properties calculated and stored
        as _diff attributes. For further information see
        ``srfm.layer.track_regrid_diff()``.

        Returns:
            fig (obj): fig object
        """
        fig = plt.figure(figsize=(11.7, 8.3))

        plt.subplot(1, 1, 1)
        plt.plot(self.wvls, self.tau_diff, label="tau")
        plt.xlabel("wavelength (um)", fontsize="12")
        plt.ylabel(r"$\Delta \Tau (\lambda)$ (%)", fontsize="12")
        # plt.yscale("log")
        plt.title("Optical depth", fontsize="15")

        plt.suptitle(
            "Difference in optical properties arising from interpolation", fontsize="18"
        )
        plt.tight_layout()
        return fig
