"""
Name: layer
Parent package: srfm
Author: Antonin Knizek
Contributors: 
Date: 26 Mar 2025
Purpose: Defines the Layer class and all its subclasses, used to contain scattering
         information for a single layer.
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
    """Class that contains mie scattering parameters, calculations and outputs."""

    def __init__(self, name=None, **parameters):
        super().__init__(name)
        for key, val in parameters.items():
            self.key = val

    def set_spc_lim(self, lo, hi):
        """Set spectral calculation grid lower and upper limits."""
        self.low_spc = lo
        self.upp_spc = hi

    def set_spec_units(self, units):
        """Set spectral calculation grid units."""
        self.spec_units = units

    def set_res(self, res):
        """Set spectral calculation grid resolution."""
        self.res = res

    def set_mass_loading(self, load):
        """Set particle mass loading, units [g m-2]"""
        self.mass_loading = load

    def set_n(self, n):
        """Set particle number concentration, units [cm-1]."""
        self.n = n

    def set_r(self, r):
        """Set mean particle radius, units [um]."""
        self.r = r

    def set_s(self, s):
        """Set particle size distribution spread."""
        self.s = s

    def set_rho(self, rho):
        """Set particle mean density, units [kg m-3] or one of strings accepted by
        utilities.number_conc_from_mass_loading."""
        self.rho = rho

    def set_s_a_den(self, s_a_den):
        """Set surface area density of the particles, units [um2 cm-3]."""
        self.s_a_den = s_a_den

    def set_v_den(self, v_den):
        """Set volume density of particles, units [um3 cm-3]."""
        self.v_den = v_den

    def set_dist_type(self, dist_type):
        """Set particle size distribution type (lognormal, normal)."""  #
        self.dist_type = dist_type

    def set_comp(self, comp):
        """Set particle composition type, one of accepted values by
        ARIA_module.get_ri_filepathname."""
        self.comp = comp

    def set_center_alt(self, center_alt):
        """Set average (center) particle layer altitude, units [km]."""
        self.center_alt = center_alt

    def set_thick(self, thick):
        """Set layer thickness (vertical extent), units [km]."""
        self.thick = thick

    def set_alt_lim(self, alt_low, alt_upp):
        """Set lower and upper layer boundary (altitudes), units [km]."""
        self.alt_upp = alt_upp
        self.alt_low = alt_low

    def set_radii(self, radii=200):
        """Set number of radii for size disitribution calculation, default 200."""
        self.radii = radii

    def set_eta(self, eta=1e-6):
        """Set size distribution cut-off (eta value), default 1e-6."""
        self.eta = eta

    def set_phase_quad_N(self, phase_quad_N=181):
        """Set number of points for the scattering quadrature (number of scattering
        angles), default 181."""
        self.phase_quad_N = phase_quad_N

    def set_phase_quad_type(self, phase_quad_type="L"):
        """Set quadrature type that will be used to calculate the phase fucntion,
        default L (Lobatto), for other options see the quadrature module."""
        self.phase_quad_type = phase_quad_type

    def set_radii_quad_type(self, radii_quad_type="T"):
        """Set quadrature type for the calculation of the size distribution,
        default T (trapezium), for other options see the quadrature module."""
        self.radii_quad_type = radii_quad_type

    def set_leg_coeffs(self, leg_coeffs=True):
        """Toggle expansion of the phase fucntion into Legendre polynomials, type bool,
        default True."""
        self.leg_coeffs = leg_coeffs

    def set_leg_coeffs_type(self, leg_coeffs_type="normalised"):
        """Set type of requested Legendre polynomial expansion coefficients, options are
        regular or normalised (default).
        """
        self.leg_coeff_type = leg_coeffs_type
    
    def set_aria(self,aria):
        """Sets path to ARIA folder."""
        self.aria = aria

    def set_multiproccess(self, multiprocess):
        """Toggle multiprocessing, True or False."""
        self.multiprocess = multiprocess

    def set_input_from_dict(self, inp_dict):
        """Set all necessary input_parameters from an input dictionary."""
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
        self.aria = inp_dict["aria"]
        self.multiprocess = inp_dict["multiprocess"]

    def test_complete_input_format(self):
        """Checks the input for correct format before running any calculation.
        Tests all input variables."""
        passmark = False

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

    def calculate_op(self):
        """This function calculates the optical property for the MieLayer class.
        The structure is the followng:
            1. checks the inputs consistency and format and calculates what is missing
            2. checks the inputs format (this is a little logically inconsistent,
               because some checks must be performed in step one before missing
               properties are calculated).
            3. calculate size distribution
            4. calculate spectral calculation grid
            5. calculate the optical properties, which are returned as class attributes
        """

        # calcualtes layer vertical extent related properties
        self.calc_layer_extent()
        self.nsv_or_ml()
        if self.test_complete_input_format():
            pass
        else:
            raise RuntimeError("Input did not pass format test.")

        self.calc_size_distribution()
        self.calc_grids()
        self.calc_optical_properties()

    def nsv_or_ml(self):
        """Checks if either number concentration (n), surface area density (s_a_den) or
        volume density (v_den) are set or is mass loading (mass_loading) is set.
        The logic is:
        
        s_a_den
           |    \
           |     \
           |      n <---> mass_loading
           |     /
           |    /
         v_den
        """
        if not (hasattr(self, "r") and self.r != None):
            raise RuntimeError("Mean particle radius (r) must be set.")
        if not (hasattr(self, "s") and self.s != None):
            raise RuntimeError("Size distribution spread (s) must be set.")

        if (
            (hasattr(self, "n") and self.n != None)
            or (hasattr(self, "s_a_den") and self.s_a_den != None)
            or (hasattr(self, "v_den") and self.v_den != None)
        ):

            self.n_s_v()

            try:
                self.mass_loading = utils.mass_loading_from_number_conc(
                    self.n, self.thick, self.rho, self.r
                )
            except AttributeError:
                warnings.warn(
                    """Could not calculate particle loading, because some of 
                thick, rho and r are not set, calculation 
                continues without it."""
                )
            except:
                warning.warn(
                    """Could not calculate mass_loading, calculation 
                continues without it."""
                )
        elif hasattr(self, "mass_loading") and self.mass_loading != None:
            self.n = utils.number_conc_from_mass_loading(
                self.mass_loading, self.rho, self.thick, self.r
            )
            self.n_s_v()

        else:
            raise RuntimeError(
                """One of number concentration (n), surface area density (s_a_den), 
            volume density (v_den) or mass loading (mass_loading) must be set."""
            )

    def n_s_v(self):
        """Checks if particle either particle number concentration (n), surface area
        density (s_a_den) or v_den are set and calculates the missing of the three.
        Requires meand particle radius (r) and distribution spread (s) to be set.
        If none of the three are set, but particle mass loading (mass_loading) is set
        instead, n will be calculated from particle loading.
        """

        if hasattr(self, "n") and self.n != None:
            self.s_a_den = (
                self.n * 4 * np.pi * self.r**2 * np.exp(2 * np.log(self.s) ** 2)
            )
            self.v_den = (
                self.n * 4 * np.pi * self.r**3 * np.exp(9 * np.log(self.s) ** 2 / 2) / 3
            )
        elif hasattr(self, "s_a_den") and self.s_a_den != None:
            self.n = self.s_a_den / (
                4 * np.pi * self.r**2 * np.exp(2 * np.log(self.s) ** 2)
            )
            self.v_den = (
                self.n * 4 * np.pi * self.r**3 * np.exp(9 * np.log(self.s) ** 2 / 2) / 3
            )
        elif hasattr(self, "v_den") and self.v_den != None:
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
        """

        if (hasattr(self, "center_alt") and hasattr(self, "thick")) and (
            self.center_alt != None and self.thick != None
        ):
            if (hasattr(self, "alt_upp") and hasattr(self, "alt_low")) and (
                self.alt_upp != None and self.alt_low != None
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
            self.alt_upp != None and self.alt_low != None
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
        """Calculates spectral grids in wavenumbers and wavelengths for the calculation
        of optical properties."""
        self.wvnm, self.wvls = utils.calc_grids(
            lo=self.low_spc, hi=self.upp_spc, res=self.res, units=self.spec_units
        )
        return

    def calc_optical_properties(self):
        """Calls function ewp_hs from optical_properties to calculate optical
        layer properties - extinction coefficient, single scatter albedo, phase function
        and (optional) Legendre polynomial coefficients.
        multiprocess must be specified.
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
                aria=self.aria,
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
                    "aria": self.aria,
                    "return_dict": self.op_dictproxy,
                    "multiprocess": True,
                },
            )
            self.op_process.start()

        else:
            raise ValueError("multiprocess must be True or False.")
        return

    def add_op_calc_output(self):
        """Takes results from optical properties calcualtion stored in a dictionary and
        assignes them as class attributes. This is a separate function so that when
        multiprocess, this step can be done after threads join in the main script in one
        command.
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
        inputs:
            self, must have "beta_ext", "ssalb", "phase_function",
            "legendre_coefficient"
            wvls - new grid, [um]
            track_diff - bool, if True, calculates differences arising
                         from interpolation
        outputs:
            old attributes now as "_old"
            new attributes "beta_ext", "ssalb", "phase_function" and
            optional ("legendre_coefficient")
            difference between old and new attributes: "_diff" version of the new
            attributes
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
                    np.interp(
                        wvls, self.wvls, self.legendre_coefficient[:, i]
                    )
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
        """Calculates the difference in op_dict (dictionary containing ewp_hs outputs)
        before and after interpolation, called by optical_properties.regrid().
        inputs:
            self
            diff_type - "pct" (differences in percent from the old value)
                      - "abs" - absolute difference
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
        self.phase_function_diff = np.abs(temp_dict["phase_function"] - self.phase_function)
        if lc_flag == True:
            self.legendre_coefficient_diff = np.abs(temp_dict["legendre_coefficient"] - self.legendre_coefficient)
        
        if diff_type == "abs":
            return
        
        elif diff_type == "pct":
            self.beta_ext_diff = np.nan_to_num(np.divide(self.beta_ext_diff,temp_dict["beta_ext"])) * 100
            self.ssalb_diff = np.nan_to_num(np.divide(self.ssalb_diff,temp_dict["ssalb"])) * 100
            self.phase_function_diff = np.nan_to_num(np.divide(self.phase_function_diff,temp_dict["phase_function"])) * 100
            if lc_flag == True:
                self.legendre_coefficient_diff = np.nan_to_num(np.divide(self.legendre_coefficient_diff,temp_dict["legendre_coefficient"])) * 100
            return        
        else:
            raise ValueError("diff_type must be pct or abs.")

        return

    def plot_diff(self, **kwargs):
        """Plots difference in calcualted optical properties.
        inputs:
            self - dictionary with differences in optical properties
            comes e.g. from self.track_regrid_diff
        outputs:
            fig object
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
