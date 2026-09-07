"""Defines the base model class for the forward model and specific model subclasses.

- Name: forward_model
- Parent package: srfm
- Author: Antonin Knizek
- Contributors:
- Date: 18 February 2025
"""

import numpy as np
from . import disort_functions as disf
from . import rfm_functions as rf
import os
from . import utilities as utils
from . import units
import pandas as pd
from scipy.signal import convolve
import scipy.constants
from scipy.interpolate import make_interp_spline
from dataclasses import dataclass

try:
    from .DISORT import disort_module_s as dms
except ImportError:
    print("Could not import disort as a python module.")
    print('Is disort compiled?"')
    print("Is the compiled module called disort_module?")
    print("Is disort compiled in its folder?")
    print("hint: run prepare_disort.sh in the DISORT folder.")

try:
    from .DISORT_dbl import disort_module_d as dmd
except ImportError:
    print("Could not import disort as a python module.")
    print('Is disort compiled?"')
    print("Is the compiled module called disort_module?")
    print("Is disort compiled in its folder?")
    print("hint: run prepare_disort.sh in the DISORT folder.")


class Fwd_model:
    """Base class for containing forward models.

    Specific forward models are defined as subclasses of this superclass.

    """

    def __init__(self, name=None, **parameters):
        self.name = name
        self.parameters = {}


@dataclass(frozen=True)
class DisortResult:
    """Store the normalized outputs from one DISORT calculation.

    The object deliberately contains only one spectral point. Long-running SRFM
    calculations can therefore copy the values into their destination arrays without
    retaining a dictionary entry for every wavenumber.

    Args:
        wavenumber: Central wavenumber of the calculation in cm-1.
        wavelength: Central wavelength of the calculation in micrometres.
        rfldir: Direct-beam flux at the requested output levels.
        rfldn: Diffuse downward flux at the requested output levels.
        flup: Diffuse upward flux at the requested output levels.
        dfdt: Flux-divergence derivative at the requested output levels.
        uavg: Mean intensity at the requested output levels.
        uu: User-angle radiance.
        albmed: Medium albedo.
        trnmed: Medium transmissivity.
    """

    wavenumber: float
    wavelength: float
    rfldir: np.ndarray
    rfldn: np.ndarray
    flup: np.ndarray
    dfdt: np.ndarray
    uavg: np.ndarray
    uu: np.ndarray
    albmed: np.ndarray
    trnmed: np.ndarray | float

    def as_dict(self):
        """Return the legacy fixed-key dictionary representation.

        Returns:
            dict: DISORT outputs using the historical public key names.
        """
        return {
            "wavenumber (cm-1)": self.wavenumber,
            "wavelength (um)": self.wavelength,
            "rfldir": self.rfldir,
            "rfldn": self.rfldn,
            "flup": self.flup,
            "dfdt": self.dfdt,
            "uavg": self.uavg,
            "uu": self.uu,
            "albmed": self.albmed,
            "trnmed": self.trnmed,
        }


class RFM(Fwd_model):
    """Class that contains the RFM forward model.

    Is subclass of Fwd_model.
    """

    def __init__(
        self,
        name="RFM",
        rfm_fldr=None,
        status="RFM model object created.",
        **parameters,
    ):
        super().__init__(name)
        self.rfm_fldr = rfm_fldr
        self.status = status
        for key, val in parameters.items():
            setattr(self, key, val)

    @utils.show_runtime
    def run_rfm(self, fldr, wipe=True):
        """Runs RFM from python.

        Args:
            fldr (str): Path to the RFM folder (relative or absolute path).
            wipe (bool): If True, removes rfm results from the previous run.
                Default is True.

        Returns:
            Runs RFM, write its outputs to the fldr folder, updates object status.

        """

        cwd = os.getcwd()
        try:
            os.chdir(f"{fldr}")
        except OSError:
            print("Invalid directory.")
            return

        if wipe == True:
            try:
                _ = [os.remove(str(i)) for i in os.listdir() if i.endswith(".asc")]
                _ = [os.remove(str(i)) for i in os.listdir() if i.endswith(".log")]
            except FileNotFoundError:
                pass
            except PermissionError:
                print("Do not have the permission to remove the file.")

        try:
            os.system("./source/rfm")
        except:
            print("Could not run rfm, check manually for errors.")

        os.chdir(cwd)
        self.status = "RFM run completed, result not yet loaded."
        return

    def add_rfm_opt_output(self, fldr, levels):
        """Load RFM output from files to the object as method.

        Args:
            fldr (str): Path to RFM folder.
            levels (array-like): Specifies levels at which to load output.

        Returns:
            Loads RFM results to the object, updates status.

        """
        self.rfm_output = rf.get_rfm_optical_depths(fldr=fldr, levels=levels)
        self.status = "RFM run completed and result loaded."
        return

    def get_wnos_from_RFM(self):
        """Determine wavenumbers from RFM output.

        Raises:
            ValueError: Raised when rfm output is empty.

        """
        if self.rfm_output is not None:
            try:
                [
                    float(i[i.rfind("_") + 1 :])
                    for i in rfm_df.columns
                    if i.startswith("dOD")
                ]
            except:
                print("Could not generate wavenumbers from rfm.")
        else:
            raise ValueError("rfm_output is empty (a NoneType object).")
        return

    def load_output_prf(self, fldr):
        """ "Loads the output profile (default prf.asc) file from RFM.

        Args:
            fldr (str): Path to RFM.

        Returns:
            Assings output_prf.

        """
        self.output_prf = rf.read_output_prf(f"{fldr}/prf.asc")
        return

    def calc_col_dens_and_mass(self, species, M=None):
        r"""Determines total column density of a species in the atmosphere.

        Integrates the species mixing ratio throughout the atmosphere to obtain total
        column density in units [molecules m\ :math:`^{-2}` \] and Dobson units [DU].
        If molar masses are provided, also calculates column masses in units
        [g m\ :math:`^{-2}` \].

        Function logic:
            1. checks if object has output_prf loaded.
            2. If not, tries to load it from a default directory.
            3. If 1 and 2 fail, error is raised.
            4. Checks if species are present as key in output_prf (and if not, then
                tries to add ppmv and check again).
            5. Integrates the profile for a given species (if they exist), returns dict
            6. If molar masses are provided, also calculates column masses.

        Args:
            species (list of str): list of species. The code first checks if the species
                are present in the profile. Mind that the rfm profile uses ppmv.
                If a key is not found, e.g. *CO*, and attempt will be made to transform
                it to *CO [ppmv]*.
            M (list of floats, optional): list of species' molar masses. If None, column
                 masses are not calculated. Default is None. Assumed units
                 [g mol\ :math:`^{-1}` \].

        Returns:
            col_den (dict): Dictionary containing species and its column density, units
                [molecules m\ :math:`^{-2}` \].
            col_den_DU (dict): Dictionary containing species and its column density,
                units [DU] (Dobson units).
            col_mass (dict, optional): Dictionary containing species and its column
                mass, units [g m\ :math:`^{-2}` \].

        Raises:
            ValueError: Raised when requested species not found.

        """
        # load output_prf (profile to calculate column density from
        if hasattr(self, "output_prf"):
            pass
        else:
            try:
                self.output_prf = rf.read_output_prf(f"./srfm/RFM//prf.asc")
            except FileNotFoundError:
                try:
                    self.output_prf = rf.read_atm_file(f"./srfm/RFM/rfm_files/day.atm")
                except:
                    print(
                        """prf not loaded and not found in default directory.
                        Load profile first through model_RFM.load_output_prf()."""
                    )

        # if species is str, transform into list
        if isinstance(species, str):
            species = [species]

        # check if species in keys
        for i_s, s in enumerate(species):
            key_matches = []
            for key in self.output_prf.keys():
                if key.startswith(s) or key.startswith(s.lower()):
                    key_matches.append(key)
            if len(key_matches) == 0:
                raise ValueError(
                    f"""Requested specie {s} not in atmospheric
                    profile."""
                )
            elif len(key_matches) > 1:
                raise ValueError(
                    f"""Requested species {s} ambiguous. Found the 
                    following possible matches: {key_matches}. Please specify you 
                    specie better."""
                )
            elif len(key_matches) == 1:
                species[i_s] = key_matches[0]  # replace specie with matched key

        # transfom output_prf lists into arrays (for vectorization)
        for key in self.output_prf.keys():
            self.output_prf[key] = np.array(self.output_prf[key])

        # lyr bounds is an array with layer boundaries. The RFM profile is specified in levels and at each level, atmospheric compositino is given.
        # to integrate the amount of gas in the atmosphere we construct atmospheric layers whose bounds in the middle of two adjacent levels

        # get layer bounds
        lyr_bounds = [
            (self.output_prf["HGT [km]"][i - 1] + self.output_prf["HGT [km]"][i - 1])
            / 2
            for i, _ in enumerate(self.output_prf["HGT [km]"][1:])
        ]
        lyr_bounds.insert(0, self.output_prf["HGT [km]"][0])
        lyr_bounds.insert(-1, self.output_prf["HGT [km]"][-1])

        lyr_bounds = np.array(lyr_bounds)

        lyr_thick = np.diff(lyr_bounds) * 1e3  # layer thicknesses in m

        # calculate the total amount of species in the atmosphere
        # layer are assumed homogeneous
        # N_tot = p/(kb*T) atmoshperic number density
        N_tot = (
            self.output_prf["PRE [mb]"]
            * 1e2
            / (self.output_prf["TEM [K]"] * scipy.constants.k)
        )  # molecules m-3

        self.col_den = {}  # dictionary to store results in
        self.col_den_DU = {}  # dictionary to store DU results in
        for s in species:
            N_species = (
                N_tot * self.output_prf[s] * 1e-6
            )  # molecules m-3 of species in the layer
            col_den_tot = np.sum(
                N_species * lyr_thick
            )  # molecules m-2 in the atmosphere
            self.col_den[s] = col_den_tot.item()
            self.col_den_DU[s] = (
                col_den_tot.item() / 2.69e20
            )  # 1 DU = 2.69e20 molec m-2

        if M is not None:
            self.col_mass = {}
            for s in species:
                self.col_mass[s] = (
                    self.col_den[s] / scipy.constants.N_A * M[species.index(s)]
                )

        return


class DISORT(Fwd_model):
    """Contain DISORT configuration, current output, and optional history.

    Args:
        name (str): Human-readable model name.
        disort_fldr (path-like | None): Optional DISORT working directory.
        disort_input (dict | None): Initial input dictionary; a fresh dictionary
            is created when omitted.
        disort_out (dict | None): Initial historical output dictionary; a fresh
            dictionary is created when omitted.
        retain_history (bool): Add every returned result to ``disort_out`` when
            true. Low-memory spectrum runners set this to false.
        disort_fmt_passmark (bool): Initial format-validation state.
        disort_integrity_passmark (bool): Initial integrity-validation state.
        status (str): Initial status text.
        **parameters: Additional attributes assigned to the instance.
    """

    def __init__(
        self,
        name="DISORT",
        disort_fldr=None,
        disort_input=None,
        disort_out=None,
        retain_history=True,
        disort_fmt_passmark=True,
        disort_integrity_passmark=True,
        status="DISORT object created.",
        **parameters,
    ):
        super().__init__(name)
        self.disort_fldr = disort_fldr
        self.disort_input = {} if disort_input is None else disort_input
        self.disort_out = {} if disort_out is None else disort_out
        self.retain_history = retain_history
        self.current_output = None
        self.disort_fmt_passmark = disort_fmt_passmark
        self.disort_integrity_passmark = disort_integrity_passmark
        self.status = status
        for key, val in parameters.items():
            setattr(self, key, val)

    def add_disort_input(self, d):
        """Add disort input parameters as a dictionary.

        Args:
            d (dict): Dictionary with DISORT input parameters.

        """
        self.disort_in = d

    def add_disort_empty_input(self):
        """Adds a dictionary with default and zero values as DISORT input.

        The dictionary contains zeros where possible and minimum or basic values in
        other instances (such as maxcly, which is 1 for one computational layer - the
        simplest possible case), etc.
        """

        maxcly = 1
        maxmom = 2
        maxcmu = 2
        maxumu = 1
        maxphi = 1
        maxulv = 1
        self.disort_input = {
            "maxcly": 1,
            "maxmom": 2,
            "maxcmu": 2,
            "maxumu": 1,
            "maxphi": 1,
            "maxulv": 1,
            "usrang": True,
            "usrtau": True,
            "ibcnd": 0,
            "onlyfl": False,
            "prnt": [True, True, True, False, True],  # print flags
            "plank": False,
            "lamber": True,
            "deltamplus": False,
            "do_pseudo_sphere": False,
            "dtauc": [0] * maxcly,
            "ssalb": [0] * maxcly,
            "pmom": np.zeros(shape=(maxmom + 1, maxcly)),
            "temper": np.zeros(shape=(maxcly + 1)),
            "wvnmlo": 0,
            "wvnmhi": 0,
            "utau": [0],
            "umu0": 0.1,
            "phi0": 0,
            "umu": [1],
            "phi": [0],
            "fbeam": 0,
            "fisot": 0,
            "albedo": 0,
            "btemp": 0,
            "ttemp": 0,
            "temis": 0,
            "earth_radius": 6371,
            "h_lyr": np.zeros(shape=(maxcly + 1)),
            "rhoq": np.zeros(shape=(int(maxcmu / 2), int(maxcmu / 2 + 1), int(maxcmu))),
            "rhou": np.zeros(shape=(maxumu, int(maxcmu / 2 + 1), maxcmu)),
            "rho_accurate": np.zeros(shape=(maxumu, maxphi)),
            "bemst": np.zeros(shape=(int(maxcmu / 2))),
            "emust": np.zeros(shape=(maxumu)),
            "accur": 0,
            "header": "",
            "rfldir": np.zeros(shape=(maxulv)),
            "rfldn": np.zeros(shape=(maxulv)),
            "flup": np.zeros(shape=(maxulv)),
            "dfdt": np.zeros(shape=(maxulv)),
            "uavg": np.zeros(shape=(maxulv)),
            "uu": np.zeros(shape=(maxumu, maxulv, maxphi)),
            "albmed": np.zeros(shape=(maxumu)),
            "trnmed": np.zeros(shape=(maxumu)),
        }

    def test_disort_input_format(self):
        """Test disort input for correct formatting."""
        if self.disort_input != {}:
            self.disort_fmt_passmark = disf.test_disort_input_format(
                maxcly=self.disort_input["maxcly"],
                maxmom=self.disort_input["maxmom"],
                maxcmu=self.disort_input["maxcmu"],
                maxumu=self.disort_input["maxumu"],
                maxphi=self.disort_input["maxphi"],
                maxulv=self.disort_input["maxulv"],
                usrang=self.disort_input["usrang"],
                usrtau=self.disort_input["usrtau"],
                ibcnd=self.disort_input["ibcnd"],
                onlyfl=self.disort_input["onlyfl"],
                prnt=self.disort_input["prnt"],
                plank=self.disort_input["plank"],
                lamber=self.disort_input["lamber"],
                deltamplus=self.disort_input["deltamplus"],
                do_pseudo_sphere=self.disort_input["do_pseudo_sphere"],
                dtauc=self.disort_input["dtauc"],
                ssalb=self.disort_input["ssalb"],
                pmom=self.disort_input["pmom"],
                temper=self.disort_input["temper"],
                wvnmlo=self.disort_input["wvnmlo"],
                wvnmhi=self.disort_input["wvnmhi"],
                utau=self.disort_input["utau"],
                umu0=self.disort_input["umu0"],
                phi0=self.disort_input["phi0"],
                umu=self.disort_input["umu"],
                phi=self.disort_input["phi"],
                fbeam=self.disort_input["fbeam"],
                fisot=self.disort_input["fisot"],
                albedo=self.disort_input["albedo"],
                btemp=self.disort_input["btemp"],
                ttemp=self.disort_input["ttemp"],
                temis=self.disort_input["temis"],
                earth_radius=self.disort_input["earth_radius"],
                h_lyr=self.disort_input["h_lyr"],
                rhoq=self.disort_input["rhoq"],
                rhou=self.disort_input["rhou"],
                rho_accurate=self.disort_input["rho_accurate"],
                bemst=self.disort_input["bemst"],
                emust=self.disort_input["emust"],
                accur=self.disort_input["accur"],
                header=self.disort_input["header"],
                rfldir=self.disort_input["rfldir"],
                rfldn=self.disort_input["rfldn"],
                flup=self.disort_input["flup"],
                dfdt=self.disort_input["dfdt"],
                uavg=self.disort_input["uavg"],
                uu=self.disort_input["uu"],
                albmed=self.disort_input["albmed"],
                trnmed=self.disort_input["trnmed"],
            )
        else:
            raise ValueError("DISORT.disort_input is not defined.")

    def test_disort_input_integrity(self):
        """Tests if DISORT input satisfies basic logical integrity and code constraints."""
        if self.disort_fmt_passmark == True:
            self.disort_integrity_passmark = disf.test_disort_input_integrity(
                maxmom=self.disort_input["maxmom"],
                maxcmu=self.disort_input["maxcmu"],
                maxumu=self.disort_input["maxumu"],
                maxphi=self.disort_input["maxphi"],
                ibcnd=self.disort_input["ibcnd"],
                onlyfl=self.disort_input["onlyfl"],
                dtauc=self.disort_input["dtauc"],
                ssalb=self.disort_input["ssalb"],
                temper=self.disort_input["temper"],
                wvnmlo=self.disort_input["wvnmlo"],
                wvnmhi=self.disort_input["wvnmhi"],
                utau=self.disort_input["utau"],
                umu0=self.disort_input["umu0"],
                phi0=self.disort_input["phi0"],
                umu=self.disort_input["umu"],
                phi=self.disort_input["phi"],
                btemp=self.disort_input["btemp"],
                ttemp=self.disort_input["ttemp"],
                temis=self.disort_input["temis"],
            )
        else:
            raise ValueError(
                "disort_format_passmark == False, the input has not passed the format test."
            )

    def set_maxcly(self, maxcly):
        """Assigns the value of maxcly.

        Args:
            maxcly (int): Number of computational layers.

        """
        self.disort_input["maxcly"] = maxcly

    def set_maxcly_from_rfm(self, rfm):
        """Sets maxcly from RFM.

        Args:
            rfm (obj): Class RFM object and has attribute rfm_output.

        Returns:
            None: ``maxcly`` is updated in ``disort_input``.

        Raises:
            AttributeError: Raised when rfm does not have rfm_output attribute.
        """
        try:
            output = rfm.rfm_output
            self.disort_input["maxcly"] = (
                output.layer_count
                if hasattr(output, "layer_count")
                else len(output["layer no."])
            )
        except (AttributeError, KeyError) as exc:
            raise AttributeError(
                "RFM must contain a populated rfm_output attribute."
            ) from exc

    def set_maxmom(self, maxmom):
        """Assigns the value of maxmom.

        Args:
            maxmom (int): Number of phase function moments.

        """
        self.disort_input["maxmom"] = maxmom

    def set_maxcmu(self, maxcmu):
        """Assigns the value of maxcmu.

        Args:
            maxcmu (int): Number of computational streams.

        """
        self.disort_input["maxcmu"] = maxcmu

    def set_maxumu(self, maxumu):
        """Assigns the value of maxumu.

        Args:
            maxumu (int): Number of output polar angles.

        """
        self.disort_input["maxumu"] = maxumu

    def set_maxphi(self, maxphi):
        """Assigns the value of maxphi.

        Args:
            maxphi (int): Number of output azimuthal angles.

        """
        self.disort_input["maxphi"] = maxphi

    def set_maxulv(self, maxulv):
        """Assigns the value of maxulv.

        Args:
            maxulv (int): Number of output optical depths.

        """
        self.disort_input["maxulv"] = maxulv

    def set_usrang(self, usrang):
        """Assigns the value of usrang.

        Args:
            usrang (bool): If True, output at user-defined polar angles requested.

        """
        self.disort_input["usrang"] = usrang

    def set_usrtau(self, usrtau):
        """Assigns the value of usrtau.

        Args:
            usrtau (bool): If True, output at user-defined optical depths requested.

        """
        self.disort_input["usrtau"] = usrtau

    def set_ibcnd(self, ibcnd):
        """Assigns the value of ibcnd.

        Args:
            ibcnd (int): Specifices a combination of boundary conditions. For full
                documentation see DISORT documentation.

        """
        self.disort_input["ibcnd"] = ibcnd

    def set_onlyfl(self, onlyfl):
        """Assigns the value of onlyfl.

        Args:
            onlyfl (bool): If True, only fluxes output, elif False, full output.

        """
        self.disort_input["onlyfl"] = onlyfl

    def set_prnt(self, prnt):
        """Assigns the value of prnt.

        Args:
            prnt (array-like): Array of shape (5,) bool values, controls printing to terminal
                from DISORT.

        """
        self.disort_input["prnt"] = prnt

    def set_plank(self, planck):
        """Assigns the value of plank.

        Args:
            planck (bool): Thermal radiation in DISORT. True - on, False - off.

        """
        self.disort_input["plank"] = planck

    def set_lamber(self, lamber):
        """Assigns the value of lamber.

        Args:
            lamber (bool): Bottom boundary in DISORT treated as Lambertian or not.

        """
        self.disort_input["lamber"] = lamber

    def set_deltamplus(self, deltamplus):
        """Assigns the value of deltamplus.

        Args:
            deltamplus: If True, use delta-M-plus approximation to calculate strongly
                forward-peaked phase functions. If False, use delta-M method.

        """
        self.disort_input["deltamplus"] = deltamplus

    def set_do_pseudo_sphere(self, d_p_s):
        """Assigns the value of do_pseudo_sphere.

        Args:
            d_p_s (bool): Do a spheric correction on the calculation.

        """
        self.disort_input["do_pseudo_sphere"] = d_p_s

    def set_dtauc_manually(self, dtauc):
        """Assigns the value of dtauc.

        Args:
            dtauc (array-like): Atmospheric optical depth structure (layers' optical
                depth).

        """
        self.disort_input["dtauc"] = dtauc

    def set_ssalb_manually(self, ssalb):
        """Assigns the value of ssalb.

        Args:
            ssalb (array-like): Layers' single scatter albedo.

        """
        self.disort_input["ssalb"] = ssalb

    def set_pmom_manually(self, pmom):
        """Assigns the value of pmom.

        Args:
            pmom (array-like): Scattering phase function Legendre polynomial expansion
                coefficients (normalised coefficients expected).

        """
        self.disort_input["pmom"] = pmom

    def set_temper(self, temper):
        """Assigns the value of temper.

        Args:
            temper (array-like): Atmospheric tempeature structure, defined in terms of
                levels.

        """
        self.disort_input["temper"] = temper

    def set_temper_from_rfm(self, rfm):
        """Sets layer temperature for disort from the RFM output.

        Takes upper temperatures from rfm_output for every layer.
        Adds the lowest layer lower temperature (adjacent to surface) as the last one.

        Args:
            rfm (obj): Class RFM object which has attribute rfm_output.

        Returns:
            None: ``temper`` is updated in ``disort_input``.

        Raises:
            AttributeError: If the RFM output has no temperature profile.

        """
        try:
            output = rfm.rfm_output
            if hasattr(output, "temperature_upper"):
                temperatures = output.temperature_upper.tolist()
                temperatures.append(float(output.temperature_lower[-1]))
            else:
                temperatures = output["T_upper (K)"].tolist()
                temperatures.append(output["T_lower (K)"].tolist()[-1])
            self.disort_input["temper"] = temperatures
        except (AttributeError, KeyError) as exc:
            raise AttributeError(
                "RFM must contain temperatures in rfm_output."
            ) from exc
        return

    def set_wvnm_range(self, lo, hi):
        """Assigns the value of wvnmlo and wvnmhi.

        Args:
            lo (int, float): Lower wavenumber for Planck function calculation.
            hi (int, float): Upper wavenumber for Planck function calculation.

        """
        self.disort_input["wvnmlo"] = lo
        self.disort_input["wvnmhi"] = hi

    def set_utau(self, utau):
        """Assigns the value of utau.

        Args:
             utau (array-like): User requested output optical depths.

        """
        self.disort_input["utau"] = utau

    def set_umu0(self, umu0):
        """Assigns the value of umu0.

        Args:
            umu0 (int, float): Incoming direct beam polar angle.

        """
        self.disort_input["umu0"] = umu0

    def set_phi0(self, phi0):
        """Assigns the value of phi0.

        Args:
            phi0 (int,float): Incoming direct beam azimuthal angle.

        """
        self.disort_input["phi0"] = phi0

    def set_umu(self, umu):
        """Assigns the value of umu.

        Args:
            umu (array-like): User requested output polar angles.

        """
        self.disort_input["umu"] = umu

    def set_phi(self, phi):
        """Assigns the value of phi.

        Args:
            phi (array-like): User requested output azimuthal angles.

        """
        self.disort_input["phi"] = phi

    def set_fbeam(self, fbeam):
        """Assigns the value of fbeam.

        Args:
            fbeam (int, float): Incoming direct beam intensity.

        """
        self.disort_input["fbeam"] = fbeam

    def set_fisot(self, fisot):
        """Assigns the value of fisot.

        Args:
            fisot (int, float): Incoming diffuse radiation intensity.

        """
        self.disort_input["fisot"] = fisot

    def set_albedo(self, albedo):
        """Assigns the value of albedo.

        Args:
            albedo (int, float): Surface albedo.

        """
        self.disort_input["albedo"] = albedo

    def set_btemp(self, btemp):
        """Assigns the value of btemp.

        Args:
            btemp (int, float): Bottom boundary temperature.

        """
        self.disort_input["btemp"] = btemp

    def set_ttemp(self, ttemp):
        """Assigns the value of ttemp.

        Args:
            ttemp (int, float): Top boundary temperature.

        """
        self.disort_input["ttemp"] = ttemp

    def set_temis(self, temis):
        """Assigns the value of temis.

        Args:
            temis (int, float): Top boundary emmisivity.

        """
        self.disort_input["temis"] = temis

    def set_earth_radius(self, earth_radius):
        """Assigns the value of earth_radius.

        Args:
            earth_radius (int, float): Earth radius, [km].

        """
        self.disort_input["earth_radius"] = earth_radius

    def set_h_lyr(self, h_lyr):
        """Assigns the value of h_lyr.

        Args:
            h_lyr (array-like): layer vertical extent.

        """
        self.disort_input["h_lyr"] = h_lyr

    def set_h_lyr_from_rfm(self, rfm):
        """Assigns the value of h_lyr from the RFM output.

        Args:
            rfm (obj): Class RFM object which has attribute rfm_output.

        Returns:
            None: ``h_lyr`` is updated in ``disort_input``.

        Raises:
            AttributeError: Raised when RFM object does not have rfm_output attribute.
        """

        try:
            output = rfm.rfm_output
            if hasattr(output, "altitude_upper"):
                heights = output.altitude_upper.tolist()
                heights.append(float(output.altitude_lower[-1]))
            else:
                heights = output["h_upper (km)"].tolist()
                heights.append(output["h_lower (km)"].tolist()[-1])
            self.disort_input["h_lyr"] = heights
        except (AttributeError, KeyError) as exc:
            raise AttributeError("RFM must contain heights in rfm_output.") from exc

    def set_rhoq(self, rhoq):
        """Assigns the value of rhoq.

        Args:
            rhoq (array-like): Something to do with BDREF in DISORT?.

        """
        self.disort_input["rhoq"] = rhoq

    def set_rhou(self, rhou):
        """Assigns the value of rhou.

        Args:
            rhou (array-like): Something to do with BDREF in DISORT?

        """

        self.disort_input["rhou"] = rhou

    def set_rho_accurate(self, rho_accurate):
        """Assigns the value of rho_accurate.

        Args:
            rhou (array-like): Something to do with BDREF in DISORT?

        """
        self.disort_input["rho_accurate"] = rho_accurate

    def set_bemst(self, bemst):
        """Assigns the value of bemst.

        Args:
            bemst (array-like): Something to do with BDREF in DISORT?

        """
        self.disort_input["bemst"] = bemst

    def set_emust(self, emust):
        """Assigns the value of emust.

        Args:
            emust (array-like): Something to do with BDREF in DISORT?

        """
        self.disort_input["emust"] = emust

    def set_accur(self, accur):
        """Assigns the value of accur.

        Args:
            accur (int, float): Convergence criterion for azimuthal (Fourier cosine)
                series.

        """
        self.disort_input["accur"] = accur

    def set_header(self, header):
        """Assigns the value of header.

        Args:
            header (str): Header for terminal output printing. The string "NO HEADER"
                will cause nothing to be printed out. Warning: If the input is an empty
                string (""), a blank line will be printed. The string must have
                len < 127.

        """
        self.disort_input["header"] = header

    def initialize_disort_output_arrays(self):
        """Initializes empty output arrays for a single disort run."""
        self.disort_input["rfldir"] = np.zeros(self.disort_input["maxulv"])
        self.disort_input["rfldn"] = np.zeros(self.disort_input["maxulv"])
        self.disort_input["flup"] = np.zeros(self.disort_input["maxulv"])
        self.disort_input["dfdt"] = np.zeros(self.disort_input["maxulv"])
        self.disort_input["uavg"] = np.zeros(self.disort_input["maxulv"])
        self.disort_input["uu"] = np.zeros(
            (
                self.disort_input["maxumu"],
                self.disort_input["maxulv"],
                self.disort_input["maxphi"],
            )
        )
        self.disort_input["albmed"] = np.zeros(self.disort_input["maxumu"])
        self.disort_input["trnmed"] = np.zeros(self.disort_input["maxumu"])

    def run_disort(self, prec="double", adjust_maxcmu=True):
        """Run DISORT at the requested precision and return its current output.

        Args:
            prec (str): Determines Fortran precision to be used (single vs double).
                Default is double precision calculation.
            adjust_maxcmu (bool): If True, check the output intensity in the downward
                direction. If that turns out negative and close to zero, it's likely
                caused by roundoff errors in the Delta-M+ algorithm. Increasing the
                number of computational streams usually fixes it, so in case this
                happens, the number of streams is automatically adjusted and DISORT run
                again.

        Returns:
            DisortResult: Normalized outputs for the current spectral point.

        Raises:
            ValueError: If ``prec`` is neither ``"single"`` nor ``"double"``.
        """
        if prec == "double":
            result = self.run_disort_double(adjust_maxcmu)
        elif prec == "single":
            result = self.run_disort_single(adjust_maxcmu)
        else:
            raise ValueError("prec must be 'single' or 'double'.")
        return result

    def _record_result(self, native_result):
        """Normalize and optionally retain one native DISORT result.

        Args:
            native_result (tuple): Values returned by the f2py DISORT wrapper.

        Returns:
            DisortResult: Structured output for the current spectral point.

        Raises:
            ValueError: If the configured wavenumber interval has zero width.
        """
        interval = self.disort_input["wvnmhi"] - self.disort_input["wvnmlo"]
        if interval == 0:
            raise ValueError("DISORT wavenumber interval must have non-zero width.")
        result = DisortResult(
            wavenumber=self.wvnm,
            wavelength=self.wvl,
            rfldir=native_result[0] / interval,
            rfldn=native_result[1] / interval,
            flup=native_result[2] / interval,
            dfdt=native_result[3] / interval,
            uavg=native_result[4] / interval,
            uu=native_result[5] / interval,
            albmed=native_result[6],
            trnmed=native_result[7] if len(native_result) == 8 else 0,
        )
        self.current_output = result
        if self.retain_history:
            self.disort_out[self.wvnm] = result.as_dict()
        self.status = "DISORT run completed."
        return result

    def set_wvnm(self, wvnm):
        """Sets wavenumber of the current run.

        Args:
            wvnm (int, float): wavenumber [cm-1].

        """
        self.wvnm = wvnm

    def set_wvl(self, wvl):
        """Sets wavelength of the current run.

        Args:
            wvl (int, float): wavelength [cm-1].

        """
        self.wvl = wvl

    def run_disort_single(self, adjust_maxcmu):
        """Run the single-precision DISORT wrapper.

        Args:
            adjust_maxcmu (bool): Retry small negative radiances with more streams.

        Returns:
            DisortResult: Normalized outputs for the current spectral point.
        """
        # run DISORT
        res = dms.disort(
            maxcly=self.disort_input["maxcly"],
            maxmom=self.disort_input["maxmom"],
            maxcmu=self.disort_input["maxcmu"],
            maxumu=self.disort_input["maxumu"],
            maxphi=self.disort_input["maxphi"],
            maxulv=self.disort_input["maxulv"],
            usrang=self.disort_input["usrang"],
            usrtau=self.disort_input["usrtau"],
            ibcnd=self.disort_input["ibcnd"],
            onlyfl=self.disort_input["onlyfl"],
            prnt=self.disort_input["prnt"],
            plank=self.disort_input["plank"],
            lamber=self.disort_input["lamber"],
            deltamplus=self.disort_input["deltamplus"],
            do_pseudo_sphere=self.disort_input["do_pseudo_sphere"],
            dtauc=self.disort_input["dtauc"],
            ssalb=self.disort_input["ssalb"],
            pmom=self.disort_input["pmom"],
            temper=self.disort_input["temper"],
            wvnmlo=self.disort_input["wvnmlo"],
            wvnmhi=self.disort_input["wvnmhi"],
            utau=self.disort_input["utau"],
            umu0=self.disort_input["umu0"],
            phi0=self.disort_input["phi0"],
            umu=self.disort_input["umu"],
            phi=self.disort_input["phi"],
            fbeam=self.disort_input["fbeam"],
            fisot=self.disort_input["fisot"],
            albedo=self.disort_input["albedo"],
            btemp=self.disort_input["btemp"],
            ttemp=self.disort_input["ttemp"],
            temis=self.disort_input["temis"],
            earth_radius=self.disort_input["earth_radius"],
            h_lyr=self.disort_input["h_lyr"],
            rhoq=self.disort_input["rhoq"],
            rhou=self.disort_input["rhou"],
            rho_accurate=self.disort_input["rho_accurate"],
            bemst=self.disort_input["bemst"],
            emust=self.disort_input["emust"],
            accur=self.disort_input["accur"],
            header=self.disort_input["header"],
            rfldir=self.disort_input["rfldir"],
            rfldn=self.disort_input["rfldn"],
            flup=self.disort_input["flup"],
            dfdt=self.disort_input["dfdt"],
            uavg=self.disort_input["uavg"],
            uu=self.disort_input["uu"],
            albmed=self.disort_input["albmed"],
            trnmed=self.disort_input["trnmed"],
        )

        #        res[5][res[5]<1e-6] = 0

        if adjust_maxcmu:
            # check if intensity is negative and potentially rerun DISORT
            # use carefully as there may be cases where intensity is negative (see DISORT docs)
            maxcmu_multiplier = 1
            old_maxcmu = self.disort_input["maxcmu"]  # save old maxcmu
            old_maxmom = self.disort_input["maxmom"]  # save old maxmom
            while not np.all(
                res[5] / (self.disort_input["wvnmhi"] - self.disort_input["wvnmlo"]) > 0
            ):
                self.disort_input["maxcmu"] = old_maxcmu * 2 * maxcmu_multiplier
                if self.disort_input["maxmom"] < self.disort_input["maxcmu"]:
                    self.disort_input["maxmom"] = self.disort_input["maxcmu"]
                self.disort_input["rhoq"] = np.zeros(
                    shape=(
                        int(self.disort_input["maxcmu"] / 2),
                        int(self.disort_input["maxcmu"] / 2 + 1),
                        int(self.disort_input["maxcmu"]),
                    )
                )
                self.disort_input["rhou"] = np.zeros(
                    shape=(
                        self.disort_input["maxumu"],
                        int(self.disort_input["maxcmu"] / 2 + 1),
                        self.disort_input["maxcmu"],
                    )
                )
                self.disort_input["bemst"] = np.zeros(
                    shape=(int(self.disort_input["maxcmu"] / 2))
                )
                res = dms.disort(
                    maxcly=self.disort_input["maxcly"],
                    maxmom=self.disort_input["maxmom"],
                    maxcmu=self.disort_input["maxcmu"],
                    maxumu=self.disort_input["maxumu"],
                    maxphi=self.disort_input["maxphi"],
                    maxulv=self.disort_input["maxulv"],
                    usrang=self.disort_input["usrang"],
                    usrtau=self.disort_input["usrtau"],
                    ibcnd=self.disort_input["ibcnd"],
                    onlyfl=self.disort_input["onlyfl"],
                    prnt=self.disort_input["prnt"],
                    plank=self.disort_input["plank"],
                    lamber=self.disort_input["lamber"],
                    deltamplus=self.disort_input["deltamplus"],
                    do_pseudo_sphere=self.disort_input["do_pseudo_sphere"],
                    dtauc=self.disort_input["dtauc"],
                    ssalb=self.disort_input["ssalb"],
                    pmom=self.disort_input["pmom"],
                    temper=self.disort_input["temper"],
                    wvnmlo=self.disort_input["wvnmlo"],
                    wvnmhi=self.disort_input["wvnmhi"],
                    utau=self.disort_input["utau"],
                    umu0=self.disort_input["umu0"],
                    phi0=self.disort_input["phi0"],
                    umu=self.disort_input["umu"],
                    phi=self.disort_input["phi"],
                    fbeam=self.disort_input["fbeam"],
                    fisot=self.disort_input["fisot"],
                    albedo=self.disort_input["albedo"],
                    btemp=self.disort_input["btemp"],
                    ttemp=self.disort_input["ttemp"],
                    temis=self.disort_input["temis"],
                    earth_radius=self.disort_input["earth_radius"],
                    h_lyr=self.disort_input["h_lyr"],
                    rhoq=self.disort_input["rhoq"],
                    rhou=self.disort_input["rhou"],
                    rho_accurate=self.disort_input["rho_accurate"],
                    bemst=self.disort_input["bemst"],
                    emust=self.disort_input["emust"],
                    accur=self.disort_input["accur"],
                    header=self.disort_input["header"],
                    rfldir=self.disort_input["rfldir"],
                    rfldn=self.disort_input["rfldn"],
                    flup=self.disort_input["flup"],
                    dfdt=self.disort_input["dfdt"],
                    uavg=self.disort_input["uavg"],
                    uu=self.disort_input["uu"],
                    albmed=self.disort_input["albmed"],
                    trnmed=self.disort_input["trnmed"],
                )

                maxcmu_multiplier += 1
                if self.disort_input["maxcmu"] > 128:
                    self.disort_input["maxcmu"] = old_maxcmu
                    self.disort_input["maxmom"] = old_maxmom
                    self.disort_input["rhoq"] = np.zeros(
                        shape=(
                            int(self.disort_input["maxcmu"] / 2),
                            int(self.disort_input["maxcmu"] / 2 + 1),
                            int(self.disort_input["maxcmu"]),
                        )
                    )
                    self.disort_input["rhou"] = np.zeros(
                        shape=(
                            self.disort_input["maxumu"],
                            int(self.disort_input["maxcmu"] / 2 + 1),
                            self.disort_input["maxcmu"],
                        )
                    )
                    self.disort_input["bemst"] = np.zeros(
                        shape=(int(self.disort_input["maxcmu"] / 2))
                    )
                    print(
                        f"maxcmu increased to {old_maxcmu * 2 * maxcmu_multiplier} and intensity still negative, skipping."
                    )
                    break

        return self._record_result(res)

    def run_disort_double(self, adjust_maxcmu):
        """Run the double-precision DISORT wrapper.

        Args:
            adjust_maxcmu (bool): Retry small negative radiances with more streams.

        Returns:
            DisortResult: Normalized outputs for the current spectral point.
        """
        res = dmd.disort(
            maxcly=self.disort_input["maxcly"],
            maxmom=self.disort_input["maxmom"],
            maxcmu=self.disort_input["maxcmu"],
            maxumu=self.disort_input["maxumu"],
            maxphi=self.disort_input["maxphi"],
            maxulv=self.disort_input["maxulv"],
            usrang=self.disort_input["usrang"],
            usrtau=self.disort_input["usrtau"],
            ibcnd=self.disort_input["ibcnd"],
            onlyfl=self.disort_input["onlyfl"],
            prnt=self.disort_input["prnt"],
            plank=self.disort_input["plank"],
            lamber=self.disort_input["lamber"],
            deltamplus=self.disort_input["deltamplus"],
            do_pseudo_sphere=self.disort_input["do_pseudo_sphere"],
            dtauc=self.disort_input["dtauc"],
            ssalb=self.disort_input["ssalb"],
            pmom=self.disort_input["pmom"],
            temper=self.disort_input["temper"],
            wvnmlo=self.disort_input["wvnmlo"],
            wvnmhi=self.disort_input["wvnmhi"],
            utau=self.disort_input["utau"],
            umu0=self.disort_input["umu0"],
            phi0=self.disort_input["phi0"],
            umu=self.disort_input["umu"],
            phi=self.disort_input["phi"],
            fbeam=self.disort_input["fbeam"],
            fisot=self.disort_input["fisot"],
            albedo=self.disort_input["albedo"],
            btemp=self.disort_input["btemp"],
            ttemp=self.disort_input["ttemp"],
            temis=self.disort_input["temis"],
            earth_radius=self.disort_input["earth_radius"],
            h_lyr=self.disort_input["h_lyr"],
            rhoq=self.disort_input["rhoq"],
            rhou=self.disort_input["rhou"],
            rho_accurate=self.disort_input["rho_accurate"],
            bemst=self.disort_input["bemst"],
            emust=self.disort_input["emust"],
            accur=self.disort_input["accur"],
            header=self.disort_input["header"],
            rfldir=self.disort_input["rfldir"],
            rfldn=self.disort_input["rfldn"],
            flup=self.disort_input["flup"],
            dfdt=self.disort_input["dfdt"],
            uavg=self.disort_input["uavg"],
            uu=self.disort_input["uu"],
            albmed=self.disort_input["albmed"],
            trnmed=self.disort_input["trnmed"],
        )

        return self._record_result(res)

    def calc_bbt(self):
        """Converts radiance to brightness temperature."""
        for key in self.disort_out.keys():
            self.disort_out[key]["uu_bbt"] = utils.convert_spectral_radiance_to_bbt(
                self.disort_out[key]["uu"], self.disort_out[key]["wavenumber (cm-1)"]
            )

    def save_model_pickle(self, filename=None, folder=None):
        """Saves the model with to a file with pickle.

        Args:
            filename (str): Name of file to be saved. If None, then default name
                "model.pkl" is used. Default is None.
            folder (str): Location for the file to be saved to. If None, then save in
                the current folder. Default is None.

        """
        if isinstance(filename, (NoneType, str)):
            if filename == None:
                fl = "model"
            else:
                fl = filename
        else:
            raise TypeError("filename must be None or str.")

        if isinstance(folder, (NoneType, str)):
            if folder == None:
                fldr = "."
            else:
                fldr = folder
        else:
            raise TypeError("folder must be None or str.")

        with open(f"{fldr}/{fl}.pkl", "wb") as f:
            pickle.dump(self, f)
        f.close()
        return

    def set_dtauc(self, tau_g, tau_R, tau_p):
        """Calculate optical depth of model layers, delta tau (dtau).

        Calculation according to the formula :math:`dtau = tau_g + tau_R + tau_p`
        Value is assigned to the object.

        Args:
            tau_g (array-like): layer optical depth from gas absorption
            tau_R (array-like): layer optical depth from Rayleigh scattering
            tau_p (array-like): layer optical depth from particle scattering

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
                    f"inputs must be np.ndarrays, pd.Series, lists, ints or floats."
                )

            return obj

        tau_g = check_convert_dtype(tau_g)
        tau_R = check_convert_dtype(tau_R)
        tau_p = check_convert_dtype(tau_p)

        self.disort_input["dtauc"] = tau_g + tau_R + tau_p
        return

    def set_ssalb(self, tau_g, tau_R, tau_p, w_p):
        """Calculates and adds the single scatter albedo of the model layers.

        The calculation formula is

        .. math::

            \\begin{eqnarray}
                w = \\frac{w_g*\\tau_g + w_R*\\tau_R + w_p*\\tau_p}
                {\\tau_g + \\tau_R + \\tau_p}
            \\end{eqnarray}

        where:
            - w: layer single scatter albedo
            - w_g: layer single scatter albedo from gas absorption, == 0
            - w_R: layer Rayleigh scattering single scatter albedo, == 1
            - w_p: layer particle scattering single scatter albedo
            - tau_g: layer gas absorption optical depth
            - tau_R: layer Rayleigh scattering optical depth
            - tau_p: layer particle scattering optical depth

        Args:
            tau_g (array-like): layer optical depth from gas absorption
            tau_R (array-like): layer optical depth from Rayleigh scattering
            tau_p (array-like): layer optical depth from particle scattering
            w_p (array_like): layer particle scattering single scatter albedo

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
                    f"inputs must be np.ndarrays, pd.Series, lists, ints or floats."
                )

            return obj

        tau_g = check_convert_dtype(tau_g)
        tau_R = check_convert_dtype(tau_R)
        tau_p = check_convert_dtype(tau_p)
        w_p = check_convert_dtype(w_p)

        ssalb = np.nan_to_num(((tau_R + w_p * tau_p) / (tau_R + tau_g + tau_p)))

        self.disort_input["ssalb"] = ssalb
        return

    def set_pmom(self, pmom_R, tau_R, w_p, tau_p, pmom_p):
        """Calculate phase function coefficients according to Don (from ORAC).

        The calculation formula is

        .. math::

            \\begin{eqnarray}
                x_i = \\frac{w_g*\\tau_g*x_{i,g} + w_R*\\tau_R*x_{i,r} +
                w_p*\\tau_p*x_{i,p}}{w_g*\\tau_g + w_r*\\tau_R + w_p*\\tau_p}
            \\end{eqnarray}

        where:
            - :math:`x_i`: Legendre polynomial coefficient
            - :math:`w_g`: = 0 - gas single scatter albedo
            - :math:`x_{i,g}`: = 0, gas phase function moment
            - :math:`w_R`: = 1, Rayleigh scattering single scatter albedo

        which makes the formula

        .. math::

            \\begin{eqnarray}
                x_i = \\frac{tau_R*x_{i,r} + w_p*\\tau_p*x_{i,p}}{tau_R + w_p*\\tau_p}
            \\end{eqnarray}

        The function works on arrays, so :math:`x_i` is replaced by array pmom.

        Args:
            pmom_R (array-like): Array of phase function moments for Rayleigh
                scattering. 2D array, where *columns* are the atmospheric layers
                *rows* are the Legendre polynomial coefficients :math:`x_{1,R}` to
                :math:`x_{n,R}`. Has shape (model_DISORT.disort_input["maxmom"] + 1,
                model_DISORT.disort_input["maxcly"]).
            pmom_p (array-like): Array of phase function coefficients for particle
                scattering. Same shape and meaning as pmom_R.
            w_p (array-like): Particle scattering single scatter albedo for each layer.
            tau_p (array-like): Particle scattering optical depth for each layer.
            tau_R (array-like): Rayleigh optical depths.

        Raises:
            ValueError: When inputs are incorrectly shaped.
            TypeError: When inputs are not np.ndarrays.

        """

        maxmom = self.disort_input["maxmom"]
        maxcly = self.disort_input["maxcly"]

        pmom_R = np.asarray(pmom_R)
        pmom_p = np.asarray(pmom_p)
        tau_R = np.asarray(tau_R, dtype=float)
        tau_p = np.asarray(tau_p, dtype=float)
        w_p = np.asarray(w_p, dtype=float)

        expected = (maxmom + 1, maxcly)
        if pmom_R.shape != expected or pmom_p.shape != expected:
            raise ValueError("pmom_p and pmom_R must have shape (maxmom+1, maxcly)")
        if (
            tau_R.shape != (maxcly,)
            or tau_p.shape != (maxcly,)
            or w_p.shape != (maxcly,)
        ):
            raise ValueError("tau_R, tau_p, w_p must have length maxcly")

        # Broadcast 1D optical-depth vectors across the Legendre dimension.
        tau_R_2d = tau_R[np.newaxis, :]
        tau_p_2d = tau_p[np.newaxis, :]
        w_p_2d = w_p[np.newaxis, :]

        # precompute numerator and denominator
        numer = tau_R_2d * pmom_R + (w_p_2d * tau_p_2d) * pmom_p
        denom = tau_R_2d + w_p_2d * tau_p_2d

        pmom = np.zeros_like(numer)
        pmom = np.divide(numer, denom, out=pmom, where=denom != 0.0)

        self.disort_input["pmom"] = pmom
        return

    def set_mixed_pmom(
        self,
        tau_R,
        w_p,
        tau_p,
        particle_moments=None,
        prec="double",
    ):
        """Build the mixed Rayleigh-particle moments in one reusable workspace.

        Rayleigh scattering has coefficients 1.0 and 0.1 at indices zero and
        two, respectively, with all other coefficients equal to zero. Particle
        coefficients are supplied only for layers that contain particles, avoiding
        the two mostly redundant dense input matrices required by :meth:`set_pmom`.

        Args:
            tau_R (array-like): Rayleigh optical depth for each retained layer.
            w_p (array-like): Particle single-scattering albedo for each layer.
            tau_p (array-like): Particle optical depth for each retained layer.
            particle_moments (Mapping[int, array-like] | None): Particle Legendre
                vectors keyed by retained atmospheric-layer index.
            prec (str): DISORT precision, either ``"single"`` or ``"double"``.

        Returns:
            numpy.ndarray: C-contiguous moment workspace passed to DISORT. The
                bundled f2py wrapper otherwise misinterprets production-sized
                Fortran-contiguous inputs, so its unavoidable conversion is left
                explicit and documented.

        Raises:
            ValueError: If precision, vector lengths, or particle-layer indices are
                invalid.
        """
        if prec not in {"single", "double"}:
            raise ValueError("prec must be 'single' or 'double'.")
        dtype = np.float32 if prec == "single" else np.float64
        maxmom = self.disort_input["maxmom"]
        maxcly = self.disort_input["maxcly"]
        tau_R = np.asarray(tau_R, dtype=dtype)
        tau_p = np.asarray(tau_p, dtype=dtype)
        w_p = np.asarray(w_p, dtype=dtype)
        if any(vector.shape != (maxcly,) for vector in (tau_R, tau_p, w_p)):
            raise ValueError("tau_R, tau_p, w_p must have length maxcly")

        workspace_key = (dtype, maxmom + 1, maxcly)
        workspaces = getattr(self, "_pmom_workspaces", {})
        pmom = workspaces.get(workspace_key)
        if pmom is None:
            # The bundled f2py DISORT interface interprets a directly supplied
            # Fortran-contiguous array incorrectly for production-sized PMOM inputs.
            # C order triggers its safe input conversion and is therefore required.
            pmom = np.empty((maxmom + 1, maxcly), dtype=dtype, order="C")
            workspaces[workspace_key] = pmom
            self._pmom_workspaces = workspaces
        pmom.fill(0)

        denominator = tau_R + w_p * tau_p
        nonzero = denominator != 0
        pmom[0, nonzero] = tau_R[nonzero] / denominator[nonzero]
        if maxmom >= 2:
            pmom[2, nonzero] = 0.1 * tau_R[nonzero] / denominator[nonzero]

        for layer_index, coefficients in (particle_moments or {}).items():
            if layer_index < 0 or layer_index >= maxcly:
                raise ValueError("Particle moment layer index is outside maxcly.")
            if not nonzero[layer_index]:
                continue
            coefficient_array = np.asarray(coefficients, dtype=dtype)
            coefficient_count = min(coefficient_array.size, maxmom + 1)
            particle_weight = (
                w_p[layer_index] * tau_p[layer_index] / denominator[layer_index]
            )
            pmom[:coefficient_count, layer_index] += (
                particle_weight * coefficient_array[:coefficient_count]
            )

        # Avoid a one-ulp overshoot from adding separately weighted terms. DISORT
        # requires the zeroth moment to be exactly normalized and rejects values
        # even infinitesimally above one.
        pmom[0, nonzero] = 1.0

        self.disort_input["pmom"] = pmom
        return pmom

    def calc_pmom(self, iphas, prec="double", gg=0):
        """Calculates phase function moments from disort using the getmom function.

        Args:
            iphas (int): phase function option. Can be:
                - 1: Isotropic
                - 2: Rayleigh
                - 3: Henyey-Greenstein with asymmetry factor GG
                - 4: Haze L as specified by Garcia/Siewert
                - 5: Cloud C.1 as specified by Garcia/Siewert
                - 6: Aerosol as specified by Kokhanovsky
                - 7: Cloud as specified by Kokhanovsky
            gg (int,float): Assymetry factor for Heyney-Greenstein case. Default is 0.
            prec (str): Required precision mode, accepted values "single" or "double".
                Default is double.

        Raises:
            ValueError: When invalid precision values is used.

        """

        if prec == "single":
            pmom = self.calc_pmom_single(iphas, gg)
        elif prec == "double":
            pmom = self.calc_pmom_double(iphas, gg)
        else:
            raise ValueError("prec must be 'single' or 'double'.")
        return pmom

    def calc_pmom_single(self, iphas, gg=0):
        """Calculates phase function moments from disort using the getmom function.

        This is a class method because it contains a loop which is inconvenient
        in the main code (getmom calcualtes phase function moments
        for a 1D array/list, not a 2D array.
        This function uses the single precision version of DISORT.

        Args:
            iphas (int): phase function option. Can be:
                - 1: Isotropic
                - 2: Rayleigh
                - 3: Henyey-Greenstein with asymmetry factor GG
                - 4: Haze L as specified by Garcia/Siewert
                - 5: Cloud C.1 as specified by Garcia/Siewert
                - 6: Aerosol as specified by Kokhanovsky
                - 7: Cloud as specified by Kokhanovsky
            gg (int,float): Assymetry factor for Heyney-Greenstein case. Default is 0.
            prec (str): Required precision mode, accepted values "single" or "double".
                Default is double.

        Returns:
            pmom (array-like): Legendre coefficients of the phase function (moments).

        """
        # check input and raise warning if necessary
        if iphas == 3 and gg == 0:
            print(
                (
                    "Assymetry factor for the Heyney-Greenstein phase function"
                    " is 0 (default). If you wish to use a different value, pass"
                    " it to the function as gg = [your value]."
                )
            )

        # initialize empty pmom array of required shape
        pmom = np.zeros(
            shape=(self.disort_input["maxmom"] + 1, self.disort_input["maxcly"])
        )

        # fill the array with phase function moments one layer at a time
        for i in range(self.disort_input["maxcly"]):
            pmom[:, i] = dms.getmom(
                iphas=iphas, gg=gg, nmom=self.disort_input["maxmom"], pmom=pmom[:, i]
            )
        return pmom

    def calc_pmom_double(self, iphas, gg=0):
        """Calculates phase function moments from disort using the getmom function.

        This is a class method because it contains a loop which is inconvenient
        in the main code (getmom calcualtes phase function moments
        for a 1D array/list, not a 2D array.
        This function uses the double precision version of DISORT.

        Args:
            iphas (int): phase function option. Can be:
                - 1: Isotropic
                - 2: Rayleigh
                - 3: Henyey-Greenstein with asymmetry factor GG
                - 4: Haze L as specified by Garcia/Siewert
                - 5: Cloud C.1 as specified by Garcia/Siewert
                - 6: Aerosol as specified by Kokhanovsky
                - 7: Cloud as specified by Kokhanovsky
            gg (int,float): Assymetry factor for Heyney-Greenstein case. Default is 0.
            prec (str): Required precision mode, accepted values "single" or "double".
                Default is double.

        Returns:
            pmom (array-like): Legendre coefficients of the phase function (moments).

        """
        # check input and raise warning if necessary
        if iphas == 3 and gg == 0:
            print(
                (
                    "Assymetry factor for the Heyney-Greenstein phase function"
                    " is 0 (default). If you wish to use a different value, pass"
                    " it to the function as gg = [your value]."
                )
            )

        # initialize empty pmom array of required shape
        pmom = np.zeros(
            shape=(self.disort_input["maxmom"] + 1, self.disort_input["maxcly"])
        )

        # fill the array with phase function moments one layer at a time
        for i in range(self.disort_input["maxcly"]):
            pmom[:, i] = dmd.getmom(
                iphas=iphas, gg=gg, nmom=self.disort_input["maxmom"], pmom=pmom[:, i]
            )
        return pmom


class SRFM(Fwd_model):
    """This class represents the final forward model.

    The class object serves as a container for the outputs from various forward models,
    be it RFM + DISORT or other.

    """

    def __init__(self, name="SRFM", **parameters):
        super().__init__(name)
        for key, val in parameters.items():
            setattr(self, key, val)

    def initialize_srfm_output_arrays_from_disort(self, DISORT, retain_outputs=None):
        """Initialize requested spectral arrays for DISORT outputs.

        By default every historical output is allocated.  Passing a collection lets
        memory-sensitive callers retain only the values that they will return, plot,
        or write. ``"radiance"`` is accepted as an alias for ``"uu"``.

        Args:
            DISORT (obj): instance of srfm.forward_model.DISORT
            retain_outputs (collection[str] | None): Output names to allocate, or
                ``None`` to preserve the historical all-output behavior.

        Returns:
            None: Requested arrays are allocated as SRFM attributes.

        Raises:
            RuntimeError: Raised when the SRFM object doesn't have wavenumber or
                wavelengths grid first.
            ValueError: Raised when an output name is not recognized.

        """
        if hasattr(self, "wvnm") and self.wvnm is not None:
            dim = len(self.wvnm)
        elif hasattr(self, "wvls") and self.wvls is not None:
            dim = len(self.wvls)
        else:
            raise RuntimeError("SRFM must have wvnm or wvls grids first.")

        output_names = {
            "rfldir",
            "rfldn",
            "flup",
            "dfdt",
            "uavg",
            "uu",
            "albmed",
            "trnmed",
        }
        if retain_outputs is None:
            retained = output_names
        else:
            aliases = {"radiance": "uu"}
            requested = {aliases.get(name, name) for name in retain_outputs}
            requested.discard("bbt")
            unknown = requested - output_names
            if unknown:
                raise ValueError(
                    "Unknown SRFM output name(s): " + ", ".join(sorted(unknown))
                )
            retained = requested

        level_shape = (dim, DISORT.disort_input["maxulv"])
        angle_shape = (dim, DISORT.disort_input["maxumu"])
        shapes = {
            "rfldir": level_shape,
            "rfldn": level_shape,
            "flup": level_shape,
            "dfdt": level_shape,
            "uavg": level_shape,
            "uu": (
                dim,
                DISORT.disort_input["maxumu"],
                DISORT.disort_input["maxulv"],
                DISORT.disort_input["maxphi"],
            ),
            "albmed": angle_shape,
            "trnmed": angle_shape,
        }
        self.retained_outputs = frozenset(retained)
        for output_name in retained:
            setattr(self, output_name, np.zeros(shapes[output_name]))
        return

    def set_wvnm(self, wvnm):
        """Assigns wavenumber grid to SRFM.

        Args:
            wvnm (array-like): Wavenumber grid.

        """
        if hasattr(self, "wvls") and self.wvls is not None:
            assert len(self.wvls) == len(
                wvnm
            ), """wvls and wvnm do not have equal 
            number of points."""
        self.wvnm = wvnm
        return

    def set_wvls(self, wvls):
        """Assigns wavelength grid to SRFM.

        Args:
            wvls (array-like): Wavelength grid.

        """
        if hasattr(self, "wvnm") and self.wvnm is not None:
            assert len(wvls) == len(
                self.wvnm
            ), """wvls and wvnm do not have equal 
            number of points."""
        self.wvls = wvls
        return

    def store_disort_result(self, result, wvl_idx):
        """Stores results from a single DISORT run into the SRFM object.

        DISORT returns results for a given wavenumber/wavelength. If the overarching
        idea is to calculate a spectrum, which is it, DISORT is run in a lopp. This
        function takes results from a single DISORT run and inserts them in SRFM arrays
        in appropriate places (at appropriate indices).

        Args:
            result (DisortResult | DISORT): Structured current result. A DISORT
                instance is also accepted for compatibility and uses its
                ``current_output`` or retained legacy dictionary entry.
            wvl_idx (int): Values are inserted into SRFM arrays at this index. The idea
                is that the DISORT calculation is performed at a certain wavenunmber.
                SRFM has initialized arrays of size matching the overall wavenumber grid
                and results from each DISORT run are inserted into the arrays at the
                index corresponding to the respective wavenumber.

        Returns:
            None: Retained arrays are updated in place.

        """
        if isinstance(result, DISORT):
            if result.current_output is not None:
                result = result.current_output
            else:
                result = result.disort_out[result.wvnm]

        for output_name in getattr(self, "retained_outputs", ()):
            source = (
                getattr(result, output_name)
                if isinstance(result, DisortResult)
                else result[output_name]
            )
            getattr(self, output_name)[wvl_idx] = source
        return

    def calc_bbt(self):
        """Converts radiance to brightness temperature.

        Expects SRFM to have wavenumber/wavelength grid and radiances (uu) as
        attributes. Calculates brightness temperature array which matches the radiance
        array in shape.

        """
        wvnm = np.asarray(self.wvnm, dtype=float)
        if self.uu.ndim < 1 or wvnm.ndim != 1 or wvnm.size != self.uu.shape[0]:
            raise ValueError(
                "The wavenumber grid must be one-dimensional and match the first "
                "radiance dimension."
            )

        # Broadcast the spectral grid over every DISORT output dimension
        # without allocating repeated copies.
        wvnm = wvnm.reshape((wvnm.size,) + (1,) * (self.uu.ndim - 1))

        self.bbt = utils.convert_spectral_radiance_to_bbt(self.uu, wvnm)
        return

    def convolve_with_iasi(self, filename):
        """Convolve radiance (uu) with IASI instrument line shape.

        The SRFM object must contain radiances (uu) and a wavenumber grid.
        Assumes regular grid.
        Note that the convolved spectrum suffers from boundary effects (scipy
        interpolate does zero padding of the data at the boundaries). Best avoided by
        calculating your original spectra at a wider interval and then interpolating/
        truncating. For interpolation (best used to get spectra at your simulated
        satellite grid) see the interp() function of this module.

        Args:
            filename (str): filename of the iasi.ils file (IASI instrument line shape
                kindly provided by Anu Dudhia, in RFM format.)

        """
        uu_unconvolved = np.asarray(self.uu)
        if uu_unconvolved.ndim < 1 or uu_unconvolved.shape[0] != len(self.wvnm):
            raise ValueError(
                "The first radiance dimension must match the wavenumber grid."
            )

        # read instrument line shape
        ils_x, ils_y, ils_lo, ils_hi = utils.read_ils(filename)

        # check if model grid is regular
        spacing = np.diff(np.asarray(self.wvnm, dtype=float))
        if spacing.size == 0 or not np.allclose(
            spacing, spacing[0], rtol=1e-10, atol=1e-12
        ):
            raise ValueError("Wavenumber grid is not regular.")

        # determine resolution from model wavenumber grid
        num = len(self.wvnm)
        lo = self.wvnm.min()
        hi = self.wvnm.max()
        res = np.round((hi - lo) / (num - 1), decimals=8)
        # this inadvertently introduces a limit
        # on the minimum resolution used in the code as 1e-8 cm-1, which should be
        # enough though, and also this may not be the numerically most stable way to go

        # generate new grid for ils
        npts = (
            int(np.floor((ils_hi - ils_lo) / res)) + 1
        )  # expected number of points in the grid
        new_x = ils_lo + np.arange(npts) * res
        #        new_x = np.linspace(ils_lo, ils_hi, int((ils_hi - ils_lo) / res + 1))

        # interpolate ils to new grid
        new_y = np.interp(new_x, ils_x, ils_y)

        # calculate sum of instrument line shape for normalization later
        norm = np.sum(new_y)
        if not np.isfinite(norm) or np.isclose(norm, 0.0):
            raise ValueError("Instrument line shape has zero or non-finite normalization.")

        # determine shape of uu from DISORT (a set of output spectra at
        # different optical depths, polar angles, and azimuthal angles)
        uu_shape = self.uu.shape  # tuple
        nwv = uu_shape[0]  # first dimension size
        rest = int(np.prod(uu_shape[1:], dtype=int))
        # rest basically gives a number of stored spectra in the variable

        # reshape uu (view)
        uu_flat = uu_unconvolved.reshape(nwv, rest)
        out_flat = np.empty_like(uu_flat)

        # loop over columns (each column is a spectrum)
        for j in range(rest):
            out_flat[:, j] = convolve(uu_flat[:, j], new_y, mode="same") / norm

        # reshape back
        self.uu = out_flat.reshape(uu_shape)

        if hasattr(self, "bbt"):
            self.calc_bbt()

        return

    def interp(self, new_wvnm):
        """Interpolate every retained spectral output to a new grid.

        Original intended use is to interpolate the calculated and already convolved
        spectra (i.e. at a lower effective resolution) to a satellite lower resolution
        grid.

        All retained raw DISORT fields are interpolated consistently. If brightness
        temperature already exists, it is recalculated from interpolated radiance
        because that conversion is nonlinear.

        Args:
            new_wvnm (array-like): New wavenumber grid in cm-1.

        Returns:
            None: Spectral arrays and coordinate grids are updated in place.

        Raises:
            ValueError: If either grid is invalid, extrapolation would be required,
                or an output's spectral dimension is inconsistent.

        """

        old_wvnm = np.asarray(self.wvnm, dtype=float)
        new_wvnm = np.asarray(new_wvnm, dtype=float)
        if old_wvnm.ndim != 1 or new_wvnm.ndim != 1:
            raise ValueError("Wavenumber grids must be one-dimensional.")
        if old_wvnm.size < 2 or not np.all(np.diff(old_wvnm) > 0):
            raise ValueError("The source wavenumber grid must be strictly increasing.")
        if not np.all(np.isfinite(new_wvnm)):
            raise ValueError("The new wavenumber grid must contain only finite values.")
        if new_wvnm.size and (
            new_wvnm.min() < old_wvnm[0] or new_wvnm.max() > old_wvnm[-1]
        ):
            raise ValueError("The new wavenumber grid must lie within the source grid.")
        spectral_outputs = set(getattr(self, "retained_outputs", ()))
        if not spectral_outputs:
            spectral_outputs = {
                name
                for name in (
                    "rfldir",
                    "rfldn",
                    "flup",
                    "dfdt",
                    "uavg",
                    "uu",
                    "albmed",
                    "trnmed",
                )
                if hasattr(self, name)
            }
        for output_name in spectral_outputs:
            old_values = np.asarray(getattr(self, output_name))
            if old_values.shape[0] != old_wvnm.size:
                raise ValueError(
                    f"The first {output_name} dimension must match the source "
                    "wavenumber grid."
                )
            interpolator = make_interp_spline(
                old_wvnm, old_values, axis=0, k=1
            )
            setattr(self, output_name, interpolator(new_wvnm))

        # calculate new wavelengths [um]
        new_wvls = (1 / new_wvnm) * 1e4

        # assign grids to class
        self.wvnm = new_wvnm
        self.wvls = new_wvls

        # Brightness temperature is nonlinear in radiance, so recalculate an
        # existing attribute instead of interpolating it independently.
        if hasattr(self, "bbt"):
            self.calc_bbt()

        return
