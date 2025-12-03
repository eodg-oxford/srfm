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
            self.key = val

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
        """Determines total column density of a species in the atmosphere.

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
    """Class that contains the DISORT forward model.

    Is subclass of Fwd_model.
    """

    def __init__(
        self,
        name="DISORT",
        disort_fldr=None,
        disort_input={},
        disort_out={},
        disort_fmt_passmark=True,
        disort_integrity_passmark=True,
        status="DISORT object created.",
        **parameters,
    ):
        super().__init__(name)
        self.disort_fldr = disort_fldr
        self.disort_input = disort_input
        self.disort_out = disort_out
        self.disort_fmt_passmark = disort_fmt_passmark
        self.disort_integrity_passmark = disort_integrity_passmark
        self.status = status
        for key, val in parameters.items():
            self.key = val

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

        Raises:
            AttributeError: Raised when rfm does not have rfm_output attribute.
        """
        try:
            self.disort_input["maxcly"] = len(rfm.rfm_output["layer no."])
        except AttributeError:
            print(
                "Class RFM does not have attribute rfm_output. Make sure to set up class RFM properly first."
            )

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

        """
        try:
            self.disort_input["temper"] = rfm.rfm_output["T_upper (K)"].tolist()
            self.disort_input["temper"].append(
                rfm.rfm_output["T_lower (K)"].tolist()[-1]
            )
        except AttributeError:
            print(
                """Class RFM does not have attribute rfm_output. Make sure to set-up 
                class RFM properly first."""
            )
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

        Raises:
            AttributeError: Raised when RFM object does not have rfm_output attribute.
        """

        try:
            self.disort_input["h_lyr"] = rfm.rfm_output["h_upper (km)"].tolist()
            self.disort_input["h_lyr"].append(
                rfm.rfm_output["h_lower (km)"].tolist()[-1]
            )
        except AttributeError:
            print(
                "Class RFM does not have attribute rfm_output. Make sure to set-up class RFM properly first."
            )

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

    def run_disort(self, prec="double"):
        """Calls function to run disort with required precision."""
        if prec == "double":
            self.run_disort_double()
        elif prec == "single":
            self.run_disort_single()
        else:
            raise ValueError("prec must be 'single' or 'double'.")
        return

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

    def run_disort_single(self):
        """Runs disort, single precision."""
        if self.disort_fmt_passmark == True:
            if self.disort_integrity_passmark == True:
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
                self.status = "DISORT run completed."
            else:
                raise ValueError(
                    "Disort input passmark is not True, run input integrity test first."
                )
        else:
            raise ValueError(
                "Disort format passmark is not True, run format test first."
            )

        return res

    def run_disort_double(self):
        """Runs disort, double precision."""
        if self.disort_fmt_passmark == True:
            if self.disort_integrity_passmark == True:
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

                self.disort_out[self.wvnm] = {}
                self.disort_out[self.wvnm]["wavenumber (cm-1)"] = self.wvnm
                self.disort_out[self.wvnm]["wavelength (um)"] = self.wvl
                self.disort_out[self.wvnm]["rfldir"] = res[0] / (
                    self.disort_input["wvnmhi"] - self.disort_input["wvnmlo"]
                )
                self.disort_out[self.wvnm]["rfldn"] = res[1] / (
                    self.disort_input["wvnmhi"] - self.disort_input["wvnmlo"]
                )
                self.disort_out[self.wvnm]["flup"] = res[2] / (
                    self.disort_input["wvnmhi"] - self.disort_input["wvnmlo"]
                )
                self.disort_out[self.wvnm]["dfdt"] = res[3] / (
                    self.disort_input["wvnmhi"] - self.disort_input["wvnmlo"]
                )
                self.disort_out[self.wvnm]["uavg"] = res[4] / (
                    self.disort_input["wvnmhi"] - self.disort_input["wvnmlo"]
                )
                self.disort_out[self.wvnm]["uu"] = res[5] / (
                    self.disort_input["wvnmhi"] - self.disort_input["wvnmlo"]
                )
                self.disort_out[self.wvnm]["albmed"] = res[6]
                if len(res) == 8:
                    self.disort_out[self.wvnm]["trnmed"] = res[7]
                else:
                    self.disort_out[self.wvnm]["trnmed"] = 0

                self.status = "DISORT run completed."
            else:
                raise ValueError(
                    "Disort input passmark is not True, run input integrity test first."
                )
        else:
            raise ValueError(
                "Disort format passmark is not True, run format test first."
            )

        return

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
        if tau_R.shape != (maxcly,) or tau_p.shape != (maxcly,) or w_p.shape != (maxcly,):
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
            self.key = val

    def initialize_srfm_output_arrays_from_disort(self, DISORT):
        """Initializes srfm output arrays to which disort outputs are appended.

        Args:
            DISORT (obj): instance of srfm.forward_model.DISORT

        Raises:
            RuntimeError: Raised when the SRFM object doesn't have wavenumber or
                wavelengths grid first.

        """
        if hasattr(self, "wvnm") and self.wvnm is not None:
            dim = len(self.wvnm)
        elif hasattr(self, "wvls") and self.wvls is not None:
            dim = len(self.wvls)
        else:
            raise RuntimeError("SRFM must have wvnm or wvls grids first.")

        self.rfldir = np.zeros((dim, DISORT.disort_input["maxulv"]))
        self.rfldn = np.zeros((dim, DISORT.disort_input["maxulv"]))
        self.flup = np.zeros((dim, DISORT.disort_input["maxulv"]))
        self.dfdt = np.zeros((dim, DISORT.disort_input["maxulv"]))
        self.uavg = np.zeros((dim, DISORT.disort_input["maxulv"]))

        self.uu = np.zeros(
            (
                dim,
                DISORT.disort_input["maxumu"],
                DISORT.disort_input["maxulv"],
                DISORT.disort_input["maxphi"],
            )
        )
        self.albmed = np.zeros((dim, DISORT.disort_input["maxumu"]))
        self.trnmed = np.zeros((dim, DISORT.disort_input["maxumu"]))
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

    def store_disort_result(self, DISORT, wvl_idx):
        """Stores results from a single DISORT run into the SRFM object.

        DISORT returns results for a given wavenumber/wavelength. If the overarching
        idea is to calculate a spectrum, which is it, DISORT is run in a lopp. This
        function takes results from a single DISORT run and inserts them in SRFM arrays
        in appropriate places (at appropriate indices).

        Args:
            DISORT (obj): instance of srfm.forward_model.DISORT
            wvl_idx (int): Values are inserted into SRFM arrays at this index. The idea
                is that the DISORT calculation is performed at a certain wavenunmber.
                SRFM has initialized arrays of size matching the overall wavenumber grid
                and results from each DISORT run are inserted into the arrays at the
                index corresponding to the respective wavenumber.

        """
        wvnm = DISORT.wvnm
        self.rfldir[wvl_idx] = DISORT.disort_out[wvnm]["rfldir"]
        self.rfldn[wvl_idx] = DISORT.disort_out[wvnm]["rfldn"]
        self.flup[wvl_idx] = DISORT.disort_out[wvnm]["flup"]
        self.dfdt[wvl_idx] = DISORT.disort_out[wvnm]["dfdt"]
        self.uavg[wvl_idx] = DISORT.disort_out[wvnm]["uavg"]
        self.uu[wvl_idx] = DISORT.disort_out[wvnm]["uu"]
        self.albmed[wvl_idx] = DISORT.disort_out[wvnm]["albmed"]
        self.trnmed[wvl_idx] = DISORT.disort_out[wvnm]["trnmed"]
        return

    def calc_bbt(self):
        """Converts radiance to brightness temperature.

        Expects SRFM to have wavenumber/wavelength grid and radiances (uu) as
        attributes. Calculates brightness temperature array which matches the radiance
        array in shape.

        """
        # strech wvnm to correct shape to be broadcastable.
        wvnm = self.wvnm[:, np.newaxis, np.newaxis, np.newaxis]  # add new axis to wvnm
        # to match the number of uu dimensions, 0th dimension (axis 0) are the same

        assert (
            wvnm.shape[0] == self.uu.shape[0]
        ), """wvnm and uu don't have the same
        shape of the first axis???"""

        uu_shape = self.uu.shape  # tuple

        for num, i in enumerate(uu_shape[1:]):
            wvnm = np.repeat(wvnm, i, axis=num + 1)

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
        # save a copy of unconvolved spectrum
        uu_unconvolved = self.uu.copy()

        # read instrument line shape
        ils_x, ils_y, ils_lo, ils_hi = utils.read_ils(filename)

        # check if model grid is regular
        a = np.diff(self.wvnm, n=2)  # calculate 2nd discrete difference
        a[a < 1e12] = 0  # remove small numbers (arising from computer precision limits)
        assert (
            np.all(a) == False
        ), """Wavenumber grid is 
        not regular."""  # check is all values in a are 0 (0 evaluates to False)

        # determine resolution from model wavenumber grid
        num = len(self.wvnm)
        lo = self.wvnm.min()
        hi = self.wvnm.max()
        res = np.round(
            (hi - lo) / num, decimals=8
        )  # this is inadvertedly introduces a limit
        # on the minimum resolution used in the code as 1e-8 cm-1, which should be
        # enough though, and also this may not be the numerically most stable way to go

        # generate new grid for ils
        new_x = np.linspace(ils_lo, ils_hi, int((ils_hi - ils_lo) / res + 1))

        # interpolate ils to new grid
        new_y = np.interp(new_x, ils_x, ils_y)

        # calculate sum of instrument line shape for normalization later
        norm = np.sum(new_y)

        # determine shape of uu from DISORT (basically a set of output spectra a
        # different optical dpeths, polar and azimuthal angles
        uu_shape = self.uu.shape  # tuple

        # determine all combinations of indices of uu
        combs = []
        for i in range(uu_shape[1]):
            for ii in range(uu_shape[2]):
                for iii in range(uu_shape[3]):
                    combs.append([i, ii, iii])

        # convolve spectra in a loop
        for c in combs:
            self.uu[:, c[0], c[1], c[2]] = (
                convolve(uu_unconvolved[:, c[0], c[1], c[2]], new_y, mode="same") / norm
            )

        return

    def interp(self, new_wvnm):
        """Interpolates radiances (uu) and brightness temperatures to a new grid.

        Original intended use is to interpolate the calculated and already convolved
        spectra (i.e. at a lower effective resolution) to a satellite lower resolution
        grid.

        The SRFM object must contain radiances (uu). If besides radiances also contains
        brightness temperatures (bbt), then these are interpolated as well.

        Output is returned as updated attributes - new wavenumber and wavelength grids
        as well as new radiances and brightness temperatures (if originally present).

        Args:
            new_wvnm: new wavenumber grid to interpolate to, units [cm-1]

        """

        uu_shape = self.uu.shape  # tuple

        new_uu_shape = list(uu_shape)
        new_uu_shape[0] = len(new_wvnm)
        new_uu = np.zeros(tuple(new_uu_shape))

        # determine all combinations of indices of uu
        combs = []
        for i in range(uu_shape[1]):
            for ii in range(uu_shape[2]):
                for iii in range(uu_shape[3]):
                    combs.append([i, ii, iii])

        # interpolate spectra in a loop
        for c in combs:
            new_uu[:, c[0], c[1], c[2]] = np.interp(
                new_wvnm, self.wvnm, self.uu[:, c[0], c[1], c[2]]
            )

        self.uu = new_uu

        # calculate new wavelengths [um]
        new_wvls = (1 / new_wvnm) * 1e4

        # assign grids to class
        self.wvnm = new_wvnm
        self.wvls = new_wvls

        return
