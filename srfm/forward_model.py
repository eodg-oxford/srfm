"""
Name: forward_model
Parent package: srfm
Author: Antonin Knizek
Contributors: 
Date: 18 February 2025
Purpose: Defines the base model class for the forward model and specific model 
subclasses, currently RFM and DISORT, along with their methods.
""" 
import numpy as np
from . import disort_functions as disf
from . import rfm_functions as rf
import os
from . import utilities as utils
from . import units
import pandas as pd
from scipy.signal import convolve

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
    def __init__(self, name=None, **parameters):
        self.name = name
        self.parameters = {}


class RFM(Fwd_model):
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
        """Attempts to run rfm from python. Can also be run manually.
        inputs are;
        fldr: the RFM folder (relative or absolute path)
        wipe: if True, removes rfm results from the previous run, default=True"""
        
        cwd = os.getcwd()
        try:
            os.chdir(f"{fldr}")
        except OSError:
            print("Invalid directory.")
            return
        
        if wipe==True:
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
        self.rfm_output = rf.get_rfm_optical_depths(fldr=fldr, levels=levels)
        self.status = "RFM run completed and result loaded."
        return

    def get_wnos_from_RFM(self):
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
    
    def load_output_prf(self,fldr):
        """"Loads the output profile (default prf.asc) file from RFM.
        fldr - path to RFM"""
        self.output_prf = rf.read_output_prf(f"{fldr}/prf.asc")
        return

class DISORT(Fwd_model):
    def __init__(
        self,
        name="DISORT",
        disort_fldr=None,
        disort_input={},
        disort_out={},
        disort_fmt_passmark=False,
        disort_integ_passmark=False,
        status="DISORT object created.",
        **parameters,
    ):
        super().__init__(name)
        self.disort_fldr = disort_fldr
        self.disort_input = disort_input
        self.disort_out = disort_out
        self.disort_fmt_passmark = disort_fmt_passmark
        self.disort_integ_passmark = disort_integ_passmark
        self.status = status
        for key, val in parameters.items():
            self.key = val

    def add_disort_input(self, d):
        """add disort input parameters as a dictionary."""
        self.disort_in = d

    def add_disort_empty_input(self):
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
        """test disort input for correct formatting."""
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
        self.disort_input["maxcly"] = maxcly

    def set_maxcly_from_rfm(self, rfm):
        """sets maxcly from rfm
        rfm is a class RFM object and has attribute rfm_output"""
        try:
            self.disort_input["maxcly"] = len(rfm.rfm_output["layer no."])
        except AttributeError:
            print(
                "Class RFM does not have attribute rfm_output. Make sure to set up class RFM properly first."
            )

    def set_maxmom(self, maxmom):
        self.disort_input["maxmom"] = maxmom

    def set_maxcmu(self, maxcmu):
        self.disort_input["maxcmu"] = maxcmu

    def set_maxumu(self, maxumu):
        self.disort_input["maxumu"] = maxumu

    def set_maxphi(self, maxphi):
        self.disort_input["maxphi"] = maxphi

    def set_maxulv(self, maxulv):
        self.disort_input["maxulv"] = maxulv

    def set_usrang(self, usrang):
        self.disort_input["usrang"] = usrang

    def set_usrtau(self, usrtau):
        self.disort_input["usrtau"] = usrtau

    def set_ibcnd(self, ibcnd):
        self.disort_input["ibcnd"] = ibcnd

    def set_onlyfl(self, onlyfl):
        self.disort_input["onlyfl"] = onlyfl

    def set_prnt(self, prnt):
        self.disort_input["prnt"] = prnt

    def set_plank(self, planck):
        self.disort_input["plank"] = planck

    def set_lamber(self, lamber):
        self.disort_input["lamber"] = lamber

    def set_deltamplus(self, deltamplus):
        self.disort_input["deltamplus"] = deltamplus

    def set_do_pseudo_sphere(self, d_p_s):
        self.disort_input["do_pseudo_sphere"] = d_p_s

    def set_dtauc_manually(self, dtauc):
        self.disort_input["dtauc"] = dtauc

    def set_ssalb_manually(self, ssalb):
        self.disort_input["ssalb"] = ssalb

    def set_pmom_manually(self, pmom):
        self.disort_input["pmom"] = pmom

    def set_temper(self, temper):
        self.disort_input["temper"] = temper

    def set_temper_from_rfm(self, rfm):
        """Sets layer temperature for disort from the RFM output.
        Takes upper temperatures from rfm_output for every layer.
        Adds the lowest layer lower temperature (adjacent to surface) as the last one.
        rfm is a class RFM object and has attribute rfm_output
        """
        try:
            self.disort_input["temper"] = rfm.rfm_output["T_upper (K)"].tolist()
            self.disort_input["temper"].append(
                rfm.rfm_output["T_lower (K)"].tolist()[-1]
            )
        except AttributeError:
            print(
                "Class RFM does not have attribute rfm_output. Make sure to set-up class RFM properly first."
            )
        return

    def set_wvnm_range(self, lo, hi):
        self.disort_input["wvnmlo"] = lo
        self.disort_input["wvnmhi"] = hi

    def set_utau(self, utau):
        self.disort_input["utau"] = utau

    def set_umu0(self, umu0):
        self.disort_input["umu0"] = umu0

    def set_phi0(self, phi0):
        self.disort_input["phi0"] = phi0

    def set_umu(self, umu):
        self.disort_input["umu"] = umu

    def set_phi(self, phi):
        self.disort_input["phi"] = phi

    def set_fbeam(self, fbeam):
        self.disort_input["fbeam"] = fbeam

    def set_fisot(self, fisot):
        self.disort_input["fisot"] = fisot

    def set_albedo(self, albedo):
        self.disort_input["albedo"] = albedo

    def set_btemp(self, btemp):
        self.disort_input["btemp"] = btemp

    def set_ttemp(self, ttemp):
        self.disort_input["ttemp"] = ttemp

    def set_temis(self, temis):
        self.disort_input["temis"] = temis

    def set_earth_radius(self, earth_radius):
        self.disort_input["earth_radius"] = earth_radius

    def set_h_lyr(self, h_lyr):
        self.disort_input["h_lyr"] = h_lyr

    def set_h_lyr_from_rfm(self):
        """Sets h_lyr for disort from the RFM output.
        rfm is a class RFM object and has attribute rfm_output
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
        self.disort_input["rhoq"] = rhoq

    def set_rhou(self, rhou):
        self.disort_input["rhou"] = rhou

    def set_rho_accurate(self, rho_accurate):
        self.disort_input["rho_accurate"] = rho_accurate

    def set_bemst(self, bemst):
        self.disort_input["bemst"] = bemst

    def set_emust(self, emust):
        self.disort_input["emust"] = emust

    def set_accur(self, accur):
        self.disort_input["accur"] = accur

    def set_header(self, header):
        self.disort_input["header"] = header

    def initialize_disort_output_arrays(self):
        """Initializes output arrays for a single disort run."""
        self.disort_input["rfldir"] = np.zeros(self.disort_input["maxulv"])
        self.disort_input["rfldn"] = np.zeros(self.disort_input["maxulv"])
        self.disort_input["flup"] = np.zeros(self.disort_input["maxulv"])
        self.disort_input["dfdt"] = np.zeros(self.disort_input["maxulv"])
        self.disort_input["uavg"] = np.zeros(self.disort_input["maxulv"])
        self.disort_input["uu"] = np.zeros(
                                        (
                                    self.disort_input["maxumu"],
                                    self.disort_input["maxulv"],
                                    self.disort_input["maxphi"]
                                )
                            )
        self.disort_input["albmed"] = np.zeros(self.disort_input["maxumu"])
        self.disort_input["trnmed"] = np.zeros(self.disort_input["maxumu"])
    
    def run_disort(self,prec="double"):
        """Calls function to run disort with required precision."""
        if prec == "double":
            self.run_disort_double()
        elif prec == "single":
            self.run_disort_single()
        else:
            raise ValueError("prec must be 'single' or 'double'.")
        return
    
    def set_wvnm(self,wvnm):
        """Sets wavenumber of the current run."""
        self.wvnm = wvnm
    
    def set_wvl(self,wvl):
        """Set wavelength of the current run."""
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

    def set_dtauc(self, tau_g,tau_R,tau_p):
        """Calculate and add optical depth of model layers, delta tau (dtau).
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
        
        self.disort_input["dtauc"] = tau_g + tau_R + tau_p
        return 
        
        
    def set_ssalb(self, tau_g,tau_R,tau_p,w_p):
        """Calculates and adds the single scatter albedo of the model layers
        according to:
        w = (w_g*tau_g + w_R*tau_R + w_p*tau_p) / (tau_g + tau_R + tau_p)
        where:
        w - layer single scatter albedo
        w_g - layer single scatter albedo from gas absorption, == 0
        w_R - layer Rayleigh scattering single scatter albedo, == 1
        w_p - layer particle scattering single scatter albedo
        tau_g - layer gas absorption optical depth
        tau_R - layer Rayleigh scattering optical depth
        tau_p - layer particle scattering optical depth
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
        w_p = check_convert_dtype(w_p)        
        
        ssalb = np.nan_to_num(((tau_R + w_p*tau_p) / (tau_R + tau_g + tau_p)))

        self.disort_input["ssalb"] = ssalb
        return
        
                
        
    def set_pmom(self, pmom_R, tau_R, w_p, tau_p, pmom_p):
        """Calculate phase function coefficients according to Don (from ORAC).
        Formula:
            xi = (w_g*tau_g*xi,g + w_R*tau_R*xi,r + w_p*tau_p*xi,p) 
                 / (w_g*tau_g + w_r*tau_R + w_p*tau_p)
            where xi - Legendre polynomial coefficient 
                  w_g = 0 - gas single scatter albedo
                  xi,g = 0, gas phase function moment
                  w_R = 1, Rayleigh scattering single scatter albedo
            which makes the formula:
            xi = (tau_R*xi,r + w_p*tau_p*xi,p) 
                 / (tau_R + w_p*tau_p)
        The function works on arrays, so xi is replaced by array pmom:
         
        pmom_R - array of phase function moments for Rayleigh scattering
               - 2D array, where 'columns' are the atmospheric layers
                                 'rows' are the Legendre polynomial 
                                        coefficients x1,R to xn,R
               - has shape (model_DISORT.disort_input["maxmom"] + 1,
                            model_DISORT.disort_input["maxcly"])
        pmom_p - array of phase function coefficients for particle scattering
               - same shape and meaning as pmom_R       
        w_p - particle scattering single scatter albedo for each layer, list/1D array
        tau_p - particle scattering optical depth for each layer, list/1D array
        tau_R - Rayleigh optical depths
        """
        def test_array(obj):
            if isinstance(obj, np.ndarray):
                if obj.shape == (self.disort_input["maxmom"] + 1,
                            self.disort_input["maxcly"]):
                    return
                else:
                    raise ValueError('Inputs must have shape (model_DISORT.disort_input["maxmom"] + 1, model_DISORT.disort_input["maxcly"]).')
            else:
                raise TypeError(f'Inputs must be np.ndarrays.')
        
        tau_R = np.stack((tau_R,)*(self.disort_input["maxmom"]+1),axis=0)
        tau_p = np.stack((tau_p,)*(self.disort_input["maxmom"]+1),axis=0)
        w_p = np.stack((w_p,)*(self.disort_input["maxmom"]+1),axis=0)
        
        for i in [pmom_R, tau_R, w_p, tau_p, pmom_p]:
            test_array(i)
        
        pmom = (tau_R*pmom_R + w_p*tau_p*pmom_p) / (tau_R + w_p*tau_p)
        self.disort_input["pmom"] = np.nan_to_num(pmom)
        return
        
    def calc_pmom(self, iphas, prec="double", gg=0):
        """Calls function to calculate phase function moments from disort using 
        the getmom fucntion.
        iphas - phase function option
                1 : Isotropic
                2 : Rayleigh
                3 : Henyey-Greenstein with asymmetry factor GG
                4 : Haze L as specified by Garcia/Siewert
                5 : Cloud C.1 as specified by Garcia/Siewert
                6 : Aerosol as specified by Kokhanovsky 
                7 : Cloud as specified by Kokhanovsky
        gg - assymetry factor for Heyney-Greenstein case
        prec - required precision mode, accepted values "single" or "double" (default)
        """
        
        if prec == "single":
            pmom = self.calc_pmom_single(iphas,gg)
        elif prec == "double":
            pmom = self.calc_pmom_double(iphas,gg)
        else:
            raise ValueError("prec must be 'single' or 'double'.")
        return pmom
    
    def calc_pmom_single(self, iphas, gg=0):
        """Calculates phase function moments from disort using the getmom
        function.
        This is a class method because it contains a loop which is inconvenient
        in the main code (getmom calcualtes phase function moments 
        for a 1D array/list, not a 2D array.
        pmom - empty pmom array of disort-required shape
        iphas - phase function option
                1 : Isotropic
                2 : Rayleigh
                3 : Henyey-Greenstein with asymmetry factor GG
                4 : Haze L as specified by Garcia/Siewert
                5 : Cloud C.1 as specified by Garcia/Siewert
                6 : Aerosol as specified by Kokhanovsky 
                7 : Cloud as specified by Kokhanovsky
        gg - assymetry factor for Heyney-Greenstein case
        nmom - index of the highest Legendre coefficient needed, 
             - set automatically
        single precision
        """
        # check input and raise warning if necessary
        if iphas == 3 and gg == 0:
            print(("Assymetry factor for the Heyney-Greenstein phase function"
                   " is 0 (default). If you wish to use a different value, pass"
                   " it to the function as gg = [your value]."))
        
        # initialize empty pmom array of required shape
        pmom = np.zeros(
            shape=(self.disort_input["maxmom"] + 1,
            self.disort_input["maxcly"])
        )
        
        # fill the array with phase function moments one layer at a time
        for i in range(self.disort_input["maxcly"]):
            pmom[:,i] = dms.getmom(iphas=iphas,
                                  gg=gg,
                                  nmom=self.disort_input["maxmom"],
                                  pmom=pmom[:,i]
                              )
        return pmom

    def calc_pmom_double(self, iphas, gg=0):
        """Calculates phase function moments from disort using the getmom
        function.
        This is a class method because it contains a loop which is inconvenient
        in the main code (getmom calcualtes phase function moments 
        for a 1D array/list, not a 2D array.
        pmom - empty pmom array of disort-required shape
        iphas - phase function option
                1 : Isotropic
                2 : Rayleigh
                3 : Henyey-Greenstein with asymmetry factor GG
                4 : Haze L as specified by Garcia/Siewert
                5 : Cloud C.1 as specified by Garcia/Siewert
                6 : Aerosol as specified by Kokhanovsky 
                7 : Cloud as specified by Kokhanovsky
        gg - assymetry factor for Heyney-Greenstein case
        nmom - index of the highest Legendre coefficient needed, 
             - set automatically
        double precision
        """
        # check input and raise warning if necessary
        if iphas == 3 and gg == 0:
            print(("Assymetry factor for the Heyney-Greenstein phase function"
                   " is 0 (default). If you wish to use a different value, pass"
                   " it to the function as gg = [your value]."))
        
        # initialize empty pmom array of required shape
        pmom = np.zeros(
            shape=(self.disort_input["maxmom"] + 1,
            self.disort_input["maxcly"])
        )
        
        # fill the array with phase function moments one layer at a time
        for i in range(self.disort_input["maxcly"]):
            pmom[:,i] = dmd.getmom(iphas=iphas,
                                  gg=gg,
                                  nmom=self.disort_input["maxmom"],
                                  pmom=pmom[:,i]
                              )
        return pmom

class SRFM(Fwd_model):
    """This class represents the final forward model.
    The class object serves as a container for the outputs from various forward models,
    be it RFM + DISORT or other.
    """
    def __init__(self,name="SRFM",**parameters):
        super().__init__(name)        
        for key, val in parameters.items():
            self.key = val
            
    def initialize_srfm_output_arrays_from_disort(self,DISORT):
        """Initializes srfm output arrays to which disort outputs are appended.
        Input:
            self
            DISORT - srfm.forward_model.DISORT object
        Outputs:
            self.rfldir - float64 array, shape(wvnm/wvls,maxulv)
        """  
        if hasattr(self, "wvnm") and self.wvnm is not None:
            dim = len(self.wvnm)
        elif hasattr(self, "wvls") and self.wvls is not None:
            dim = len(self.wvls)
        else:
            raise RuntimeError("SRFM must have wvnm or wvls grids first.")
        
        
        self.rfldir = np.zeros((dim,DISORT.disort_input["maxulv"]))
        self.rfldn = np.zeros((dim,DISORT.disort_input["maxulv"]))
        self.flup = np.zeros((dim,DISORT.disort_input["maxulv"]))
        self.dfdt = np.zeros((dim,DISORT.disort_input["maxulv"]))
        self.uavg = np.zeros((dim,DISORT.disort_input["maxulv"]))
        
        self.uu = np.zeros((dim,DISORT.disort_input["maxumu"],DISORT.disort_input["maxulv"],DISORT.disort_input["maxphi"]))
        self.albmed = np.zeros((dim,DISORT.disort_input["maxumu"]))
        self.trnmed = np.zeros((dim,DISORT.disort_input["maxumu"]))
        return
        
    def set_wvnm(self, wvnm):
        """Assigns wavenumber grid."""
        if hasattr(self,"wvls") and self.wvls is not None:
            assert len(self.wvls) == len(wvnm), """wvls and wvnm do not have equal 
            number of points."""
        self.wvnm = wvnm
        return
    
    def set_wvls(self,wvls):
        """Assigns wavelength grid."""
        if hasattr(self, "wvnm") and self.wvnm is not None:
            assert len(wvls) == len(self.wvnm), """wvls and wvnm do not have equal 
            number of points."""
        self.wvls = wvls
        return
    
    def store_disort_result(self,DISORT,wvl_idx):
        """Stores results from a single DISORT run to the SRFM() object."""
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
        """Converts radiance to brightness temperature."""
        # strech wvnm to correct shape to be broadcastable.
        wvnm = self.wvnm[:,np.newaxis,np.newaxis,np.newaxis] # add new axis to wvnm
        # to match the number of uu dimensions, 0th dimension (axis 0) are the same
        
        assert wvnm.shape[0] == self.uu.shape[0], """wvnm and uu don't have the same
        shape of the first axis???"""
        
        uu_shape = self.uu.shape # tuple
        
        for num,i in enumerate(uu_shape[1:]):
            wvnm = np.repeat(wvnm,i,axis=num+1)
        
            
        self.bbt = utils.convert_spectral_radiance_to_bbt(
            self.uu, wvnm
        )
        return
    
    def convolve_with_iasi(self,filename):
        """Convolve radiance (uu) with iasi instrument line shape.
        Inputs:
            self, must contain uu
            filename - file name of the iasi.ils file (iasi instrument line shape kindly
            provided by Anu Dudhia, in RFM format.)
            assumes regular grid
        """        
        # save a copy of unconvolved spectrum
        uu_unconvolved = self.uu.copy()
        
        # read instrument line shape
        ils_x, ils_y, ils_lo, ils_hi = utils.read_ils(filename)
        
        # check if model grid is regular
        a = np.diff(self.wvnm, n=2) # calculate 2nd discrete difference
        a[a < 1e12] = 0 # remove small numbers (arising from computer precision limits)
        assert np.all(a) == False, """Wavenumber grid is 
        not regular.""" # check is all values in a are 0 (0 evaluates to False)
        
        # determine resolution from model wavenumber grid
        num = len(self.wvnm)
        lo = self.wvnm.min()
        hi = self.wvnm.max()
        res = np.round((hi-lo)/num,decimals=8) # this is inadvertedly introduces a limit
        # on the minimum resolution used in the code as 1e-8 cm-1, which should be 
        # enough though, and also this may not be the numerically most stable way to go
        
        # generate new grid for ils
        new_x = np.linspace(ils_lo,
                            ils_hi,
                            int((ils_hi-ils_lo)/res+1)
                            )
        
        # interpolate ils to new grid
        new_y = np.interp(new_x, ils_x, ils_y)
        
        # calculate sum of instrument line shape for normalization later
        norm = np.sum(new_y)
       
        # determine shape of uu from DISORT (basically a set of output spectra a 
        # different optical dpeths, polar and azimuthal angles
        uu_shape = self.uu.shape # tuple
        
        # determine all combinations of indices of uu
        combs = []
        for i in range(uu_shape[1]):
            for ii in range(uu_shape[2]):
                for iii in range(uu_shape[3]):
                    combs.append([i,ii,iii])
        
        # convolve spectra in a loop
        for c in combs:
            self.uu[:,c[0],c[1],c[2]] = convolve(uu_unconvolved[:,c[0],c[1],c[2]],new_y,mode="same")/norm

        return
        
    def interp(self,new_wvnm):
        """Interpolates radiances (uu) and brightness temperatures to a new grid.
        inputs:
            self - must contain uu
                 - if also contains bbt, bbt is interpolated as well
            new_wvnm - new wavenumber grid to interpolate to, units [cm-1]
        outputs:
            self.uu_interp
            (conditional) self.bb_interp
            self.wvnm, self.wvls - new grid
        """
        # strech wvnm to correct shape to be broadcastable.
#        old_wvnm = self.wvnm[:,np.newaxis,np.newaxis,np.newaxis] # add new axis to wvnm
        # to match the number of uu dimensions, 0th dimension (axis 0) are the same
        
#        n_wvnm = new_wvnm[:,np.newaxis,np.newaxis,np.newaxis]
        
        uu_shape = self.uu.shape # tuple
        
        new_uu_shape = list(uu_shape)
        new_uu_shape[0] = len(new_wvnm)
        new_uu = np.zeros(tuple(new_uu_shape))
        
        # determine all combinations of indices of uu
        combs = []
        for i in range(uu_shape[1]):
            for ii in range(uu_shape[2]):
                for iii in range(uu_shape[3]):
                    combs.append([i,ii,iii])
        
        # interpolate spectra in a loop
        for c in combs:
            new_uu[:,c[0],c[1],c[2]] = np.interp(new_wvnm,self.wvnm,self.uu[:,c[0],c[1],c[2]])
        
        self.uu = new_uu
#        self.uu = np.interp(n_wvnm,old_wvnm,self.uu)
        
        # calculate new wavelengths [um]
        new_wvls = (1/new_wvnm)*1e4
        
        # assign grids to class
        self.wvnm = new_wvnm
        self.wvls = new_wvls
        
        return
        
