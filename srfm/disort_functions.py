"""
Name: disort_functions
Parent package: srfm
Author: Antonin Knizek
Contributors: 
Date: 18 February 2025
Purpose: Provides functions that enable the user to work with DISORT (which itself is
in Fortran 77/90) from a python interface.
""" 
import numpy as np
import warnings
from . import utilities as utils


def update_maxcmu(maxcmu):
    if (maxcmu % 2) != 0:
        maxcmu = maxcmu + 1
        print(f"maxcmu has been updated to {maxcmu}.")
    return maxcmu


def update_maxulv(maxulv, usrtau, maxcly):
    if usrtau == False:
        maxulv = maxcly + 1
        print(f"maxulv has been updated to {maxulv}.")
    return maxulv


def update_maxumu(maxumu, usrang, maxcmu):
    if usrang == False:
        maxumu = maxcmu
        print(f"maxumu has been updated to {maxumu}")
    return maxumu

def test_disort_input_format(
    maxcly,
    maxmom,
    maxcmu,
    maxumu,
    maxphi,
    maxulv,
    usrang,
    usrtau,
    ibcnd,
    onlyfl,
    prnt,
    plank,
    lamber,
    deltamplus,
    do_pseudo_sphere,
    dtauc,
    ssalb,
    pmom,
    temper,
    wvnmlo,
    wvnmhi,
    utau,
    umu0,
    phi0,
    umu,
    phi,
    fbeam,
    fisot,
    albedo,
    btemp,
    ttemp,
    temis,
    earth_radius,
    h_lyr,
    rhoq,
    rhou,
    rho_accurate,
    bemst,
    emust,
    accur,
    header,
    rfldir,
    rfldn,
    flup,
    dfdt,
    uavg,
    uu,
    albmed,
    trnmed,
):

    passmark = False
    # test maxcly
    if not isinstance(maxcly, int):
        raise TypeError("maxcly (or nlyr) must be type integer.")
    # test maxmom
    if not isinstance(maxmom, int):
        raise TypeError("maxmom (or nmom) must be type integer.")
    # test maxcmu
    if not isinstance(maxcmu, int):
        raise TypeError("maxcmu (or nstr) must be type integer.")

    # test maxumu
    if not isinstance(maxumu, int):
        raise TypeError("maxumu (or numu) must be type integer.")

    # test maxphi
    if not isinstance(maxphi, int):
        raise TypeError("maxphi (or nphi) must be type integer.")

    # test maxulv
    if not isinstance(maxulv, int):
        raise TypeError("maxulv (or ntau) must be type integer.")

    # test usrang
    if not isinstance(usrang, bool):
        raise TypeError("usrang must be type boolean.")

    # test usrtau
    if not isinstance(usrtau, bool):
        raise TypeError("usrtau must be type boolean.")

    # test ibcnd
    if not isinstance(ibcnd, int):
        raise TypeError("ibcnd must be type integer.")

    # test onlyfl
    if not isinstance(onlyfl, bool):
        raise TypeError("onlyfl must be type bool.")

    # test prnt
    if not isinstance(prnt, list):
        raise TypeError("prnt must be a list.")
    if not len(prnt) == 5:
        raise ValueError("prnt must a a list of length 5.")
    for i in prnt:
        if not isinstance(i, bool):
            raise ValueError("All elements in prnt must be type boolean.")

    # test plank
    if not isinstance(plank, bool):
        raise TypeError("plank must be type bool.")

    # test lamber
    if not isinstance(lamber, bool):
        raise TypeError("lamber must be type bool.")

    # test deltamplus
    if not isinstance(deltamplus, bool):
        raise TypeError("deltamplus must be type bool.")

    # test do_pseudo_sphere
    if not isinstance(do_pseudo_sphere, bool):
        raise TypeError("do_pseudo_sphere must be type bool.")

    # test dtauc
    if not isinstance(dtauc, (list, np.ndarray)):
        raise TypeError("dtauc must be type list or np.ndarray.")
    if isinstance(dtauc, list):
        if len(dtauc) != maxcly:
            raise ValueError(
                "if dtauc is a list, it must satisfy len(dtauc) == maxcly."
            )
        for i in dtauc:
            if type(i) != int and type(i) != float:
                raise ValueError(
                    "if dtauc is a list, it can only contain integers and floats"
                )

    if isinstance(dtauc, np.ndarray):
        if dtauc.shape != (maxcly,):
            raise ValueError("if dtauc is a np.ndarray, it must have shape (maxcly,).")
        for i in dtauc:
            if type(i.item()) != float and type(i.item()) != int:
                raise ValueError(
                    "if dtauc is a np.ndarray, its elements must be integers or floats."
                )

    if isinstance(dtauc, list):
        max_dtauc = max(dtauc)
    elif isinstance(dtauc, np.ndarray):
        max_dtauc = np.max(dtauc).item()

    # test ssalb
    if not isinstance(ssalb, (list, np.ndarray)):
        raise TypeError("ssalb must be type list or np.ndarray.")
    if isinstance(ssalb, list) and not len(ssalb) == maxcly:
        raise ValueError("if ssalb is a list, it must satisfy len(ssalb) == maxcly.")
    if isinstance(ssalb, np.ndarray) and not ssalb.shape == (maxcly,):
        raise ValueError("if ssalb is a np.ndarray, it must have shape (maxcly,).")

    # test pmom
    if not isinstance(pmom, np.ndarray):
        raise TypeError("pmom must be type np.ndarray.")
    if not pmom.shape == (maxmom + 1, maxcly):
        raise ValueError("pmom must have shape (maxmom+1,maxcly).")
    for i in pmom:
        for j in i:
            if not isinstance(j.item(), (int, float)):
                raise ValueError("pmom may only contain integers and floats.")

    # test temper
    if isinstance(temper, list):
        if len(temper) != maxcly + 1:
            raise ValueError(
                "if temper is a list, it must satisfy len(temper) == maxcly+1"
            )
        for i in temper:
            if not isinstance(i, (int, float)):
                raise ValueError(
                    "if temper is a list, it may contain only integers and floats."
                )
    elif isinstance(temper, np.ndarray):
        if temper.shape != (maxcly + 1,):
            raise ValueError(
                "if temper is a np.ndarray, it must have shape (maxcly+1,)"
            )
        for i in temper:
            if not isinstance(i.item(), (int, float)):
                raise ValueError(
                    "if temper is a np.ndarray, it may contain only integers and floats."
                )
    else:
        raise TypeError("temper must be a list or a np.ndarray.")
    #
    #    #test temper
    #    if not isinstance(temper, np.ndarray):
    #        raise TypeError('temper must be type nd.ndarray.')
    #    if not temper.shape == (maxcly+1,):
    #        raise ValueError('temper must have shape (maxcly+1,).')

    # test wvnmlo
    if not isinstance(wvnmlo, (int, float)):
        raise TypeError("wvnmlo must be type integer.")

    # test wvnmhi
    if not isinstance(wvnmhi, (int, float)):
        raise TypeError("wvnmhi must be type integer.")

    # test utau
    if not isinstance(utau, (list, np.ndarray)):
        raise TypeError("utau must be type list or np.ndarray.")
    if isinstance(utau, list):
        if len(utau) != maxulv:
            raise ValueError("if utau is a list, it must satisfy len(utau) == maxulv.")
        for i in utau:
            if type(i) != int and type(i) != float:
                raise ValueError(
                    "if utau is a list, it must contain only integers and floats."
                )

    if isinstance(utau, np.ndarray):
        if utau.shape != (maxulv,):
            raise ValueError(
                "if utau is a np.ndarray, it must have shape == (maxulv,)"
            )
        for i in utau:
            if type(i.item()) != int and type(i.item()) != float:
                raise ValueError(
                    "if utau is a np.ndarray, its elements must be integers or floats."
                )

    # test umu0
    if not isinstance(umu0, (float, int)):
        raise TypeError("umu0 must be an integer or a float.")

    # test phi0
    if not isinstance(phi0, (float, int)):
        raise TypeError("phi0 must be an integer or a float.")

    # test umu
    if isinstance(umu, list):
        if len(umu) != maxumu:
            raise ValueError("if umu is a list, it must satisfy len(umu) == maxumu")
        for i in umu:
            if not isinstance(i, (int, float)):
                raise ValueError(
                    "if umu is a list, it may contain only integers and floats."
                )
    elif isinstance(umu, np.ndarray):
        if umu.shape != (maxumu,):
            raise ValueError("if umu is a np.ndarray, it must have shape (maxumu,)")
        for i in umu:
            if not isinstance(i.item(), (int, float)):
                raise ValueError(
                    "if umu is a np.ndarray, it may contain only integers and floats."
                )
    else:
        raise TypeError("umu must be a list or a np.ndarray.")

    # test phi
    if isinstance(phi, list):
        if len(phi) != maxphi:
            raise ValueError("if phi is a list, it must satisfy len(phi) == maxphi")
        for i in phi:
            if not isinstance(i, (int, float)):
                raise ValueError(
                    "if phi is a list, it may contain only integers and floats."
                )
    elif isinstance(phi, np.ndarray):
        if phi.shape != (maxphi,):
            raise ValueError("if phi is a np.ndarray, it must have shape (maxphi,)")
        for i in phi:
            if not isinstance(i.item(), (int, float)):
                raise ValueError(
                    "if phi is a np.ndarray, it may contain only integers and floats."
                )
    else:
        raise TypeError("phi must be a list or a np.ndarray.")

    # test fbeam
    if not isinstance(fbeam, (float, int)):
        raise TypeError("fbeam must be an integer or a float.")

    # test fisot
    if not isinstance(fisot, (float, int)):
        raise TypeError("fisot must be an integer or a float.")

    # test albedo
    if not isinstance(albedo, (float, int)):
        raise TypeError("albedo must be an integer or a float.")

    # test btemp
    if not isinstance(btemp, (float, int)):
        raise TypeError("btemp must be an integer or a float.")

    # test ttemp
    if not isinstance(ttemp, (float, int)):
        raise TypeError("ttemp must be an integer or a float.")

    # test temis
    if not isinstance(temis, (float, int)):
        raise TypeError("albedo must be an integer or a float.")

    # test earth_radius
    if not isinstance(earth_radius, (float, int)):
        raise TypeError("earth_radius must be an integer or a float.")
    if not earth_radius == 6371.0:
        print("Are you sure that your Earth radius is correct?")

    # test h_lyr
    if isinstance(h_lyr, list):
        if len(h_lyr) != maxcly + 1:
            raise ValueError(
                "if h_lyr is a list, it must satisfy len(h_lyr) == maxcly+1"
            )
        for i in h_lyr:
            if not isinstance(i, (int, float)):
                raise ValueError(
                    "if h_lyr is a list, it may contain only integers and floats."
                )
    elif isinstance(h_lyr, np.ndarray):
        if h_lyr.shape != (maxcly + 1,):
            raise ValueError("if h_lyr is a np.ndarray, it must have shape (maxcly+1,)")
        for i in h_lyr:
            if not isinstance(i.item(), (int, float)):
                raise ValueError(
                    "if h_lyr is a np.ndarray, it may contain only integers and floats."
                )
    else:
        raise TypeError("h_lyr must be a list or a np.ndarray.")

    # test rhoq
    if isinstance(rhoq, np.ndarray):
        if rhoq.shape != (int(maxcmu / 2), int(maxcmu / 2 + 1), maxcmu):
            raise ValueError("rhoq must have shape (maxcmu/2,maxcmu/2+1,maxcmu)")
        for i in rhoq:
            for j in i:
                for k in j:
                    if not isinstance(k.item(), (int, float)):
                        raise ValueError("rhoq may contain only integers and floats.")
    else:
        raise TypeError("rhoq must be a np.ndarray.")

    # test rhou
    if isinstance(rhou, np.ndarray):
        if rhou.shape != (maxumu, int(maxcmu / 2 + 1), maxcmu):
            raise ValueError("rhou must have shape (maxumu,maxcmu/2+1,maxcmu)")
        for i in rhou:
            for j in i:
                for k in j:
                    if not isinstance(k.item(), (int, float)):
                        raise ValueError("rhou may contain only integers and floats.")
    else:
        raise TypeError("rhou must be a np.ndarray.")

    # test rho_accurate
    if isinstance(rho_accurate, np.ndarray):
        if rho_accurate.shape != (maxumu, maxphi):
            raise ValueError("rho_accurate must have shape (maxumu,maxphi)")
        for i in rho_accurate:
            for j in i:
                if not isinstance(j.item(), (int, float)):
                    raise ValueError(
                        "rho_accurate may contain only integers and floats."
                    )
    else:
        raise TypeError("rho_accurate must be a np.ndarray.")

    # test bemst
    if isinstance(bemst, list):
        if len(bemst) != int(maxcmu / 2):
            raise ValueError(
                "if bemst is a list, it must satisfy len(bemst) == maxcmu/2"
            )
        for i in bemst:
            if not isinstance(i, (int, float)):
                raise ValueError(
                    "if bemst is a list, it may contain only integers and floats."
                )
    elif isinstance(bemst, np.ndarray):
        if bemst.shape != (int(maxcmu / 2),):
            raise ValueError("if bemst is a np.ndarray, it must have shape (maxcmu/2,)")
        for i in bemst:
            if not isinstance(i.item(), (int, float)):
                raise ValueError(
                    "if bemst is a np.ndarray, it may contain only integers and floats."
                )
    else:
        raise TypeError("bemst must be a list or a np.ndarray.")

    # test emust
    if isinstance(emust, list):
        if len(emust) != maxumu:
            raise ValueError("if emust is a list, it must satisfy len(emust) == maxumu")
        for i in emust:
            if not isinstance(i, (int, float)):
                raise ValueError(
                    "if emust is a list, it may contain only integers and floats."
                )
    elif isinstance(emust, np.ndarray):
        if emust.shape != (maxumu,):
            raise ValueError("if emust is a np.ndarray, it must have shape (maxumu,)")
        for i in emust:
            if not isinstance(i.item(), (int, float)):
                raise ValueError(
                    "if emust is a np.ndarray, it may contain only integers and floats."
                )
    else:
        raise TypeError("emust must be a list or a np.ndarray.")

    # test accur
    if not isinstance(accur, (float, int)):
        raise TypeError("accur must be an integer or a float.")

    # test header
    if isinstance(header, str):
        if len(header) > 127:
            raise ValueError("header must ahve a maximum 127 characters.")
    else:
        raise TypeError("header must be a string.")

    # test rfldir
    if isinstance(rfldir, list):
        if len(rfldir) != maxulv:
            raise ValueError(
                "if rfldir is a list, it must satisfy len(rfldir) == maxulv. "
                + "Warning: rfldir is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in rfldir:
            if not isinstance(i, (int, float)):
                raise ValueError(
                    "if rfldir is a list, it may contain only integers and floats. "
                    + "Warning: rfldir is an output variable, it's contents are not used by disort and will be overwritten."
                )
    elif isinstance(rfldir, np.ndarray):
        if rfldir.shape != (maxulv,):
            raise ValueError(
                "if rfldir is a np.ndarray, it must have shape (maxulv,). "
                + "Warning: rfldir is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in rfldir:
            if not isinstance(i.item(), (int, float)):
                raise ValueError(
                    "if rfldir is a np.ndarray, it may contain only integers and floats. "
                    + "Warning: rfldir is an output variable, it's contents are not used by disort and will be overwritten."
                )
    else:
        raise TypeError(
            "rfldir must be a list or a np.ndarray. "
            + "Warning: rfldir is an output variable, it's contents are not used by disort and will be overwritten."
        )

    # test rfldn
    if isinstance(rfldn, list):
        if len(rfldn) != maxulv:
            raise ValueError(
                "if rfldn is a list, it must satisfy len(rfldn) == maxulv. "
                + "Warning: rfldn is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in rfldn:
            if not isinstance(i, (int, float)):
                raise ValueError(
                    "if rfldn is a list, it may contain only integers and floats. "
                    + "Warning: rfldn is an output variable, it's contents are not used by disort and will be overwritten."
                )
    elif isinstance(rfldn, np.ndarray):
        if rfldn.shape != (maxulv,):
            raise ValueError(
                "if rfldn is a np.ndarray, it must have shape (rfldn,). "
                + "Warning: rfldn is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in rfldn:
            if not isinstance(i.item(), (int, float)):
                raise ValueError(
                    "if rfldn is a np.ndarray, it may contain only integers and floats. "
                    + "Warning: rfldn is an output variable, it's contents are not used by disort and will be overwritten."
                )
    else:
        raise TypeError(
            "rfldn must be a list or a np.ndarray. "
            + "Warning: rfldn is an output variable, it's contents are not used by disort and will be overwritten."
        )

    # test flup
    if isinstance(flup, list):
        if len(flup) != maxulv:
            raise ValueError(
                "if flup is a list, it must satisfy len(flup) == maxulv. "
                + "Warning: flup is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in flup:
            if not isinstance(i, (int, float)):
                raise ValueError(
                    "if flup is a list, it may contain only integers and floats. "
                    + "Warning: flup is an output variable, it's contents are not used by disort and will be overwritten."
                )
    elif isinstance(flup, np.ndarray):
        if flup.shape != (maxulv,):
            raise ValueError(
                "if flup is a np.ndarray, it must have shape (flup,). "
                + "Warning: flup is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in flup:
            if not isinstance(i.item(), (int, float)):
                raise ValueError(
                    "if flup is a np.ndarray, it may contain only integers and floats. "
                    + "Warning: flup is an output variable, it's contents are not used by disort and will be overwritten."
                )
    else:
        raise TypeError(
            "flup must be a list or a np.ndarray. "
            + "Warning: flup is an output variable, it's contents are not used by disort and will be overwritten."
        )

    # test dfdt
    if isinstance(dfdt, list):
        if len(dfdt) != maxulv:
            raise ValueError(
                "if dfdt is a list, it must satisfy len(dfdt) == maxulv. "
                + "Warning: dfdt is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in dfdt:
            if not isinstance(i, (int, float)):
                raise ValueError(
                    "if dfdt is a list, it may contain only integers and floats. "
                    + "Warning: dfdt is an output variable, it's contents are not used by disort and will be overwritten."
                )
    elif isinstance(dfdt, np.ndarray):
        if dfdt.shape != (maxulv,):
            raise ValueError(
                "if dfdt is a np.ndarray, it must have shape (dfdt,). "
                + "Warning: dfdt is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in dfdt:
            if not isinstance(i.item(), (int, float)):
                raise ValueError(
                    "if dfdt is a np.ndarray, it may contain only integers and floats. "
                    + "Warning: dfdt is an output variable, it's contents are not used by disort and will be overwritten."
                )
    else:
        raise TypeError(
            "dfdt must be a list or a np.ndarray. "
            + "Warning: dfdt is an output variable, it's contents are not used by disort and will be overwritten."
        )

    # test uavg
    if isinstance(uavg, list):
        if len(uavg) != maxulv:
            raise ValueError(
                "if uavg is a list, it must satisfy len(uavg) == maxulv. "
                + "Warning: uavg is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in uavg:
            if not isinstance(i, (int, float)):
                raise ValueError(
                    "if uavg is a list, it may contain only integers and floats. "
                    + "Warning: uavg is an output variable, it's contents are not used by disort and will be overwritten."
                )
    elif isinstance(uavg, np.ndarray):
        if uavg.shape != (maxulv,):
            raise ValueError(
                "if uavg is a np.ndarray, it must have shape (uavg,). "
                + "Warning: uavg is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in uavg:
            if not isinstance(i.item(), (int, float)):
                raise ValueError(
                    "if uavg is a np.ndarray, it may contain only integers and floats. "
                    + "Warning: uavg is an output variable, it's contents are not used by disort and will be overwritten."
                )
    else:
        raise TypeError(
            "uavg must be a list or a np.ndarray. "
            + "Warning: uavg is an output variable, it's contents are not used by disort and will be overwritten."
        )

    # test uu
    if isinstance(uu, np.ndarray):
        if uu.shape != (maxumu, maxulv, maxphi):
            raise ValueError("uu must have shape (maxumu,maxulv,maxphi)")
        for i in uu:
            for j in i:
                for k in j:
                    if not isinstance(k.item(), (int, float)):
                        raise ValueError("uu may contain only integers and floats.")
    else:
        raise TypeError("uu must be a np.ndarray.")

    # test albmed
    if isinstance(albmed, list):
        if len(albmed) != maxumu:
            raise ValueError(
                "if albmed is a list, it must satisfy len(albmed) == maxumu. "
                + "Warning: albmed is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in albmed:
            if not isinstance(i, (int, float)):
                raise ValueError(
                    "if albmed is a list, it may contain only integers and floats. "
                    + "Warning: albmed is an output variable, it's contents are not used by disort and will be overwritten."
                )
    elif isinstance(albmed, np.ndarray):
        if albmed.shape != (maxumu,):
            raise ValueError(
                "if albmed is a np.ndarray, it must have shape (maxumu,). "
                + "Warning: albmed is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in albmed:
            if not isinstance(i.item(), (int, float)):
                raise ValueError(
                    "if albmed is a np.ndarray, it may contain only integers and floats. "
                    + "Warning: albmed is an output variable, it's contents are not used by disort and will be overwritten."
                )
    else:
        raise TypeError(
            "albmed must be a list or a np.ndarray. "
            + "Warning: albmed is an output variable, it's contents are not used by disort and will be overwritten."
        )

    # test trnmed
    if isinstance(trnmed, list):
        if len(trnmed) != maxumu:
            raise ValueError(
                "if trnmed is a list, it must satisfy len(trnmed) == maxumu. "
                + "Warning: trnmed is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in trnmed:
            if not isinstance(i, (int, float)):
                raise ValueError(
                    "if trnmed is a list, it may contain only integers and floats. "
                    + "Warning: trnmed is an output variable, it's contents are not used by disort and will be overwritten."
                )
    elif isinstance(trnmed, np.ndarray):
        if trnmed.shape != (maxumu,):
            raise ValueError(
                "if trnmed is a np.ndarray, it must have shape (maxumu,). "
                + "Warning: trnmed is an output variable, it's contents are not used by disort and will be overwritten."
            )
        for i in trnmed:
            if not isinstance(i.item(), (int, float)):
                raise ValueError(
                    "if trnmed is a np.ndarray, it may contain only integers and floats. "
                    + "Warning: trnmed is an output variable, it's contents are not used by disort and will be overwritten."
                )
    else:
        raise TypeError(
            "trnmed must be a list or a np.ndarray. "
            + "Warning: trnmed is an output variable, it's contents are not used by disort and will be overwritten."
        )
    passmark = True
    return passmark


def test_disort_input_integrity(
    maxmom,
    maxcmu,
    maxumu,
    maxphi,
    ibcnd,
    dtauc,
    ssalb,
    temper,
    wvnmlo,
    wvnmhi,
    utau,
    umu0,
    phi0,
    umu,
    phi,
    btemp,
    ttemp,
    temis,
):

    passmark = False

    # test maxmom
    if maxmom < maxcmu:
        raise ValueError("Maxmom must be >= maxcmu.")

    # test maxcmu
    if maxcmu < 2:
        raise ValueError("Maxcmu must be >= 2.")
    elif (maxcmu % 2) != 0:
        raise ValueError("Maxcmu must be an even number.")

    # test maxumu
    if maxumu == 0 and onlyfl == False:
        raise ValueError("Maxumu can be 0 only when onlyfl == True.")

    # test maxphi
    if maxphi == 0 and onlyfl == False:
        raise ValueError("Maxphi can be 0 only when onlyfl == True.")

    # test ibcnd
    if ibcnd != 0 and ibcnd != 1:
        raise ValueError("Ibcnd can only be 0 or 1.")

    # test ssalb
    if isinstance(ssalb, list):
        for i in ssalb:
            if i < 0 or i > 1:
                raise ValueError("ssalb must be between 0 and 1.")
    elif isinstance(ssalb, np.ndarray):
        for i in ssalb:
            if i.item() < 0 or i.item() > 1:
                raise ValueError("ssalb must be between 0 and 1.")

    # test temper
    if isinstance(temper, list):
        for i in temper:
            if i < 0:
                raise ValueError("Sorry, you can't have negative temperature.")
    elif isinstance(temper, np.ndarray):
        for i in temper:
            if i.item() < 0:
                raise ValueError("Sorry, you can't have negative temperature.")

    # test wvnmlo and wvnmhi
    if wvnmlo > wvnmhi:
        raise ValueError("Wvnmlo must be < wvnmhi.")
    elif wvnmlo == wvnmhi:
        print("Wnvmlo and wvnmhi are equal, Planck function contribution will be 0.")

    # test utau
    if isinstance(utau, list):
        max_utau = max(utau)
    elif isinstance(utau, np.ndarray):
        max_utau = np.max(utau).item()

    if isinstance(dtauc, list):
        max_dtauc = max(dtauc)
    elif isinstance(dtauc, np.ndarray):
        max_dtauc = np.max(dtauc).item()

    if max_utau > max_dtauc:
        raise ValueError("Values in utau must not exceed the max value of dtauc.")

    # test umu0
    if umu0 < -1 or umu0 > 1:
        raise ValueError("Umu0 must be between -1 and 1.")

    if isinstance(umu, list):
        for i in umu:
            if i == -umu0:
                warnings.warn("Asking for output at incident beam angle, risk of 0/0, consider changing umu or umu0.",UserWarning)

    elif isinstance(umu, np.ndarray):
        for i in umu:
            if i.item() == -umu0:
                warnings.warn("Asking for output at incident beam angle, risk of 0/0, consider changing umu or umu0.",UserWarning)

    # test phi0
    if phi0 < 0 or phi0 > 360:
        raise ValueError("phi0 must be between 0 and 360.")

    # test umu
    if isinstance(umu, list):
        tumu = umu
    elif isinstance(umu, np.ndarray):
        tumu = umu.tolist()

    if 0 in tumu:
        raise ValueError("Umu must not contain 0.")
    x = []
    x.extend(tumu)
    x.sort()
    if x != tumu:
        raise ValueError("Umu must be stricly increasing.")

    # test phi
    if isinstance(phi, list):
        temp_phi = phi
    elif isinstance(phi, np.ndarray):
        temp_phi = phi.tolist()

    for i in temp_phi:
        if i < 0 or i > 360:
            raise ValueError("Phi must contain values between 0 and 360.")

    # test btemp
    if btemp < 0:
        raise ValueError("Btemp must be >= 0.")

    # test ttemp
    if ttemp < 0:
        raise ValueError("Ttemp must be >= 0.")

    # test temis
    if temis < 0:
        raise ValueError("Temis must be >= 0.")

    passmark = True
    return passmark

def get_pmom_from_disort(self, pmom, iphas, gg, nmom, lyrs=None):
    """Calculate phase function moments from disort.
    pmom - pmom array, can be all zeros
    iphas - disort getmom parameter, phase function option, 1-7
    gg - disort getmom parameter, assymetry parameter
         for Henyey_greenstein function
    nmom - index of the highest Legendre coefficient needed
           should equal maxmom (number of phase function moments)
    lyrs - layers to perform calculation for
        if None, calculate for all layers
    Returns array with shape pmom.shape with Legendre polynomial 
    coefficients for the selected phase function.
    """
    if lyrs == None:
        lyrs = list(range(pmom.shape[1]))
    if isinstance(lyrs, (list, tuple, set)):
        for i in lyrs:
            if not isinstance(i, int):
                try:
                    i = int(i)
                except TypeError:
                    print("Could not convert values in lyrs to ints.")

            if i not in range(pmom.shape[1]):
                raise ValueError("Layer number not in layers.")
            pmom[:, i] = dm.getmom(iphas=iphas, gg=gg, nmom=nmom, pmom=pmom[:, i])
            return pmom
    else:
        raise TypeError("Layers must be type list, tuple, or set.")
