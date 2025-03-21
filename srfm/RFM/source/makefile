# 'makefile' for RFM created by IDL program make.pro
# To create/update RFM executable: just type 'make'
# To clean up (by deleting all .o,.mod files): type 'make clean'

# Edit next 2 lines to change compiler and options
F90 = gfortran
FLAGS = -fbounds-check -Wall
# if you don't need any flags, just change to: FLAGS = 

# change anything beyond this point and you're on your own!

# use makefile wildcard to construct list of all *.f90 files
SOURCES = $(wildcard *.f90)
# use makefile substitution reference to convert *.f90 list to *.o list
OBJECTS = $(SOURCES:.f90=.o)

# rfm executable depends on every .o file
rfm : $(OBJECTS)
	$(F90) -o rfm *.o

# set pattern for files with no additional dependencies
%.o %.mod : %.f90
	$(F90) $< -c $(FLAGS)

# set patterns for files which depend on other .mod files
# note that program rfm has no associated .mod file

rfm.o : rfm.f90 kind_dat.mod flgcom_dat.mod run_id_dat.mod spccom_dat.mod \
          hdrcom_dat.mod rfmlun_dat.mod rfmdal_sub.mod rfmdrv_sub.mod \
          rfmprf_sub.mod rfmpth_sub.mod rfmspc_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

addair_sub.o addair_sub.mod : addair_sub.f90 kind_dat.mod atmcom_dat.mod \
          flgcom_dat.mod idxcon_dat.mod addgas_sub.mod addvmr_sub.mod \
          idxgas_fnc.mod prfgra_sub.mod
	$(F90) $< -c $(FLAGS)

addatm_sub.o addatm_sub.mod : addatm_sub.f90 kind_dat.mod atmcom_dat.mod \
          gracom_dat.mod jaccom_dat.mod tancom_dat.mod obscom_dat.mod \
          qfncom_dat.mod ibrakt_gen.mod val1di_gen.mod
	$(F90) $< -c $(FLAGS)

addgas_sub.o addgas_sub.mod : addgas_sub.f90 kind_dat.mod flgcom_dat.mod \
          gascom_dat.mod idxcon_dat.mod shpcon_dat.mod c11int_fnc.mod \
          isolst_sub.mod molidx_sub.mod
	$(F90) $< -c $(FLAGS)

addgfl_sub.o addgfl_sub.mod : addgfl_sub.f90 kind_dat.mod gflcom_dat.mod
	$(F90) $< -c $(FLAGS)

addgra_sub.o addgra_sub.mod : addgra_sub.f90 kind_dat.mod atmcom_dat.mod \
          gracom_dat.mod
	$(F90) $< -c $(FLAGS)

addnte_sub.o addnte_sub.mod : addnte_sub.f90 kind_dat.mod gascom_dat.mod \
          ntecom_dat.mod atmcom_dat.mod addvib_sub.mod c11int_fnc.mod \
          idxgas_fnc.mod idxnte_fnc.mod idxqfn_fnc.mod molidx_sub.mod
	$(F90) $< -c $(FLAGS)

addpth_sub.o addpth_sub.mod : addpth_sub.f90 kind_dat.mod flgcom_dat.mod \
          jaccom_dat.mod pthcom_dat.mod iclpth_fnc.mod
	$(F90) $< -c $(FLAGS)

addqal_sub.o addqal_sub.mod : addqal_sub.f90 kind_dat.mod qalcom_dat.mod
	$(F90) $< -c $(FLAGS)

addqfn_sub.o addqfn_sub.mod : addqfn_sub.f90 kind_dat.mod ntecom_dat.mod \
          qfncom_dat.mod atmcom_dat.mod idxqfn_fnc.mod
	$(F90) $< -c $(FLAGS)

addtan_sub.o addtan_sub.mod : addtan_sub.f90 kind_dat.mod tancom_dat.mod
	$(F90) $< -c $(FLAGS)

addvib_sub.o addvib_sub.mod : addvib_sub.f90 kind_dat.mod atmcom_dat.mod \
          gracom_dat.mod
	$(F90) $< -c $(FLAGS)

addvmr_sub.o addvmr_sub.mod : addvmr_sub.f90 kind_dat.mod atmcom_dat.mod \
          gascom_dat.mod gracom_dat.mod flgcom_dat.mod
	$(F90) $< -c $(FLAGS)

adjcom_dat.o adjcom_dat.mod : adjcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

adjust_sub.o adjust_sub.mod : adjust_sub.f90 kind_dat.mod adjcom_dat.mod \
          hitcom_dat.mod setcom_dat.mod idxcon_dat.mod phycon_dat.mod \
          flgcom_dat.mod nteclc_sub.mod qtfct_fnc.mod ymix_fnc.mod
	$(F90) $< -c $(FLAGS)

atmaux_sub.o atmaux_sub.mod : atmaux_sub.f90 kind_dat.mod atmcom_dat.mod \
          gracom_dat.mod flgcom_dat.mod phyadj_dat.mod phycon_dat.mod \
          refrac_fnc.mod
	$(F90) $< -c $(FLAGS)

atmchk_sub.o atmchk_sub.mod : atmchk_sub.f90 kind_dat.mod atmcom_dat.mod \
          flgcom_dat.mod tancom_dat.mod c9real_gen.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

atmcom_dat.o atmcom_dat.mod : atmcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

atmfil_sub.o atmfil_sub.mod : atmfil_sub.f90 kind_dat.mod atmcom_dat.mod \
          rfmlun_dat.mod atmgrd_sub.mod atmprf_sub.mod atmqal_sub.mod \
          c11int_fnc.mod c9real_gen.mod chkprf_sub.mod nxtprf_sub.mod \
          opnfil_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

atmgra_sub.o atmgra_sub.mod : atmgra_sub.f90 kind_dat.mod gracom_dat.mod \
          atmcom_dat.mod intgra_sub.mod movgra_sub.mod
	$(F90) $< -c $(FLAGS)

atmgrd_sub.o atmgrd_sub.mod : atmgrd_sub.f90 kind_dat.mod flgcom_dat.mod \
          atmini_sub.mod atmlay_sub.mod
	$(F90) $< -c $(FLAGS)

atmini_sub.o atmini_sub.mod : atmini_sub.f90 kind_dat.mod atmcom_dat.mod \
          gascom_dat.mod flgcom_dat.mod
	$(F90) $< -c $(FLAGS)

atmlay_sub.o atmlay_sub.mod : atmlay_sub.f90 kind_dat.mod c11int_fnc.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

atmlev_sub.o atmlev_sub.mod : atmlev_sub.f90 kind_dat.mod atmcom_dat.mod \
          addatm_sub.mod c11int_fnc.mod c9real_gen.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

atmpar_sub.o atmpar_sub.mod : atmpar_sub.f90 kind_dat.mod atmcom_dat.mod \
          flgcom_dat.mod atmini_sub.mod atmprf_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

atmprf_sub.o atmprf_sub.mod : atmprf_sub.f90 kind_dat.mod atmcom_dat.mod \
          flgcom_dat.mod idxcon_dat.mod sfccom_dat.mod addgas_sub.mod \
          addnte_sub.mod addvmr_sub.mod chkgas_sub.mod chkprf_sub.mod \
          idxgas_fnc.mod idxnte_fnc.mod interp_gen.mod prfgra_sub.mod \
          usemol_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

atmpsi_sub.o atmpsi_sub.mod : atmpsi_sub.f90 kind_dat.mod flgcom_dat.mod \
          gracom_dat.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

atmqal_sub.o atmqal_sub.mod : atmqal_sub.f90 kind_dat.mod upcase_fnc.mod
	$(F90) $< -c $(FLAGS)

bbfwgt_fnc.o bbfwgt_fnc.mod : bbfwgt_fnc.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

binfil_sub.o binfil_sub.mod : binfil_sub.f90 kind_dat.mod c11int_fnc.mod \
          setifp_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

bright_fnc.o bright_fnc.mod : bright_fnc.f90 kind_dat.mod phycon_dat.mod
	$(F90) $< -c $(FLAGS)

c11int_fnc.o c11int_fnc.mod : c11int_fnc.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

c9real_gen.o c9real_gen.mod : c9real_gen.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

chico2_fnc.o chico2_fnc.mod : chico2_fnc.f90 kind_dat.mod chidat_dat.mod
	$(F90) $< -c $(FLAGS)

chidat_dat.o chidat_dat.mod : chidat_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

chishp_sub.o chishp_sub.mod : chishp_sub.f90 kind_dat.mod adjcom_dat.mod \
          flgcom_dat.mod phycon_dat.mod chico2_fnc.mod humlck_sub.mod
	$(F90) $< -c $(FLAGS)

chkabs_sub.o chkabs_sub.mod : chkabs_sub.f90 kind_dat.mod gascom_dat.mod \
          lflcom_dat.mod sflcom_dat.mod spccom_dat.mod
	$(F90) $< -c $(FLAGS)

chkco2_sub.o chkco2_sub.mod : chkco2_sub.f90 kind_dat.mod flgcom_dat.mod \
          gascom_dat.mod shpcon_dat.mod idxcon_dat.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

chkflx_sub.o chkflx_sub.mod : chkflx_sub.f90 kind_dat.mod atmcom_dat.mod \
          c9real_gen.mod
	$(F90) $< -c $(FLAGS)

chkfov_sub.o chkfov_sub.mod : chkfov_sub.f90 kind_dat.mod fovcom_dat.mod \
          tancom_dat.mod flgcom_dat.mod addtan_sub.mod c11int_fnc.mod \
          chklim_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

chkgas_sub.o chkgas_sub.mod : chkgas_sub.f90 kind_dat.mod idgstr_sub.mod \
          isolst_sub.mod locase_fnc.mod molidx_sub.mod reaqal_sub.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

chkhom_sub.o chkhom_sub.mod : chkhom_sub.f90 kind_dat.mod c9real_gen.mod
	$(F90) $< -c $(FLAGS)

chklev_sub.o chklev_sub.mod : chklev_sub.f90 kind_dat.mod levcom_dat.mod \
          tancom_dat.mod atmcom_dat.mod flgcom_dat.mod obscom_dat.mod \
          addtan_sub.mod c11int_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

chklim_sub.o chklim_sub.mod : chklim_sub.f90 kind_dat.mod tancom_dat.mod \
          flgcom_dat.mod atmcom_dat.mod obscom_dat.mod phyadj_dat.mod \
          c9real_gen.mod gracnv_sub.mod ibrakt_gen.mod tancnv_sub.mod
	$(F90) $< -c $(FLAGS)

chknad_sub.o chknad_sub.mod : chknad_sub.f90 kind_dat.mod flgcom_dat.mod \
          c9real_gen.mod
	$(F90) $< -c $(FLAGS)

chknam_sub.o chknam_sub.mod : chknam_sub.f90 kind_dat.mod flgcom_dat.mod \
          namcom_dat.mod gascom_dat.mod gracom_dat.mod spccom_dat.mod \
          tancom_dat.mod locase_fnc.mod
	$(F90) $< -c $(FLAGS)

chkprf_sub.o chkprf_sub.mod : chkprf_sub.f90 kind_dat.mod c9real_gen.mod
	$(F90) $< -c $(FLAGS)

chksfc_sub.o chksfc_sub.mod : chksfc_sub.f90 kind_dat.mod tancom_dat.mod \
          atmcom_dat.mod flgcom_dat.mod sfccom_dat.mod addtan_sub.mod \
          iniqad_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

chktan_sub.o chktan_sub.mod : chktan_sub.f90 kind_dat.mod tancom_dat.mod \
          chklim_sub.mod
	$(F90) $< -c $(FLAGS)

ciachk_sub.o ciachk_sub.mod : ciachk_sub.f90 kind_dat.mod ciacom_dat.mod \
          gascom_dat.mod idxcon_dat.mod c9real_gen.mod namgas_fnc.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

ciacom_dat.o ciacom_dat.mod : ciacom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

ciadef_sub.o ciadef_sub.mod : ciadef_sub.f90 kind_dat.mod ciacom_dat.mod \
          gascom_dat.mod lenrec_dat.mod idxcon_dat.mod ciafil_sub.mod \
          lexist_fnc.mod upcase_fnc.mod usemol_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

ciafil_sub.o ciafil_sub.mod : ciafil_sub.f90 kind_dat.mod gascom_dat.mod \
          spccom_dat.mod idxcon_dat.mod rfmlun_dat.mod addair_sub.mod \
          ciachk_sub.mod ciamol_sub.mod opnfil_sub.mod reacia_sub.mod \
          usemol_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

ciaint_fnc.o ciaint_fnc.mod : ciaint_fnc.f90 kind_dat.mod ciacom_dat.mod \
          interp_gen.mod
	$(F90) $< -c $(FLAGS)

ciamol_sub.o ciamol_sub.mod : ciamol_sub.f90 kind_dat.mod idxcon_dat.mod \
          molidx_sub.mod
	$(F90) $< -c $(FLAGS)

ciapth_sub.o ciapth_sub.mod : ciapth_sub.f90 kind_dat.mod ciacom_dat.mod \
          cipcom_dat.mod pthcom_dat.mod clccom_dat.mod phycon_dat.mod \
          iclpth_fnc.mod idxpth_fnc.mod
	$(F90) $< -c $(FLAGS)

cipcom_dat.o cipcom_dat.mod : cipcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

clccom_dat.o clccom_dat.mod : clccom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

coowgt_sub.o coowgt_sub.mod : coowgt_sub.f90 kind_dat.mod atmcom_dat.mod \
          phycon_dat.mod phyadj_dat.mod
	$(F90) $< -c $(FLAGS)

ctmc25_sub.o ctmc25_sub.mod : ctmc25_sub.f90 kind_dat.mod clccom_dat.mod \
          h2oc25_dat.mod widcom_dat.mod phycon_dat.mod lkpidx_fnc.mod
	$(F90) $< -c $(FLAGS)

ctmc32_sub.o ctmc32_sub.mod : ctmc32_sub.f90 kind_dat.mod clccom_dat.mod \
          h2omtc_dat.mod widcom_dat.mod phycon_dat.mod lkpidx_fnc.mod
	$(F90) $< -c $(FLAGS)

ctmckd_sub.o ctmckd_sub.mod : ctmckd_sub.f90 kind_dat.mod clccom_dat.mod \
          h2ockd_dat.mod widcom_dat.mod phycon_dat.mod lkpidx_fnc.mod
	$(F90) $< -c $(FLAGS)

ctmco2_dat.o ctmco2_dat.mod : ctmco2_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

ctmco2_sub.o ctmco2_sub.mod : ctmco2_sub.f90 kind_dat.mod clccom_dat.mod \
          ctmco2_dat.mod widcom_dat.mod lkpidx_fnc.mod
	$(F90) $< -c $(FLAGS)

ctmh2o_sub.o ctmh2o_sub.mod : ctmh2o_sub.f90 kind_dat.mod clccom_dat.mod \
          h2octm_dat.mod widcom_dat.mod phycon_dat.mod lkpidx_fnc.mod
	$(F90) $< -c $(FLAGS)

ctmn2_sub.o ctmn2_sub.mod : ctmn2_sub.f90 kind_dat.mod clccom_dat.mod \
          widcom_dat.mod phycon_dat.mod lkpidx_fnc.mod
	$(F90) $< -c $(FLAGS)

ctmo2_sub.o ctmo2_sub.mod : ctmo2_sub.f90 kind_dat.mod clccom_dat.mod \
          widcom_dat.mod phycon_dat.mod lkpidx_fnc.mod
	$(F90) $< -c $(FLAGS)

dimfil_sub.o dimfil_sub.mod : dimfil_sub.f90 kind_dat.mod lenrec_dat.mod \
          rfmlun_dat.mod nxtfl2_sub.mod opnfil_sub.mod sgnarr_gen.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

dimgrd_sub.o dimgrd_sub.mod : dimgrd_sub.f90 kind_dat.mod nxtfld_sub.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

dimpcg_sub.o dimpcg_sub.mod : dimpcg_sub.f90 kind_dat.mod atmcom_dat.mod \
          tabcom_dat.mod c9real_gen.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

dimplv_sub.o dimplv_sub.mod : dimplv_sub.f90 kind_dat.mod atmcom_dat.mod \
          tabcom_dat.mod c9real_gen.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

dimtab_sub.o dimtab_sub.mod : dimtab_sub.f90 kind_dat.mod atmcom_dat.mod \
          flgcom_dat.mod tabcom_dat.mod phycon_dat.mod rfmlun_dat.mod \
          c9real_gen.mod interp_gen.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

dopshp_sub.o dopshp_sub.mod : dopshp_sub.f90 kind_dat.mod adjcom_dat.mod \
          phycon_dat.mod
	$(F90) $< -c $(FLAGS)

drpval_fnc.o drpval_fnc.mod : drpval_fnc.f90 kind_dat.mod gracom_dat.mod \
          atmcom_dat.mod phycon_dat.mod val2di_fnc.mod
	$(F90) $< -c $(FLAGS)

drvatm_sub.o drvatm_sub.mod : drvatm_sub.f90 kind_dat.mod lenrec_dat.mod \
          flgcom_dat.mod atmchk_sub.mod atmfil_sub.mod atmgra_sub.mod \
          atmpar_sub.mod atmpsi_sub.mod nxtfld_sub.mod parfld_sub.mod
	$(F90) $< -c $(FLAGS)

drvchk_sub.o drvchk_sub.mod : drvchk_sub.f90 flgcom_dat.mod tancom_dat.mod \
          atmaux_sub.mod chkabs_sub.mod chkfov_sub.mod chklev_sub.mod \
          chknam_sub.mod chksfc_sub.mod chktan_sub.mod coowgt_sub.mod
	$(F90) $< -c $(FLAGS)

drvcia_sub.o drvcia_sub.mod : drvcia_sub.f90 kind_dat.mod lenrec_dat.mod \
          ciacom_dat.mod ciadef_sub.mod ciafil_sub.mod nxtfld_sub.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvdim_sub.o drvdim_sub.mod : drvdim_sub.f90 kind_dat.mod tabcom_dat.mod \
          lenrec_dat.mod dimgrd_sub.mod dimfil_sub.mod dimpcg_sub.mod \
          dimplv_sub.mod dimtab_sub.mod endchk_sub.mod lexist_fnc.mod \
          nxtfld_sub.mod upcase_fnc.mod
	$(F90) $< -c $(FLAGS)

drvfin_sub.o drvfin_sub.mod : drvfin_sub.f90 kind_dat.mod lenrec_dat.mod \
          fincom_dat.mod flgcom_dat.mod phycon_dat.mod c11int_fnc.mod \
          c9real_gen.mod endchk_sub.mod nxtfld_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvflg_sub.o drvflg_sub.mod : drvflg_sub.f90 kind_dat.mod flgcom_dat.mod \
          nxtfld_sub.mod upcase_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvfov_sub.o drvfov_sub.mod : drvfov_sub.f90 kind_dat.mod lenrec_dat.mod \
          endchk_sub.mod fovfil_sub.mod nxtfld_sub.mod
	$(F90) $< -c $(FLAGS)

drvgas_sub.o drvgas_sub.mod : drvgas_sub.f90 kind_dat.mod lenrec_dat.mod \
          flgcom_dat.mod gascom_dat.mod idxcon_dat.mod addgas_sub.mod \
          chkco2_sub.mod gasall_sub.mod gaschk_sub.mod namgas_fnc.mod \
          nxtffl_sub.mod reaqal_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvgrd_sub.o drvgrd_sub.mod : drvgrd_sub.f90 kind_dat.mod lenrec_dat.mod \
          gflcom_dat.mod grddef_sub.mod grdfil_sub.mod nxtfld_sub.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvhdr_sub.o drvhdr_sub.mod : drvhdr_sub.f90 kind_dat.mod lenrec_dat.mod \
          hdrcom_dat.mod endchk_sub.mod nxtrec_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvhit_sub.o drvhit_sub.mod : drvhit_sub.f90 kind_dat.mod lenrec_dat.mod \
          hflcom_dat.mod hitchk_sub.mod hitfil_sub.mod hitreq_sub.mod \
          hittyp_sub.mod nxtfld_sub.mod parfld_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvils_sub.o drvils_sub.mod : drvils_sub.f90 kind_dat.mod lenrec_dat.mod \
          ilsfil_sub.mod ilsspc_sub.mod nxtfld_sub.mod
	$(F90) $< -c $(FLAGS)

drvjac_sub.o drvjac_sub.mod : drvjac_sub.f90 kind_dat.mod lenrec_dat.mod \
          flgcom_dat.mod jacalt_sub.mod jactan_sub.mod jactgt_sub.mod \
          locase_fnc.mod nxtffl_sub.mod upcase_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvkey_sub.o drvkey_sub.mod : drvkey_sub.f90 kind_dat.mod flgcom_dat.mod \
          lenrec_dat.mod tancom_dat.mod upcase_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvlev_sub.o drvlev_sub.mod : drvlev_sub.f90 kind_dat.mod lenrec_dat.mod \
          levchk_sub.mod nxtffl_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvlut_sub.o drvlut_sub.mod : drvlut_sub.f90 kind_dat.mod lenrec_dat.mod \
          lflcom_dat.mod gascom_dat.mod spccom_dat.mod lutfil_sub.mod \
          lutdef_sub.mod nxtfld_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvnte_sub.o drvnte_sub.mod : drvnte_sub.f90 kind_dat.mod lenrec_dat.mod \
          flgcom_dat.mod rfmlun_dat.mod atmgra_sub.mod atmpsi_sub.mod \
          ntedef_sub.mod ntefil_sub.mod nxtfld_sub.mod
	$(F90) $< -c $(FLAGS)

drvobs_sub.o drvobs_sub.mod : drvobs_sub.f90 kind_dat.mod lenrec_dat.mod \
          obscom_dat.mod flgcom_dat.mod c9real_gen.mod nxtfld_sub.mod \
          obschk_sub.mod parfld_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvout_sub.o drvout_sub.mod : drvout_sub.f90 kind_dat.mod lenrec_dat.mod \
          namcom_dat.mod nxtrec_sub.mod parfld_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvphy_sub.o drvphy_sub.mod : drvphy_sub.f90 kind_dat.mod phyadj_dat.mod \
          lenrec_dat.mod phycon_dat.mod c9real_gen.mod nxtrec_sub.mod \
          parfld_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvrej_sub.o drvrej_sub.mod : drvrej_sub.f90 kind_dat.mod gascom_dat.mod \
          lenrec_dat.mod rejcom_dat.mod chkgas_sub.mod idxgas_fnc.mod \
          nxtrec_sub.mod txtfld_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvsfc_sub.o drvsfc_sub.mod : drvsfc_sub.f90 kind_dat.mod lenrec_dat.mod \
          sfccom_dat.mod atmcom_dat.mod flgcom_dat.mod c9real_gen.mod \
          nxtfld_sub.mod parfld_sub.mod sfcems_sub.mod sfclev_sub.mod \
          upcase_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvshp_sub.o drvshp_sub.mod : drvshp_sub.f90 kind_dat.mod gascom_dat.mod \
          idxcon_dat.mod lenrec_dat.mod shpcon_dat.mod flgcom_dat.mod \
          nxtrec_sub.mod shpgas_sub.mod txtfld_sub.mod upcase_fnc.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvskp_sub.o drvskp_sub.mod : drvskp_sub.f90 kind_dat.mod lenrec_dat.mod \
          nxtrec_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvspc_sub.o drvspc_sub.mod : drvspc_sub.f90 kind_dat.mod flgcom_dat.mod \
          lenrec_dat.mod nxtrec_sub.mod spcchk_sub.mod spcfil_sub.mod \
          spclab_sub.mod spctxt_sub.mod txtfld_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvsvd_sub.o drvsvd_sub.mod : drvsvd_sub.f90 kind_dat.mod lenrec_dat.mod \
          sflcom_dat.mod gascom_dat.mod spccom_dat.mod nxtfld_sub.mod \
          svddef_sub.mod svdfil_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvtan_sub.o drvtan_sub.mod : drvtan_sub.f90 kind_dat.mod flgcom_dat.mod \
          lenrec_dat.mod tancom_dat.mod nxtffl_sub.mod parfld_sub.mod \
          tanchk_sub.mod tanflx_sub.mod tanlos_sub.mod tanmtx_sub.mod \
          tanuni_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

drvxsc_sub.o drvxsc_sub.mod : drvxsc_sub.f90 kind_dat.mod lenrec_dat.mod \
          xsccom_dat.mod nxtfld_sub.mod xscfil_sub.mod xscdef_sub.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

dshval_fnc.o dshval_fnc.mod : dshval_fnc.f90 kind_dat.mod atmcom_dat.mod \
          gracom_dat.mod intrvl_sub.mod
	$(F90) $< -c $(FLAGS)

endchk_sub.o endchk_sub.mod : endchk_sub.f90 kind_dat.mod nxtfld_sub.mod
	$(F90) $< -c $(FLAGS)

fincom_dat.o fincom_dat.mod : fincom_dat.f90 kind_dat.mod rfmcon_dat.mod
	$(F90) $< -c $(FLAGS)

flxatm_sub.o flxatm_sub.mod : flxatm_sub.f90 kind_dat.mod atmcom_dat.mod \
          fulcom_dat.mod levcom_dat.mod qadcom_dat.mod fincom_dat.mod \
          flgcom_dat.mod flxlay_sub.mod
	$(F90) $< -c $(FLAGS)

flxlay_sub.o flxlay_sub.mod : flxlay_sub.f90 kind_dat.mod fincom_dat.mod \
          pthcom_dat.mod atmcom_dat.mod flgcom_dat.mod ptbcon_dat.mod \
          idxpth_fnc.mod planck_fnc.mod srcbfx_fnc.mod
	$(F90) $< -c $(FLAGS)

flxpth_sub.o flxpth_sub.mod : flxpth_sub.f90 kind_dat.mod atmcom_dat.mod \
          pthcom_dat.mod tancom_dat.mod flgcom_dat.mod phycon_dat.mod \
          iclpth_fnc.mod vrtsum_sub.mod
	$(F90) $< -c $(FLAGS)

flxsfc_sub.o flxsfc_sub.mod : flxsfc_sub.f90 kind_dat.mod sfccom_dat.mod \
          phycon_dat.mod interp_gen.mod planck_fnc.mod
	$(F90) $< -c $(FLAGS)

flxspa_sub.o flxspa_sub.mod : flxspa_sub.f90 kind_dat.mod phyadj_dat.mod \
          planck_fnc.mod
	$(F90) $< -c $(FLAGS)

fovcom_dat.o fovcom_dat.mod : fovcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

fovfil_sub.o fovfil_sub.mod : fovfil_sub.f90 kind_dat.mod fovcom_dat.mod \
          flgcom_dat.mod rfmlun_dat.mod opnfil_sub.mod
	$(F90) $< -c $(FLAGS)

fulcom_dat.o fulcom_dat.mod : fulcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

gasall_sub.o gasall_sub.mod : gasall_sub.f90 kind_dat.mod spccom_dat.mod \
          addgas_sub.mod c11int_fnc.mod lstabs_sub.mod molidx_sub.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

gaschk_sub.o gaschk_sub.mod : gaschk_sub.f90 kind_dat.mod addgas_sub.mod \
          chkgas_sub.mod gasqal_sub.mod gastql_sub.mod locase_fnc.mod
	$(F90) $< -c $(FLAGS)

gascom_dat.o gascom_dat.mod : gascom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

gasqal_sub.o gasqal_sub.mod : gasqal_sub.f90 kind_dat.mod gascom_dat.mod \
          addqal_sub.mod c11int_fnc.mod reaqal_sub.mod useqal_fnc.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

gastql_sub.o gastql_sub.mod : gastql_sub.f90 kind_dat.mod gascom_dat.mod \
          shpcon_dat.mod flgcom_dat.mod idxcon_dat.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

gauqad_sub.o gauqad_sub.mod : gauqad_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

gflcom_dat.o gflcom_dat.mod : gflcom_dat.f90 kind_dat.mod namcom_dat.mod
	$(F90) $< -c $(FLAGS)

gracnv_sub.o gracnv_sub.mod : gracnv_sub.f90 kind_dat.mod obscom_dat.mod \
          tancom_dat.mod atmcom_dat.mod flgcom_dat.mod phyadj_dat.mod \
          phycon_dat.mod raygra_sub.mod
	$(F90) $< -c $(FLAGS)

gracom_dat.o gracom_dat.mod : gracom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

gradvs_sub.o gradvs_sub.mod : gradvs_sub.f90 kind_dat.mod phyadj_dat.mod \
          flgcom_dat.mod phycon_dat.mod drpval_fnc.mod dshval_fnc.mod \
          rfrval_fnc.mod
	$(F90) $< -c $(FLAGS)

grapth_sub.o grapth_sub.mod : grapth_sub.f90 kind_dat.mod atmcom_dat.mod \
          pthcom_dat.mod tancom_dat.mod phycon_dat.mod sfccom_dat.mod \
          grasum_sub.mod
	$(F90) $< -c $(FLAGS)

grasum_sub.o grasum_sub.mod : grasum_sub.f90 kind_dat.mod phyadj_dat.mod \
          phycon_dat.mod raygra_sub.mod valgra_sub.mod
	$(F90) $< -c $(FLAGS)

grdcom_dat.o grdcom_dat.mod : grdcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

grddef_sub.o grddef_sub.mod : grddef_sub.f90 kind_dat.mod lenrec_dat.mod \
          spccom_dat.mod grdfil_sub.mod lexist_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

grdfil_sub.o grdfil_sub.mod : grdfil_sub.f90 kind_dat.mod spccom_dat.mod \
          addgfl_sub.mod reagrd_sub.mod reaspc_sub.mod spcrng_sub.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

h2oc25_dat.o h2oc25_dat.mod : h2oc25_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

h2ockd_dat.o h2ockd_dat.mod : h2ockd_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

h2octm_dat.o h2octm_dat.mod : h2octm_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

h2omtc_dat.o h2omtc_dat.mod : h2omtc_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

hdbchk_sub.o hdbchk_sub.mod : hdbchk_sub.f90 kind_dat.mod hdbcom_dat.mod \
          setcom_dat.mod idxpar_fnc.mod idxset_fnc.mod lenrec_fnc.mod \
          locase_fnc.mod txtfld_sub.mod txtpar_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

hdbcom_dat.o hdbcom_dat.mod : hdbcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

hdbfil_sub.o hdbfil_sub.mod : hdbfil_sub.f90 kind_dat.mod hdbcom_dat.mod \
          hdbchk_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

hflcom_dat.o hflcom_dat.mod : hflcom_dat.f90 kind_dat.mod hdbcom_dat.mod \
          hitcom_dat.mod
	$(F90) $< -c $(FLAGS)

hgtstr_fnc.o hgtstr_fnc.mod : hgtstr_fnc.f90 kind_dat.mod tancom_dat.mod
	$(F90) $< -c $(FLAGS)

hitchk_sub.o hitchk_sub.mod : hitchk_sub.f90 kind_dat.mod gascom_dat.mod \
          hflcom_dat.mod hitcom_dat.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

hitcom_dat.o hitcom_dat.mod : hitcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

hitfil_sub.o hitfil_sub.mod : hitfil_sub.f90 kind_dat.mod hflcom_dat.mod \
          hdbcom_dat.mod rfmlun_dat.mod binfil_sub.mod c11int_fnc.mod \
          c9real_gen.mod hdbfil_sub.mod hitqal_sub.mod idxset_fnc.mod \
          lexist_fnc.mod parfil_sub.mod twoqal_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

hitqal_sub.o hitqal_sub.mod : hitqal_sub.f90 kind_dat.mod gascom_dat.mod \
          c11int_fnc.mod molidx_sub.mod
	$(F90) $< -c $(FLAGS)

hitrec_sub.o hitrec_sub.mod : hitrec_sub.f90 kind_dat.mod hflcom_dat.mod \
          hitcom_dat.mod recbin_sub.mod rechdb_sub.mod recpar_sub.mod
	$(F90) $< -c $(FLAGS)

hitreq_sub.o hitreq_sub.mod : hitreq_sub.f90 kind_dat.mod gascom_dat.mod \
          rfmcon_dat.mod spccom_dat.mod
	$(F90) $< -c $(FLAGS)

hittyp_sub.o hittyp_sub.mod : hittyp_sub.f90 kind_dat.mod rfmlun_dat.mod \
          lexist_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

hompth_sub.o hompth_sub.mod : hompth_sub.f90 kind_dat.mod atmcom_dat.mod \
          pthcom_dat.mod phycon_dat.mod iclpth_fnc.mod
	$(F90) $< -c $(FLAGS)

humlck_sub.o humlck_sub.mod : humlck_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

ibrakt_gen.o ibrakt_gen.mod : ibrakt_gen.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

iclpth_fnc.o iclpth_fnc.mod : iclpth_fnc.f90 kind_dat.mod clccom_dat.mod \
          flgcom_dat.mod pthcom_dat.mod rfmcon_dat.mod
	$(F90) $< -c $(FLAGS)

idgnew_fnc.o idgnew_fnc.mod : idgnew_fnc.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

idgold_fnc.o idgold_fnc.mod : idgold_fnc.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

idgstr_sub.o idgstr_sub.mod : idgstr_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

idxcon_dat.o idxcon_dat.mod : idxcon_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

idxgas_fnc.o idxgas_fnc.mod : idxgas_fnc.f90 kind_dat.mod gascom_dat.mod
	$(F90) $< -c $(FLAGS)

idxnte_fnc.o idxnte_fnc.mod : idxnte_fnc.f90 kind_dat.mod ntecom_dat.mod
	$(F90) $< -c $(FLAGS)

idxpar_fnc.o idxpar_fnc.mod : idxpar_fnc.f90 kind_dat.mod setcom_dat.mod
	$(F90) $< -c $(FLAGS)

idxpth_fnc.o idxpth_fnc.mod : idxpth_fnc.f90 kind_dat.mod pthcom_dat.mod \
          tancom_dat.mod
	$(F90) $< -c $(FLAGS)

idxqfn_fnc.o idxqfn_fnc.mod : idxqfn_fnc.f90 kind_dat.mod qfncom_dat.mod
	$(F90) $< -c $(FLAGS)

idxset_fnc.o idxset_fnc.mod : idxset_fnc.f90 kind_dat.mod setcom_dat.mod
	$(F90) $< -c $(FLAGS)

ilscom_dat.o ilscom_dat.mod : ilscom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

ilsfil_sub.o ilsfil_sub.mod : ilsfil_sub.f90 kind_dat.mod ilscom_dat.mod \
          flgcom_dat.mod phycon_dat.mod rfmlun_dat.mod nxtrec_sub.mod \
          opnfil_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

ilsgrd_sub.o ilsgrd_sub.mod : ilsgrd_sub.f90 kind_dat.mod grdcom_dat.mod
	$(F90) $< -c $(FLAGS)

ilsint_sub.o ilsint_sub.mod : ilsint_sub.f90 kind_dat.mod ilscom_dat.mod \
          interp_gen.mod
	$(F90) $< -c $(FLAGS)

ilsspc_sub.o ilsspc_sub.mod : ilsspc_sub.f90 kind_dat.mod ilscom_dat.mod \
          spccom_dat.mod
	$(F90) $< -c $(FLAGS)

inibin_sub.o inibin_sub.mod : inibin_sub.f90 kind_dat.mod setifp_sub.mod
	$(F90) $< -c $(FLAGS)

inihdb_sub.o inihdb_sub.mod : inihdb_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

inihit_sub.o inihit_sub.mod : inihit_sub.f90 kind_dat.mod hflcom_dat.mod \
          inibin_sub.mod inihdb_sub.mod inipar_sub.mod recbin_sub.mod \
          rechdb_sub.mod recpar_sub.mod
	$(F90) $< -c $(FLAGS)

inilbl_sub.o inilbl_sub.mod : inilbl_sub.f90 kind_dat.mod clccom_dat.mod \
          gascom_dat.mod hflcom_dat.mod widcom_dat.mod
	$(F90) $< -c $(FLAGS)

inilut_sub.o inilut_sub.mod : inilut_sub.f90 kind_dat.mod gascom_dat.mod \
          lflcom_dat.mod lutcom_dat.mod luttab_sub.mod
	$(F90) $< -c $(FLAGS)

inipar_sub.o inipar_sub.mod : inipar_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

iniqad_sub.o iniqad_sub.mod : iniqad_sub.f90 kind_dat.mod qadcom_dat.mod \
          flgcom_dat.mod phycon_dat.mod gauqad_sub.mod
	$(F90) $< -c $(FLAGS)

inistt_sub.o inistt_sub.mod : inistt_sub.f90 kind_dat.mod spccom_dat.mod \
          sttcom_dat.mod widcom_dat.mod gascom_dat.mod
	$(F90) $< -c $(FLAGS)

inisvd_sub.o inisvd_sub.mod : inisvd_sub.f90 kind_dat.mod gascom_dat.mod \
          sflcom_dat.mod svdcom_dat.mod rfmlun_dat.mod reasvd_sub.mod
	$(F90) $< -c $(FLAGS)

interp_gen.o interp_gen.mod : interp_gen.f90 kind_dat.mod val1di_gen.mod
	$(F90) $< -c $(FLAGS)

intgra_sub.o intgra_sub.mod : intgra_sub.f90 kind_dat.mod interp_gen.mod
	$(F90) $< -c $(FLAGS)

intrvl_sub.o intrvl_sub.mod : intrvl_sub.f90 kind_dat.mod ibrakt_gen.mod
	$(F90) $< -c $(FLAGS)

isolst_sub.o isolst_sub.mod : isolst_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

jacalt_sub.o jacalt_sub.mod : jacalt_sub.f90 kind_dat.mod atmcom_dat.mod \
          jaccom_dat.mod atmlev_sub.mod c9real_gen.mod namgas_fnc.mod
	$(F90) $< -c $(FLAGS)

jaccom_dat.o jaccom_dat.mod : jaccom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

jacfov_sub.o jacfov_sub.mod : jacfov_sub.f90 kind_dat.mod fovcom_dat.mod \
          jaccom_dat.mod tancom_dat.mod addtan_sub.mod
	$(F90) $< -c $(FLAGS)

jaciso_sub.o jaciso_sub.mod : jaciso_sub.f90 kind_dat.mod lflcom_dat.mod \
          addgas_sub.mod addvmr_sub.mod idxgas_fnc.mod
	$(F90) $< -c $(FLAGS)

jacpth_sub.o jacpth_sub.mod : jacpth_sub.f90 kind_dat.mod flgcom_dat.mod \
          jaccom_dat.mod pthcom_dat.mod tancom_dat.mod clccom_dat.mod \
          addpth_sub.mod addtan_sub.mod c11int_fnc.mod flxpth_sub.mod \
          grapth_sub.mod hompth_sub.mod jacfov_sub.mod jacsky_sub.mod \
          jtppth_fnc.mod limpth_sub.mod newpth_fnc.mod ptbatm_sub.mod \
          vrtpth_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

jacsky_sub.o jacsky_sub.mod : jacsky_sub.f90 kind_dat.mod jaccom_dat.mod \
          tancom_dat.mod addtan_sub.mod
	$(F90) $< -c $(FLAGS)

jactan_sub.o jactan_sub.mod : jactan_sub.f90 kind_dat.mod tancom_dat.mod \
          jacalt_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

jactgt_sub.o jactgt_sub.mod : jactgt_sub.f90 kind_dat.mod flgcom_dat.mod \
          jaccom_dat.mod chkgas_sub.mod idxgas_fnc.mod jaciso_sub.mod \
          usemol_fnc.mod
	$(F90) $< -c $(FLAGS)

jtppth_fnc.o jtppth_fnc.mod : jtppth_fnc.f90 jaccom_dat.mod pthcom_dat.mod \
          fovcom_dat.mod
	$(F90) $< -c $(FLAGS)

lenrec_dat.o lenrec_dat.mod : lenrec_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

lenrec_fnc.o lenrec_fnc.mod : lenrec_fnc.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

levchk_sub.o levchk_sub.mod : levchk_sub.f90 kind_dat.mod levcom_dat.mod \
          atmcom_dat.mod atmlev_sub.mod
	$(F90) $< -c $(FLAGS)

levcom_dat.o levcom_dat.mod : levcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

lflcom_dat.o lflcom_dat.mod : lflcom_dat.f90 kind_dat.mod lenrec_dat.mod
	$(F90) $< -c $(FLAGS)

limpth_sub.o limpth_sub.mod : limpth_sub.f90 kind_dat.mod atmcom_dat.mod \
          pthcom_dat.mod tancom_dat.mod phycon_dat.mod raysum_sub.mod
	$(F90) $< -c $(FLAGS)

linshp_sub.o linshp_sub.mod : linshp_sub.f90 kind_dat.mod flgcom_dat.mod \
          shpcon_dat.mod chishp_sub.mod dopshp_sub.mod lorshp_sub.mod \
          mixshp_sub.mod sub25w_sub.mod voishp_sub.mod vvwcor_sub.mod \
          vvwshp_sub.mod
	$(F90) $< -c $(FLAGS)

lkpidx_fnc.o lkpidx_fnc.mod : lkpidx_fnc.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

locase_fnc.o locase_fnc.mod : locase_fnc.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

lorshp_sub.o lorshp_sub.mod : lorshp_sub.f90 kind_dat.mod adjcom_dat.mod \
          phycon_dat.mod
	$(F90) $< -c $(FLAGS)

lstabs_sub.o lstabs_sub.mod : lstabs_sub.f90 kind_dat.mod optdat_dat.mod
	$(F90) $< -c $(FLAGS)

lutcom_dat.o lutcom_dat.mod : lutcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

lutdef_sub.o lutdef_sub.mod : lutdef_sub.f90 kind_dat.mod lenrec_dat.mod \
          spccom_dat.mod gascom_dat.mod lflcom_dat.mod lexist_fnc.mod \
          lutfil_sub.mod namgas_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

lutfil_sub.o lutfil_sub.mod : lutfil_sub.f90 kind_dat.mod lflcom_dat.mod \
          spccom_dat.mod rfmlun_dat.mod idxgas_fnc.mod lutinf_sub.mod \
          namgas_fnc.mod reaqal_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

lutinf_sub.o lutinf_sub.mod : lutinf_sub.f90 kind_dat.mod phycon_dat.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

lutpth_sub.o lutpth_sub.mod : lutpth_sub.f90 kind_dat.mod clccom_dat.mod \
          lutcom_dat.mod atmcom_dat.mod lutwgt_sub.mod
	$(F90) $< -c $(FLAGS)

luttab_sub.o luttab_sub.mod : luttab_sub.f90 kind_dat.mod lutcom_dat.mod \
          tabcom_dat.mod phycon_dat.mod idxgas_fnc.mod lutpth_sub.mod \
          opnfil_sub.mod subaxs_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

lutwgt_sub.o lutwgt_sub.mod : lutwgt_sub.f90 kind_dat.mod tabcom_dat.mod \
          phycon_dat.mod ibrakt_gen.mod lutwrn_sub.mod
	$(F90) $< -c $(FLAGS)

lutwrn_sub.o lutwrn_sub.mod : lutwrn_sub.f90 kind_dat.mod tabcom_dat.mod \
          c11int_fnc.mod c9real_gen.mod namgas_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

maknam_sub.o maknam_sub.mod : maknam_sub.f90 kind_dat.mod jaccom_dat.mod \
          levcom_dat.mod run_id_dat.mod spccom_dat.mod tancom_dat.mod \
          gracom_dat.mod hgtstr_fnc.mod namgas_fnc.mod
	$(F90) $< -c $(FLAGS)

mixdat_dat.o mixdat_dat.mod : mixdat_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

mixshp_sub.o mixshp_sub.mod : mixshp_sub.f90 kind_dat.mod adjcom_dat.mod \
          phycon_dat.mod humlck_sub.mod
	$(F90) $< -c $(FLAGS)

molidx_sub.o molidx_sub.mod : molidx_sub.f90 kind_dat.mod gascom_dat.mod \
          locase_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

movgra_sub.o movgra_sub.mod : movgra_sub.f90 kind_dat.mod atmcom_dat.mod \
          gracom_dat.mod
	$(F90) $< -c $(FLAGS)

namcom_dat.o namcom_dat.mod : namcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

namgas_fnc.o namgas_fnc.mod : namgas_fnc.f90 gascom_dat.mod
	$(F90) $< -c $(FLAGS)

nammol_fnc.o nammol_fnc.mod : nammol_fnc.f90 kind_dat.mod gascom_dat.mod \
          idxgas_fnc.mod namgas_fnc.mod
	$(F90) $< -c $(FLAGS)

newpth_fnc.o newpth_fnc.mod : newpth_fnc.f90 kind_dat.mod ptbcon_dat.mod \
          jaccom_dat.mod pthcom_dat.mod
	$(F90) $< -c $(FLAGS)

nteclc_sub.o nteclc_sub.mod : nteclc_sub.f90 kind_dat.mod atmcom_dat.mod \
          hitcom_dat.mod ntecom_dat.mod phycon_dat.mod idxqfn_fnc.mod \
          val1di_gen.mod qtfct_fnc.mod qtnte_fnc.mod
	$(F90) $< -c $(FLAGS)

ntecom_dat.o ntecom_dat.mod : ntecom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

ntedef_sub.o ntedef_sub.mod : ntedef_sub.f90 kind_dat.mod gascom_dat.mod \
          lenrec_dat.mod lexist_fnc.mod ntefil_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

ntefil_sub.o ntefil_sub.mod : ntefil_sub.f90 kind_dat.mod atmcom_dat.mod \
          flgcom_dat.mod rfmlun_dat.mod addnte_sub.mod addqfn_sub.mod \
          idxnte_fnc.mod interp_gen.mod nxtrec_sub.mod opnfil_sub.mod \
          prfgra_sub.mod usemol_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

nxtffl_sub.o nxtffl_sub.mod : nxtffl_sub.f90 kind_dat.mod lenrec_dat.mod \
          rfmlun_dat.mod lexist_fnc.mod nxtfl2_sub.mod nxtrec_sub.mod \
          opnfil_sub.mod txtfld_sub.mod
	$(F90) $< -c $(FLAGS)

nxtfl2_sub.o nxtfl2_sub.mod : nxtfl2_sub.f90 kind_dat.mod lenrec_dat.mod \
          nxtrec_sub.mod txtfld_sub.mod
	$(F90) $< -c $(FLAGS)

nxtfld_sub.o nxtfld_sub.mod : nxtfld_sub.f90 kind_dat.mod lenrec_dat.mod \
          nxtrec_sub.mod txtfld_sub.mod
	$(F90) $< -c $(FLAGS)

nxtprf_sub.o nxtprf_sub.mod : nxtprf_sub.f90 kind_dat.mod upcase_fnc.mod
	$(F90) $< -c $(FLAGS)

nxtrec_sub.o nxtrec_sub.mod : nxtrec_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

obschk_sub.o obschk_sub.mod : obschk_sub.f90 kind_dat.mod atmcom_dat.mod \
          tancom_dat.mod obscom_dat.mod flgcom_dat.mod atmlev_sub.mod \
          c9real_gen.mod
	$(F90) $< -c $(FLAGS)

obscom_dat.o obscom_dat.mod : obscom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

opnfil_sub.o opnfil_sub.mod : opnfil_sub.f90 kind_dat.mod lenrec_dat.mod \
          lexist_fnc.mod nxtrec_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

opnout_sub.o opnout_sub.mod : opnout_sub.f90 kind_dat.mod flgcom_dat.mod \
          namcom_dat.mod maknam_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

optdat_dat.o optdat_dat.mod : optdat_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

parfil_sub.o parfil_sub.mod : parfil_sub.f90 kind_dat.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

parfld_sub.o parfld_sub.mod : parfld_sub.f90 kind_dat.mod upcase_fnc.mod
	$(F90) $< -c $(FLAGS)

phyadj_dat.o phyadj_dat.mod : phyadj_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

phycon_dat.o phycon_dat.mod : phycon_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

planck_fnc.o planck_fnc.mod : planck_fnc.f90 kind_dat.mod phycon_dat.mod
	$(F90) $< -c $(FLAGS)

prfgra_sub.o prfgra_sub.mod : prfgra_sub.f90 kind_dat.mod atmcom_dat.mod \
          gracom_dat.mod addgra_sub.mod
	$(F90) $< -c $(FLAGS)

ptbatm_sub.o ptbatm_sub.mod : ptbatm_sub.f90 kind_dat.mod atmcom_dat.mod \
          gracom_dat.mod jaccom_dat.mod ptbcon_dat.mod atmaux_sub.mod \
          ptbprf_sub.mod
	$(F90) $< -c $(FLAGS)

ptbcon_dat.o ptbcon_dat.mod : ptbcon_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

ptbprf_sub.o ptbprf_sub.mod : ptbprf_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

pthclc_sub.o pthclc_sub.mod : pthclc_sub.f90 kind_dat.mod pthcom_dat.mod \
          tancom_dat.mod atmcom_dat.mod clccom_dat.mod iclpth_fnc.mod \
          idxpth_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

pthcom_dat.o pthcom_dat.mod : pthcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

pthwrt_sub.o pthwrt_sub.mod : pthwrt_sub.f90 kind_dat.mod atmcom_dat.mod \
          flgcom_dat.mod hdrcom_dat.mod namcom_dat.mod obscom_dat.mod \
          pthcom_dat.mod tancom_dat.mod phyadj_dat.mod phycon_dat.mod \
          rfmlun_dat.mod sfccom_dat.mod c9real_gen.mod idxpth_fnc.mod \
          maknam_sub.mod namgas_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

qadcom_dat.o qadcom_dat.mod : qadcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

qalcom_dat.o qalcom_dat.mod : qalcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

qfncom_dat.o qfncom_dat.mod : qfncom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

qtfct_fnc.o qtfct_fnc.mod : qtfct_fnc.f90 kind_dat.mod tpsdat_dat.mod \
          qtwarn_sub.mod
	$(F90) $< -c $(FLAGS)

qtnte_fnc.o qtnte_fnc.mod : qtnte_fnc.f90 kind_dat.mod phycon_dat.mod \
          qtfct_fnc.mod
	$(F90) $< -c $(FLAGS)

qtwarn_sub.o qtwarn_sub.mod : qtwarn_sub.f90 kind_dat.mod c9real_gen.mod \
          nammol_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

radlay_sub.o radlay_sub.mod : radlay_sub.f90 kind_dat.mod fincom_dat.mod \
          flgcom_dat.mod pthcom_dat.mod tancom_dat.mod atmcom_dat.mod \
          idxpth_fnc.mod planck_fnc.mod srcbfx_fnc.mod
	$(F90) $< -c $(FLAGS)

radlev_sub.o radlev_sub.mod : radlev_sub.f90 kind_dat.mod fulcom_dat.mod \
          levcom_dat.mod
	$(F90) $< -c $(FLAGS)

radmtx_sub.o radmtx_sub.mod : radmtx_sub.f90 kind_dat.mod flgcom_dat.mod \
          fulcom_dat.mod levcom_dat.mod qadcom_dat.mod fincom_dat.mod \
          tancom_dat.mod flxatm_sub.mod flxsfc_sub.mod flxspa_sub.mod
	$(F90) $< -c $(FLAGS)

radsfc_sub.o radsfc_sub.mod : radsfc_sub.f90 kind_dat.mod fincom_dat.mod \
          sfccom_dat.mod jaccom_dat.mod ptbcon_dat.mod interp_gen.mod \
          planck_fnc.mod
	$(F90) $< -c $(FLAGS)

radsky_fnc.o radsky_fnc.mod : radsky_fnc.f90 kind_dat.mod fincom_dat.mod \
          pthcom_dat.mod qadcom_dat.mod atmcom_dat.mod flgcom_dat.mod \
          phycon_dat.mod idxpth_fnc.mod planck_fnc.mod srcbfx_fnc.mod
	$(F90) $< -c $(FLAGS)

raygra_sub.o raygra_sub.mod : raygra_sub.f90 kind_dat.mod atmcom_dat.mod \
          phyadj_dat.mod phycon_dat.mod gradvs_sub.mod
	$(F90) $< -c $(FLAGS)

raysum_sub.o raysum_sub.mod : raysum_sub.f90 kind_dat.mod phyadj_dat.mod \
          phycon_dat.mod dshval_fnc.mod rfrval_fnc.mod snell_sub.mod \
          valatm_sub.mod
	$(F90) $< -c $(FLAGS)

reacia_sub.o reacia_sub.mod : reacia_sub.f90 kind_dat.mod ciacom_dat.mod \
          phycon_dat.mod idxgas_fnc.mod
	$(F90) $< -c $(FLAGS)

reacyc_sub.o reacyc_sub.mod : reacyc_sub.f90 kind_dat.mod hitcom_dat.mod \
          c11int_fnc.mod reahit_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

reagrd_sub.o reagrd_sub.mod : reagrd_sub.f90 kind_dat.mod grdcom_dat.mod \
          phycon_dat.mod rfmlun_dat.mod
	$(F90) $< -c $(FLAGS)

reahit_sub.o reahit_sub.mod : reahit_sub.f90 kind_dat.mod gascom_dat.mod \
          hitcom_dat.mod rejcom_dat.mod hitrec_sub.mod idxgas_fnc.mod \
          setnte_sub.mod useqal_fnc.mod valiso_fnc.mod
	$(F90) $< -c $(FLAGS)

realut_sub.o realut_sub.mod : realut_sub.f90 kind_dat.mod lutcom_dat.mod \
          phycon_dat.mod
	$(F90) $< -c $(FLAGS)

reaqal_sub.o reaqal_sub.mod : reaqal_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

reaspc_sub.o reaspc_sub.mod : reaspc_sub.f90 kind_dat.mod grdcom_dat.mod \
          fincom_dat.mod phycon_dat.mod rfmlun_dat.mod
	$(F90) $< -c $(FLAGS)

reasvd_sub.o reasvd_sub.mod : reasvd_sub.f90 kind_dat.mod svdcom_dat.mod \
          idxgas_fnc.mod svdgrd_sub.mod svdhdr_sub.mod svdlim_sub.mod \
          svdpth_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

reaxsc_sub.o reaxsc_sub.mod : reaxsc_sub.f90 kind_dat.mod xflcom_dat.mod \
          molidx_sub.mod
	$(F90) $< -c $(FLAGS)

recbin_sub.o recbin_sub.mod : recbin_sub.f90 kind_dat.mod hitcom_dat.mod \
          idxcon_dat.mod c11int_fnc.mod
	$(F90) $< -c $(FLAGS)

rechdb_sub.o rechdb_sub.mod : rechdb_sub.f90 kind_dat.mod hdbcom_dat.mod \
          setcom_dat.mod hitcom_dat.mod idxcon_dat.mod phycon_dat.mod
	$(F90) $< -c $(FLAGS)

recpar_sub.o recpar_sub.mod : recpar_sub.f90 kind_dat.mod hitcom_dat.mod \
          idxcon_dat.mod phycon_dat.mod
	$(F90) $< -c $(FLAGS)

refrac_fnc.o refrac_fnc.mod : refrac_fnc.f90 kind_dat.mod phycon_dat.mod
	$(F90) $< -c $(FLAGS)

rejcom_dat.o rejcom_dat.mod : rejcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

rexcom_dat.o rexcom_dat.mod : rexcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

rexpth_sub.o rexpth_sub.mod : rexpth_sub.f90 kind_dat.mod atmcom_dat.mod \
          clccom_dat.mod idxcon_dat.mod rexcom_dat.mod idxgas_fnc.mod \
          val1di_gen.mod
	$(F90) $< -c $(FLAGS)

rfmcon_dat.o rfmcon_dat.mod : rfmcon_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

rfmdal_sub.o rfmdal_sub.mod : rfmdal_sub.f90 kind_dat.mod ciacom_dat.mod \
          gascom_dat.mod ilscom_dat.mod xsccom_dat.mod
	$(F90) $< -c $(FLAGS)

rfmdrv_sub.o rfmdrv_sub.mod : rfmdrv_sub.f90 rfmlun_dat.mod drvatm_sub.mod \
          drvchk_sub.mod drvcia_sub.mod drvdim_sub.mod drvfin_sub.mod \
          drvflg_sub.mod drvfov_sub.mod drvgas_sub.mod drvgrd_sub.mod \
          drvhdr_sub.mod drvhit_sub.mod drvils_sub.mod drvjac_sub.mod \
          drvkey_sub.mod drvlev_sub.mod drvlut_sub.mod drvnte_sub.mod \
          drvobs_sub.mod drvout_sub.mod drvphy_sub.mod drvrej_sub.mod \
          drvsfc_sub.mod drvshp_sub.mod drvskp_sub.mod drvspc_sub.mod \
          drvsvd_sub.mod drvtan_sub.mod drvxsc_sub.mod opnfil_sub.mod
	$(F90) $< -c $(FLAGS)

rfmlun_dat.o rfmlun_dat.mod : rfmlun_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

rfmprf_sub.o rfmprf_sub.mod : rfmprf_sub.f90 kind_dat.mod atmcom_dat.mod \
          flgcom_dat.mod hdrcom_dat.mod namcom_dat.mod ntecom_dat.mod \
          gracom_dat.mod rfmlun_dat.mod maknam_sub.mod movgra_sub.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

rfmpth_sub.o rfmpth_sub.mod : rfmpth_sub.f90 flgcom_dat.mod ciapth_sub.mod \
          grapth_sub.mod hompth_sub.mod jacpth_sub.mod limpth_sub.mod \
          pthclc_sub.mod pthwrt_sub.mod rexpth_sub.mod tabpth_sub.mod \
          vrtpth_sub.mod
	$(F90) $< -c $(FLAGS)

rfmspc_sub.o rfmspc_sub.mod : rfmspc_sub.f90 kind_dat.mod flgcom_dat.mod \
          spccia_sub.mod spcctm_sub.mod spcdal_sub.mod spcfin_sub.mod \
          spcflx_sub.mod spcfov_sub.mod spcful_sub.mod spcgrd_sub.mod \
          spcils_sub.mod spcini_sub.mod spcint_sub.mod spcjac_sub.mod \
          spclut_sub.mod spcout_sub.mod spcrad_sub.mod spcrex_sub.mod \
          spcsvd_sub.mod spctab_sub.mod spcwid_sub.mod spcwng_sub.mod \
          spcxsc_sub.mod
	$(F90) $< -c $(FLAGS)

rfrval_fnc.o rfrval_fnc.mod : rfrval_fnc.f90 kind_dat.mod atmcom_dat.mod \
          gracom_dat.mod val1di_gen.mod val2di_fnc.mod
	$(F90) $< -c $(FLAGS)

run_id_dat.o run_id_dat.mod : run_id_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

savxsc_sub.o savxsc_sub.mod : savxsc_sub.f90 kind_dat.mod xflcom_dat.mod \
          xsccom_dat.mod triang_sub.mod
	$(F90) $< -c $(FLAGS)

setcom_dat.o setcom_dat.mod : setcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

setifp_sub.o setifp_sub.mod : setifp_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

setnte_sub.o setnte_sub.mod : setnte_sub.f90 kind_dat.mod hitcom_dat.mod \
          ntecom_dat.mod
	$(F90) $< -c $(FLAGS)

sfccom_dat.o sfccom_dat.mod : sfccom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

sfcems_sub.o sfcems_sub.mod : sfcems_sub.f90 kind_dat.mod sfccom_dat.mod \
          spccom_dat.mod rfmlun_dat.mod c9real_gen.mod lexist_fnc.mod \
          opnfil_sub.mod sgnarr_gen.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

sfclev_sub.o sfclev_sub.mod : sfclev_sub.f90 kind_dat.mod atmcom_dat.mod \
          tancom_dat.mod obscom_dat.mod atmlev_sub.mod c9real_gen.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

sflcom_dat.o sflcom_dat.mod : sflcom_dat.f90 kind_dat.mod lenrec_dat.mod
	$(F90) $< -c $(FLAGS)

sgnarr_gen.o sgnarr_gen.mod : sgnarr_gen.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

shpcon_dat.o shpcon_dat.mod : shpcon_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

shpgas_sub.o shpgas_sub.mod : shpgas_sub.f90 kind_dat.mod gascom_dat.mod \
          shpcon_dat.mod idxcon_dat.mod chkgas_sub.mod idxgas_fnc.mod
	$(F90) $< -c $(FLAGS)

snell_sub.o snell_sub.mod : snell_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

spcchk_sub.o spcchk_sub.mod : spcchk_sub.f90 kind_dat.mod flgcom_dat.mod \
          spccom_dat.mod phycon_dat.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

spccia_sub.o spccia_sub.mod : spccia_sub.f90 kind_dat.mod cipcom_dat.mod \
          fincom_dat.mod ciaint_fnc.mod
	$(F90) $< -c $(FLAGS)

spccom_dat.o spccom_dat.mod : spccom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

spcctm_sub.o spcctm_sub.mod : spcctm_sub.f90 kind_dat.mod clccom_dat.mod \
          gascom_dat.mod idxcon_dat.mod flgcom_dat.mod widcom_dat.mod \
          ctmc32_sub.mod ctmc25_sub.mod ctmckd_sub.mod ctmco2_sub.mod \
          ctmh2o_sub.mod ctmn2_sub.mod ctmo2_sub.mod
	$(F90) $< -c $(FLAGS)

spcdal_sub.o spcdal_sub.mod : spcdal_sub.f90 kind_dat.mod lutcom_dat.mod \
          svdcom_dat.mod
	$(F90) $< -c $(FLAGS)

spcfil_sub.o spcfil_sub.mod : spcfil_sub.f90 kind_dat.mod spccom_dat.mod \
          flgcom_dat.mod grdcom_dat.mod phycon_dat.mod addgfl_sub.mod \
          reagrd_sub.mod reaspc_sub.mod spcrng_sub.mod
	$(F90) $< -c $(FLAGS)

spcfin_sub.o spcfin_sub.mod : spcfin_sub.f90 kind_dat.mod clccom_dat.mod \
          fincom_dat.mod gascom_dat.mod hitcom_dat.mod idxcon_dat.mod \
          adjust_sub.mod inihit_sub.mod linshp_sub.mod reacyc_sub.mod
	$(F90) $< -c $(FLAGS)

spcflx_sub.o spcflx_sub.mod : spcflx_sub.f90 kind_dat.mod flgcom_dat.mod \
          qadcom_dat.mod fincom_dat.mod sfccom_dat.mod phyadj_dat.mod \
          flxatm_sub.mod flxsfc_sub.mod planck_fnc.mod radmtx_sub.mod \
          tramtx_sub.mod
	$(F90) $< -c $(FLAGS)

spcfov_sub.o spcfov_sub.mod : spcfov_sub.f90 kind_dat.mod fovcom_dat.mod \
          fulcom_dat.mod jaccom_dat.mod tancom_dat.mod
	$(F90) $< -c $(FLAGS)

spcful_sub.o spcful_sub.mod : spcful_sub.f90 kind_dat.mod flgcom_dat.mod \
          fulcom_dat.mod gflcom_dat.mod grdcom_dat.mod spccom_dat.mod \
          fincom_dat.mod tancom_dat.mod c9real_gen.mod reagrd_sub.mod \
          reaspc_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

spcgrd_sub.o spcgrd_sub.mod : spcgrd_sub.f90 kind_dat.mod fincom_dat.mod \
          fulcom_dat.mod widcom_dat.mod clccom_dat.mod rfmcon_dat.mod \
          ibrakt_gen.mod
	$(F90) $< -c $(FLAGS)

spcils_sub.o spcils_sub.mod : spcils_sub.f90 kind_dat.mod fulcom_dat.mod \
          spccom_dat.mod flgcom_dat.mod tancom_dat.mod ilsgrd_sub.mod \
          ilsint_sub.mod
	$(F90) $< -c $(FLAGS)

spcini_sub.o spcini_sub.mod : spcini_sub.f90 kind_dat.mod flgcom_dat.mod \
          gascom_dat.mod spccom_dat.mod widcom_dat.mod fulcom_dat.mod \
          rfmcon_dat.mod inilbl_sub.mod inilut_sub.mod inistt_sub.mod \
          inisvd_sub.mod
	$(F90) $< -c $(FLAGS)

spcint_sub.o spcint_sub.mod : spcint_sub.f90 kind_dat.mod fulcom_dat.mod \
          flgcom_dat.mod interp_gen.mod
	$(F90) $< -c $(FLAGS)

spcjac_sub.o spcjac_sub.mod : spcjac_sub.f90 kind_dat.mod fulcom_dat.mod \
          jaccom_dat.mod tancom_dat.mod
	$(F90) $< -c $(FLAGS)

spclab_sub.o spclab_sub.mod : spclab_sub.f90 kind_dat.mod spccom_dat.mod
	$(F90) $< -c $(FLAGS)

spclos_fnc.o spclos_fnc.mod : spclos_fnc.f90 kind_dat.mod fulcom_dat.mod \
          tancom_dat.mod ptbcon_dat.mod
	$(F90) $< -c $(FLAGS)

spclut_sub.o spclut_sub.mod : spclut_sub.f90 kind_dat.mod clccom_dat.mod \
          fincom_dat.mod lutcom_dat.mod realut_sub.mod
	$(F90) $< -c $(FLAGS)

spcout_sub.o spcout_sub.mod : spcout_sub.f90 kind_dat.mod flgcom_dat.mod \
          fulcom_dat.mod jaccom_dat.mod levcom_dat.mod namcom_dat.mod \
          phycon_dat.mod tancom_dat.mod bright_fnc.mod spclos_fnc.mod \
          spcwrt_sub.mod wrtstt_sub.mod
	$(F90) $< -c $(FLAGS)

spcrad_sub.o spcrad_sub.mod : spcrad_sub.f90 kind_dat.mod fincom_dat.mod \
          flgcom_dat.mod fulcom_dat.mod tancom_dat.mod atmcom_dat.mod \
          obscom_dat.mod phyadj_dat.mod planck_fnc.mod radlay_sub.mod \
          radlev_sub.mod radsfc_sub.mod radsky_fnc.mod
	$(F90) $< -c $(FLAGS)

spcrex_sub.o spcrex_sub.mod : spcrex_sub.f90 kind_dat.mod fincom_dat.mod \
          rexcom_dat.mod phycon_dat.mod
	$(F90) $< -c $(FLAGS)

spcrng_sub.o spcrng_sub.mod : spcrng_sub.f90 kind_dat.mod phycon_dat.mod \
          rfmlun_dat.mod
	$(F90) $< -c $(FLAGS)

spcsvd_sub.o spcsvd_sub.mod : spcsvd_sub.f90 kind_dat.mod clccom_dat.mod \
          fincom_dat.mod svdcom_dat.mod
	$(F90) $< -c $(FLAGS)

spctab_sub.o spctab_sub.mod : spctab_sub.f90 kind_dat.mod tabcom_dat.mod \
          namcom_dat.mod opnout_sub.mod tabhdr_sub.mod tabwrt_sub.mod
	$(F90) $< -c $(FLAGS)

spctxt_sub.o spctxt_sub.mod : spctxt_sub.f90 kind_dat.mod spccom_dat.mod
	$(F90) $< -c $(FLAGS)

spcwid_sub.o spcwid_sub.mod : spcwid_sub.f90 kind_dat.mod clccom_dat.mod \
          gascom_dat.mod hitcom_dat.mod widcom_dat.mod flgcom_dat.mod \
          idxcon_dat.mod rfmcon_dat.mod adjust_sub.mod inihit_sub.mod \
          linshp_sub.mod reahit_sub.mod widstt_sub.mod
	$(F90) $< -c $(FLAGS)

spcwng_sub.o spcwng_sub.mod : spcwng_sub.f90 kind_dat.mod clccom_dat.mod \
          fincom_dat.mod gascom_dat.mod sttcom_dat.mod widcom_dat.mod \
          flgcom_dat.mod
	$(F90) $< -c $(FLAGS)

spcwrt_sub.o spcwrt_sub.mod : spcwrt_sub.f90 kind_dat.mod rfmlun_dat.mod \
          opnout_sub.mod wrthdr_sub.mod wrtspc_sub.mod
	$(F90) $< -c $(FLAGS)

spcxsc_sub.o spcxsc_sub.mod : spcxsc_sub.f90 kind_dat.mod clccom_dat.mod \
          fincom_dat.mod gascom_dat.mod sttcom_dat.mod xsccom_dat.mod \
          flgcom_dat.mod phycon_dat.mod triint_sub.mod xscint_fnc.mod
	$(F90) $< -c $(FLAGS)

srcbfx_fnc.o srcbfx_fnc.mod : srcbfx_fnc.f90 kind_dat.mod tancom_dat.mod \
          atmcom_dat.mod phyadj_dat.mod bbfwgt_fnc.mod planck_fnc.mod \
          temval_fnc.mod
	$(F90) $< -c $(FLAGS)

sttcom_dat.o sttcom_dat.mod : sttcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

sub25w_sub.o sub25w_sub.mod : sub25w_sub.f90 kind_dat.mod adjcom_dat.mod \
          phycon_dat.mod
	$(F90) $< -c $(FLAGS)

subaxs_sub.o subaxs_sub.mod : subaxs_sub.f90 kind_dat.mod tabcom_dat.mod \
          c11int_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

svdcom_dat.o svdcom_dat.mod : svdcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

svddef_sub.o svddef_sub.mod : svddef_sub.f90 kind_dat.mod gascom_dat.mod \
          lenrec_dat.mod spccom_dat.mod sflcom_dat.mod idgold_fnc.mod \
          lexist_fnc.mod svdfil_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

svdfil_sub.o svdfil_sub.mod : svdfil_sub.f90 kind_dat.mod sflcom_dat.mod \
          spccom_dat.mod idxgas_fnc.mod svdinf_sub.mod namgas_fnc.mod \
          reaqal_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

svdgrd_sub.o svdgrd_sub.mod : svdgrd_sub.f90 kind_dat.mod fulcom_dat.mod
	$(F90) $< -c $(FLAGS)

svdhdr_sub.o svdhdr_sub.mod : svdhdr_sub.f90 kind_dat.mod idgnew_fnc.mod \
          upcase_fnc.mod
	$(F90) $< -c $(FLAGS)

svdinf_sub.o svdinf_sub.mod : svdinf_sub.f90 kind_dat.mod rfmlun_dat.mod \
          svdhdr_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

svdlim_sub.o svdlim_sub.mod : svdlim_sub.f90 kind_dat.mod c11int_fnc.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

svdpth_sub.o svdpth_sub.mod : svdpth_sub.f90 kind_dat.mod clccom_dat.mod \
          svdcom_dat.mod svdwgt_sub.mod
	$(F90) $< -c $(FLAGS)

svdwgt_sub.o svdwgt_sub.mod : svdwgt_sub.f90 kind_dat.mod phycon_dat.mod \
          svdwrn_sub.mod
	$(F90) $< -c $(FLAGS)

svdwrn_sub.o svdwrn_sub.mod : svdwrn_sub.f90 kind_dat.mod c11int_fnc.mod \
          c9real_gen.mod namgas_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

tabcom_dat.o tabcom_dat.mod : tabcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

tabhdr_sub.o tabhdr_sub.mod : tabhdr_sub.f90 kind_dat.mod fulcom_dat.mod \
          gascom_dat.mod hdrcom_dat.mod tabcom_dat.mod
	$(F90) $< -c $(FLAGS)

tabpth_sub.o tabpth_sub.mod : tabpth_sub.f90 kind_dat.mod pthcom_dat.mod \
          tabcom_dat.mod phycon_dat.mod iclpth_fnc.mod
	$(F90) $< -c $(FLAGS)

tabwrt_sub.o tabwrt_sub.mod : tabwrt_sub.f90 kind_dat.mod fincom_dat.mod \
          pthcom_dat.mod tabcom_dat.mod
	$(F90) $< -c $(FLAGS)

tanchk_sub.o tanchk_sub.mod : tanchk_sub.f90 kind_dat.mod tancom_dat.mod \
          flgcom_dat.mod atmcom_dat.mod phycon_dat.mod atmlev_sub.mod \
          chkflx_sub.mod chkhom_sub.mod chknad_sub.mod hgtstr_fnc.mod
	$(F90) $< -c $(FLAGS)

tancnv_sub.o tancnv_sub.mod : tancnv_sub.f90 kind_dat.mod atmcom_dat.mod \
          flgcom_dat.mod obscom_dat.mod tancom_dat.mod phyadj_dat.mod \
          phycon_dat.mod rfrval_fnc.mod tanitr_sub.mod
	$(F90) $< -c $(FLAGS)

tancom_dat.o tancom_dat.mod : tancom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

tanflx_sub.o tanflx_sub.mod : tanflx_sub.f90 kind_dat.mod tancom_dat.mod \
          atmcom_dat.mod iniqad_sub.mod
	$(F90) $< -c $(FLAGS)

tanitr_sub.o tanitr_sub.mod : tanitr_sub.f90 kind_dat.mod atmcom_dat.mod \
          flgcom_dat.mod phyadj_dat.mod dshval_fnc.mod rfrval_fnc.mod
	$(F90) $< -c $(FLAGS)

tanlos_sub.o tanlos_sub.mod : tanlos_sub.f90 kind_dat.mod jaccom_dat.mod \
          flgcom_dat.mod tancom_dat.mod
	$(F90) $< -c $(FLAGS)

tanmtx_sub.o tanmtx_sub.mod : tanmtx_sub.f90 kind_dat.mod levcom_dat.mod \
          tancom_dat.mod
	$(F90) $< -c $(FLAGS)

tanuni_sub.o tanuni_sub.mod : tanuni_sub.f90 kind_dat.mod tancom_dat.mod \
          flgcom_dat.mod locase_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

temval_fnc.o temval_fnc.mod : temval_fnc.f90 kind_dat.mod gracom_dat.mod \
          atmcom_dat.mod val1di_gen.mod val2di_fnc.mod
	$(F90) $< -c $(FLAGS)

tpsdat_dat.o tpsdat_dat.mod : tpsdat_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

tramtx_sub.o tramtx_sub.mod : tramtx_sub.f90 kind_dat.mod fulcom_dat.mod \
          levcom_dat.mod qadcom_dat.mod fincom_dat.mod phycon_dat.mod \
          tancom_dat.mod
	$(F90) $< -c $(FLAGS)

triang_sub.o triang_sub.mod : triang_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

triint_sub.o triint_sub.mod : triint_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

twoqal_sub.o twoqal_sub.mod : twoqal_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

txtfld_sub.o txtfld_sub.mod : txtfld_sub.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

txtpar_fnc.o txtpar_fnc.mod : txtpar_fnc.f90 kind_dat.mod setcom_dat.mod
	$(F90) $< -c $(FLAGS)

upcase_fnc.o upcase_fnc.mod : upcase_fnc.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

usemol_fnc.o usemol_fnc.mod : usemol_fnc.f90 kind_dat.mod gascom_dat.mod
	$(F90) $< -c $(FLAGS)

useqal_fnc.o useqal_fnc.mod : useqal_fnc.f90 kind_dat.mod qalcom_dat.mod
	$(F90) $< -c $(FLAGS)

val1di_gen.o val1di_gen.mod : val1di_gen.f90 kind_dat.mod ibrakt_gen.mod
	$(F90) $< -c $(FLAGS)

val2di_fnc.o val2di_fnc.mod : val2di_fnc.f90 kind_dat.mod intrvl_sub.mod \
          val1di_gen.mod
	$(F90) $< -c $(FLAGS)

valatm_sub.o valatm_sub.mod : valatm_sub.f90 kind_dat.mod atmcom_dat.mod \
          val1di_gen.mod
	$(F90) $< -c $(FLAGS)

valgra_sub.o valgra_sub.mod : valgra_sub.f90 kind_dat.mod gracom_dat.mod \
          atmcom_dat.mod val2di_fnc.mod
	$(F90) $< -c $(FLAGS)

valiso_fnc.o valiso_fnc.mod : valiso_fnc.f90 kind_dat.mod gascom_dat.mod \
          c11int_fnc.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

voishp_sub.o voishp_sub.mod : voishp_sub.f90 kind_dat.mod adjcom_dat.mod \
          phycon_dat.mod
	$(F90) $< -c $(FLAGS)

vrtpth_sub.o vrtpth_sub.mod : vrtpth_sub.f90 kind_dat.mod atmcom_dat.mod \
          flgcom_dat.mod pthcom_dat.mod obscom_dat.mod iclpth_fnc.mod \
          vrtsum_sub.mod
	$(F90) $< -c $(FLAGS)

vrtsum_sub.o vrtsum_sub.mod : vrtsum_sub.f90 kind_dat.mod atmcom_dat.mod \
          phyadj_dat.mod phycon_dat.mod flgcom_dat.mod
	$(F90) $< -c $(FLAGS)

vvwcor_sub.o vvwcor_sub.mod : vvwcor_sub.f90 kind_dat.mod adjcom_dat.mod \
          phycon_dat.mod
	$(F90) $< -c $(FLAGS)

vvwshp_sub.o vvwshp_sub.mod : vvwshp_sub.f90 kind_dat.mod adjcom_dat.mod \
          phycon_dat.mod
	$(F90) $< -c $(FLAGS)

widcom_dat.o widcom_dat.mod : widcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

widstt_sub.o widstt_sub.mod : widstt_sub.f90 kind_dat.mod sttcom_dat.mod
	$(F90) $< -c $(FLAGS)

wrthdr_sub.o wrthdr_sub.mod : wrthdr_sub.f90 kind_dat.mod flgcom_dat.mod \
          hdrcom_dat.mod jaccom_dat.mod levcom_dat.mod spccom_dat.mod \
          tancom_dat.mod phycon_dat.mod c9real_gen.mod
	$(F90) $< -c $(FLAGS)

wrtlog_sub.o wrtlog_sub.mod : wrtlog_sub.f90 kind_dat.mod rfmlun_dat.mod
	$(F90) $< -c $(FLAGS)

wrtspc_sub.o wrtspc_sub.mod : wrtspc_sub.f90 kind_dat.mod flgcom_dat.mod \
          phycon_dat.mod
	$(F90) $< -c $(FLAGS)

wrtstt_sub.o wrtstt_sub.mod : wrtstt_sub.f90 kind_dat.mod gascom_dat.mod \
          hdrcom_dat.mod idxcon_dat.mod sttcom_dat.mod namcom_dat.mod \
          rfmlun_dat.mod opnout_sub.mod
	$(F90) $< -c $(FLAGS)

xflcom_dat.o xflcom_dat.mod : xflcom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

xsccom_dat.o xsccom_dat.mod : xsccom_dat.f90 kind_dat.mod
	$(F90) $< -c $(FLAGS)

xscdef_sub.o xscdef_sub.mod : xscdef_sub.f90 kind_dat.mod gascom_dat.mod \
          lenrec_dat.mod shpcon_dat.mod lexist_fnc.mod xscfil_sub.mod \
          wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

xscfil_sub.o xscfil_sub.mod : xscfil_sub.f90 kind_dat.mod gascom_dat.mod \
          spccom_dat.mod rfmlun_dat.mod idxgas_fnc.mod molidx_sub.mod \
          opnfil_sub.mod reaxsc_sub.mod savxsc_sub.mod wrtlog_sub.mod
	$(F90) $< -c $(FLAGS)

xscint_fnc.o xscint_fnc.mod : xscint_fnc.f90 kind_dat.mod xsccom_dat.mod
	$(F90) $< -c $(FLAGS)

ymix_fnc.o ymix_fnc.mod : ymix_fnc.f90 kind_dat.mod hitcom_dat.mod \
          mixdat_dat.mod
	$(F90) $< -c $(FLAGS)

# define what happens for 'make clean'
.PHONY: clean
clean:
	rm *.o *.mod
