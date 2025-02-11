MODULE RFMDRV_SUB
CONTAINS
SUBROUTINE RFMDRV ( FAIL, ERRMSG )
!
! VERSION
!   22JUN24 AD Checked.
!   26JUN23 AD Remove DRVNAM - use DRVOUT instead
!   20DEC17 AD F90 conversion of rfminp.for. Tested.
!
! DESCRIPTION
!   Read RFM driver table
!   Called once by RFM.     
!
! GLOBAL DATA
    USE RFMLUN_DAT, ONLY: LUNDRV ! LUN for driver file
!
! SUBROUTINES
    USE DRVATM_SUB ! Read RFM driver table *ATM section
    USE DRVCHK_SUB ! Cross-check driver table inputs
    USE DRVCIA_SUB ! Read RFM driver table *CIA section
    USE DRVDIM_SUB ! Read RFM driver table *DIM section
    USE DRVFIN_SUB ! Read RFM driver table *FIN section
    USE DRVFLG_SUB ! Read RFM driver table *FLG section
    USE DRVFOV_SUB ! Read RFM driver table *FOV section
    USE DRVGAS_SUB ! Read RFM driver table *GAS section
    USE DRVGRD_SUB ! Read RFM driver table *GRD section
    USE DRVHDR_SUB ! Read RFM driver table *HDR section
    USE DRVHIT_SUB ! Read RFM driver table *HIT section
    USE DRVILS_SUB ! Read RFM driver table *ILS section
    USE DRVJAC_SUB ! Read RFM driver table *JAC section
    USE DRVKEY_SUB ! Check section key from driver table
    USE DRVLEV_SUB ! Read RFM driver table *LEV section
    USE DRVLUT_SUB ! Read RFM driver table *LUT section
    USE DRVNTE_SUB ! Read RFM driver table *NTE section
    USE DRVOBS_SUB ! Read RFM driver table *OBS section
    USE DRVOUT_SUB ! Read RFM driver table *OUT section
    USE DRVPHY_SUB ! Read RFM driver table *PHY section
    USE DRVREJ_SUB ! Read RFM driver table *REJ section
    USE DRVSFC_SUB ! Read RFM driver table *SFC section
    USE DRVSHP_SUB ! Read RFM driver table *SHP section
    USE DRVSKP_SUB ! Skip RFM driver table section
    USE DRVSPC_SUB ! Read RFM driver table *SPC section
    USE DRVSVD_SUB ! Read RFM driver table *SVD section
    USE DRVTAN_SUB ! Read RFM driver table *TAN section
    USE DRVXSC_SUB ! Read RFM driver table *XSC section
    USE OPNFIL_SUB ! Open input file
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    CHARACTER(4) :: KEY ! Section header
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Open driver file 
  CALL OPNFIL ( LUNDRV, 'rfm.drv', FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
! Read through each section of driver file
  DO
    CALL DRVKEY ( LUNDRV, KEY, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    SELECT CASE ( KEY ) 
! Mandatory sections, in sequence (*TAN or *DIM for #6)
    CASE ('*HDR') ; CALL DRVHDR ( LUNDRV, FAIL, ERRMSG )
    CASE ('*FLG') ; CALL DRVFLG ( LUNDRV, FAIL, ERRMSG )
    CASE ('*SPC') ; CALL DRVSPC ( LUNDRV, FAIL, ERRMSG )
    CASE ('*GAS') ; CALL DRVGAS ( LUNDRV, FAIL, ERRMSG )
    CASE ('*ATM') ; CALL DRVATM ( LUNDRV, FAIL, ERRMSG )
    CASE ('*TAN') ; CALL DRVTAN ( LUNDRV, FAIL, ERRMSG )
    CASE ('*DIM') ; CALL DRVDIM ( LUNDRV, FAIL, ERRMSG )
! Optional sections
    CASE ('*CIA') ; CALL DRVCIA ( LUNDRV, FAIL, ERRMSG )
    CASE ('*FIN') ; CALL DRVFIN ( LUNDRV, FAIL, ERRMSG )
    CASE ('*FOV') ; CALL DRVFOV ( LUNDRV, FAIL, ERRMSG )
    CASE ('*GRD') ; CALL DRVGRD ( LUNDRV, FAIL, ERRMSG )
    CASE ('*HIT') ; CALL DRVHIT ( LUNDRV, FAIL, ERRMSG )
    CASE ('*ILS') ; CALL DRVILS ( LUNDRV, FAIL, ERRMSG )
    CASE ('*JAC') ; CALL DRVJAC ( LUNDRV, FAIL, ERRMSG )
    CASE ('*LEV') ; CALL DRVLEV ( LUNDRV, FAIL, ERRMSG )
    CASE ('*LUT') ; CALL DRVLUT ( LUNDRV, FAIL, ERRMSG )
    CASE ('*NTE') ; CALL DRVNTE ( LUNDRV, FAIL, ERRMSG )
    CASE ('*OBS') ; CALL DRVOBS ( LUNDRV, FAIL, ERRMSG )
    CASE ('*OUT') ; CALL DRVOUT ( LUNDRV, FAIL, ERRMSG )
    CASE ('*PHY') ; CALL DRVPHY ( LUNDRV, FAIL, ERRMSG )
    CASE ('*REJ') ; CALL DRVREJ ( LUNDRV, FAIL, ERRMSG ) 
    CASE ('*SFC') ; CALL DRVSFC ( LUNDRV, FAIL, ERRMSG ) 
    CASE ('*SHP') ; CALL DRVSHP ( LUNDRV, FAIL, ERRMSG ) 
    CASE ('*SVD') ; CALL DRVSVD ( LUNDRV, FAIL, ERRMSG )
    CASE ('*XSC') ; CALL DRVXSC ( LUNDRV, FAIL, ERRMSG )
!
    CASE ('skip') ; CALL DRVSKP ( LUNDRV, FAIL, ERRMSG ) 
    CASE ('*END') ; EXIT
    CASE DEFAULT
      FAIL = .TRUE.
      ERRMSG = 'F-RFMDRV: Unrecognised section in driver table: ' // KEY
    END SELECT
    IF ( FAIL ) RETURN
  END DO
!
  CLOSE ( LUNDRV ) 
!
! Check all expected sections have been found (use LUNDRV=-1 as special case)
  CALL DRVKEY ( -1, KEY, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN
!
  CALL DRVCHK ( FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
END SUBROUTINE RFMDRV
END MODULE RFMDRV_SUB

