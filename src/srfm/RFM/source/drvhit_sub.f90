MODULE DRVHIT_SUB
CONTAINS
SUBROUTINE DRVHIT ( LUNDRV, FAIL, ERRMSG )
!
! VERSION 
!   06AUG24 AD Checked.
!   11AUG23 AD Tested.
!   19JUN23 AD Allow for multiple files. Remove LUNHIT.
!   22MAY23 AD Allow for PARAM=VALUE construction 
!   29JAN20 AD Use HITTYP, TYPCHK to identify type of HITRAN file. Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Read RFM driver table *HIT section
!   Called by RFMDRV once if *HIT section found in driver table.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE HFLCOM_DAT, ONLY: MAXIDM, NHFL ! Max IDM value & No. of HITRAN files
!
! SUBROUTINES
    USE HITCHK_SUB ! Check molecule assignments for HITRAN files
    USE HITFIL_SUB ! Open HITRAN line data file and check contents
    USE HITREQ_SUB ! Set Wno range and list of molecules for HITRAN line data
    USE HITTYP_SUB ! Identify type of HITRAN data file
    USE NXTFLD_SUB ! Load next field from section of driver file
    USE PARFLD_SUB ! Extract Parameter=Value string from record
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for Driver File
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL           :: GOTPAR ! T= found PARAM=VALUE fields
    INTEGER(I4)       :: LENGTH ! Length of field in Driver Table
    CHARACTER(LENREC) :: FIELD  ! Field read from driver file
    CHARACTER(LENREC) :: NAMHIT ! Name of HITRAN file 
    CHARACTER(6)      :: TYPHIT ! 'BINFIL', 'PARFIL' or 'HDBFIL' expected
    LOGICAL           :: IDMREQ(MAXIDM) ! T=use molec# from file
    REAL(R8)          :: WNOREQ(2)      ! Required wavenumber range [cm-1]
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Established required wavenumber range and list of molecules
  CALL HITREQ ( WNOREQ, IDMREQ ) 
!
  DO
! Read name of HITRAN binary file
    CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
!
! New version - check for PARAM=VALUE form
    CALL PARFLD ( FIELD, GOTPAR, TYPHIT, NAMHIT ) 
!
    IF ( GOTPAR ) THEN
      IF ( TYPHIT .NE. 'BINFIL' .AND. &
           TYPHIT .NE. 'HDBFIL' .AND. &
           TYPHIT .NE. 'PARFIL'         ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-DRVHIT: Unrecognised HITRAN file type: ' // TYPHIT
        RETURN
      END IF    
    ELSE ! from here continue with old version ...
      NAMHIT = FIELD(1:LENGTH)
!
! Establish type of HITRAN data - set in HFLCOM
      CALL HITTYP ( NAMHIT, TYPHIT, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    END IF
!
! Open HITRAN file and check contents
    CALL HITFIL ( NAMHIT, TYPHIT, WNOREQ, IDMREQ, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
!
  END DO
!
  CALL HITCHK
!
  IF ( NHFL .EQ. 0 ) &
    CALL WRTLOG ( 'W-DRVHIT: No usable HITRAN line data files' )
!
END SUBROUTINE DRVHIT
END MODULE DRVHIT_SUB
