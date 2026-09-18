MODULE DRVLEV_SUB
CONTAINS
SUBROUTINE DRVLEV ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   12JUL26 AD Checked.
!   01AUG25 AD Use ADDLEV instead of LEVCHK. Add log messages.
!   11FEB25 AD Checked.
!   01MAY17 AD F90 conversion of inplev.for. Tested.
!
! DESCRIPTION
!   Read RFM driver table *LEV section
!   Called by RFMDRV once if LEVFLG set TRUE.
!   Reads a list of monotonically increasing altitude levels [km] for 
!   intermediate spectral outputs.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE ATMCOM_DAT, ONLY: HGTSFC, HGTTOA, SETHGT ! Atmospheric profile data
!
! SUBROUTINES
    USE ADDLEV_SUB ! Check and add output altitude level
    USE C9REAL_GEN ! Write real number as C*9 string
    USE NXTFFL_SUB ! Load next field from rfm.drv, expanding filenames
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
    LOGICAL           :: ANYLEV = .FALSE. ! T= at least one value read 
    INTEGER(I4)       :: LENGTH ! No.characters in FIELD
    REAL(R4)          :: HGTLEV ! Altitude [km] read from FIELD
    CHARACTER(LENREC) :: FIELD  ! Field extracted from driver table
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Check that altitude profile has been specified
  IF ( .NOT. SETHGT ) THEN
    ERRMSG = 'F-DRVLEV: *HGT profile has not been supplied'
    FAIL = .TRUE.
    RETURN
  END IF
!  
  CALL WRTLOG ( 'I-DRVLEV: Set output levels at: ', .TRUE. ) 
  DO
    CALL NXTFFL ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
    ANYLEV = .TRUE.
    CALL ADDLEV ( FIELD, HGTSFC, HGTTOA, HGTLEV, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    CALL WRTLOG ( C9REAL(HGTLEV), .TRUE. )
  END DO
  CALL WRTLOG ( '', .FALSE. )  
!
  IF ( .NOT. ANYLEV ) THEN
    FAIL = .TRUE. 
    ERRMSG = 'F-DRVLEV: No output levels listed in *LEV section'
  END IF
!
END SUBROUTINE DRVLEV
END MODULE DRVLEV_SUB
