MODULE DRVLEV_SUB
CONTAINS
SUBROUTINE DRVLEV ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   13FEB24 AD Checked.
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
!
! SUBROUTINES
    USE LEVCHK_SUB ! Check if string is legal altitude level
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
    INTEGER(I4)       :: LENGTH           ! No.characters in FIELD
    CHARACTER(LENREC) :: LEVSTR           ! Field extracted from driver table
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL WRTLOG ( 'I-DRVLEV: ', .TRUE. ) 
!  
  DO
    CALL NXTFFL ( LUNDRV, LEVSTR, LENGTH, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
    ANYLEV = .TRUE.
    CALL LEVCHK ( LEVSTR, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    CALL WRTLOG ( LEVSTR, .TRUE. ) 
  END DO
!
  CALL WRTLOG ( '', .FALSE. ) 
!
  IF ( .NOT. ANYLEV ) THEN
    FAIL = .TRUE. 
    ERRMSG = 'F-DRVLEV: No output levels listed in *LEV section'
  END IF
!
END SUBROUTINE DRVLEV
END MODULE DRVLEV_SUB
