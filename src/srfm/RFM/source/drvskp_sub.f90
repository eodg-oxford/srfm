MODULE DRVSKP_SUB
CONTAINS
SUBROUTINE DRVSKP ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   02NOV23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Skip RFM driver table section
!   Called by RFMDRV for every section not required 
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
!
! SUBROUTINES
    USE NXTREC_SUB ! Load next record from input file
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
    LOGICAL           :: ENDSEC ! T=end of section reached
    CHARACTER(LENREC) :: RECORD ! Record read from driver file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL WRTLOG ( 'I-DRVSKP: Ignoring section - not required' )
!
  ENDSEC = .FALSE.
  DO WHILE ( .NOT. ENDSEC )
    CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
  END DO
!
END SUBROUTINE DRVSKP
END MODULE DRVSKP_SUB
