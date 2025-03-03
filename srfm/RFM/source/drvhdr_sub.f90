MODULE DRVHDR_SUB
CONTAINS
SUBROUTINE DRVHDR ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   20DEC23 AD Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Read RFM driver table *HDR section
!   Called once by RFMDRV.
! 
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE HDRCOM_DAT, ONLY: TXTHDR ! Text header from driver table
!
! SUBROUTINES
    USE ENDCHK_SUB ! Check end of Driver Table section has been reached
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
    LOGICAL           :: ENDSEC ! T=end of section detected
    CHARACTER(LENREC) :: RECORD ! Record read from Driver Table
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! First non-blank record should be text for Header
  CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
  IF ( ENDSEC ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVHDR: No header text found in *HDR section'
    RETURN
  END IF
!
  TXTHDR = '!' // RECORD(1:79)
  CALL WRTLOG ( '  ' // TXTHDR )
!
! Shouldn't be any more non-blank records in *HDR section
  CALL ENDCHK ( LUNDRV, '*HDR', FAIL, ERRMSG )
!
END SUBROUTINE DRVHDR
END MODULE DRVHDR_SUB
