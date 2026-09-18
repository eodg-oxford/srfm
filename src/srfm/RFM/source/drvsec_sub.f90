MODULE DRVSEC_SUB
CONTAINS
SUBROUTINE DRVSEC ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   13JUL26 AD Checked.
!   01AUG25 AD Original. Extracted from DRVTAN.
!
! DESCRIPTION
!   Read RFM driver table *SEC section
!   Called by RFMDRV once if paths through plane-parallel atmosphere selected.
!   Reads inputs for sec zenith angles.
!
! VARIABLE KINDS
    USE KIND_DAT
! 
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE TANCOM_DAT ! Tangent path data
    USE FLGCOM_DAT, ONLY: NADFLG !  T = nadir-viewing
!
! SUBROUTINES
    USE NXTFFL_SUB ! Load next field from rfm.drv, expanding filenames
    USE TANFLD_SUB ! Check if string is valid *TAN entry and insert in TANCOM
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
    INTEGER(I4)       :: LENGTH ! Length of field read from driver table
    REAL(R4)          :: VALUE  ! Value read from FIELD
    CHARACTER(LENREC) :: FIELD  ! String containing sec theta value
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  CALL WRTLOG ( 'I-DRVSEC: ', .TRUE. ) 
!
! Read each field in *SEC section
  DO
    CALL NXTFFL ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
    CALL TANFLD ( FIELD, VALUE, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    IF ( VALUE .LT. 1.0 ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-DRVSEC: Specified NADir/ZENith view sec(theta) < 1, value =' &
               // TRIM ( FIELD(1:20) ) 
      RETURN
    END IF
    CALL WRTLOG ( ' '//FIELD, .TRUE. ) 
  END DO
!
! Check at least one tangent height supplied
  IF ( NTAN .EQ. 0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVSEC: No entries in *SEC section'
  END IF
!
  TAN%ITN = 1 
  TAN%IAT = 1
  TAN%SFC = NADFLG
  TAN%SEC = TAN%USR
  TAN(1)%CLC = .TRUE.
!
END SUBROUTINE DRVSEC
END MODULE DRVSEC_SUB

