MODULE DRVLEN_SUB
CONTAINS
SUBROUTINE DRVLEN ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   11JUL26 AD Checked.
!   01AUG25 AD Extracted from DRVTAN.
!
! DESCRIPTION
!   Read RFM driver table *LEN section
!   Called by RFMDRV once if HOM flag is enabled.
!   Reads inputs for homogeneous path length and (optionally) units.
!
! VARIABLE KINDS
    USE KIND_DAT
! 
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE TANCOM_DAT ! Tangent path data
    USE FLGCOM_DAT, ONLY: SFCFLG ! T = Allow for opaque surface 
!
! SUBROUTINES
    USE LOCASE_FNC ! Convert text string to lower case
    USE NXTFFL_SUB ! Load next field from rfm.drv, expanding filenames
    USE PARFLD_SUB ! Extract Parameter=Value string from record
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
    LOGICAL           :: GOTPAR ! T=FIELD contains PARAM=VALUE pair
    LOGICAL           :: GOTUNI = .FALSE. ! T=UNITS already set
    INTEGER(I4)       :: LENGTH ! Length of field read from driver table
    REAL(R4)          :: LENVAL ! Length value read from FIELD
    CHARACTER(LENREC) :: FIELD  ! Tangent height or Tan.Hgt filename
    CHARACTER(LENREC) :: PARAM  ! PARAM part of PARAM=VALUE pair
    CHARACTER(LENREC) :: VALUE  ! VALUE part of PARAM=VALUE pair
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  CALL WRTLOG ( 'I-DRVLEN: ', .TRUE. ) 
!
! Read each field in *LEN section
  DO
    CALL NXTFFL ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
    CALL PARFLD ( FIELD, GOTPAR, PARAM, VALUE ) 
    IF ( GOTPAR ) THEN
      IF ( PARAM .EQ. 'UNITS' ) THEN
        IF ( GOTUNI ) THEN
          ERRMSG = 'F-DRVLEN: Repeated setting of UNITS in *TAN/*LEN section'
          FAIL = .TRUE.
          RETURN
        END IF
        GOTUNI = .TRUE.
        USRUNI = LOCASE ( VALUE ) 
        SELECT CASE ( USRUNI )
        CASE ( 'km' ) ; UNITAN = 1.0
        CASE ( 'm'  ) ; UNITAN = 1.0E-3
        CASE ( 'cm' ) ; UNITAN = 1.0E-5
        CASE ( 'mm' ) ; UNITAN = 1.0E-6
        CASE DEFAULT
          ERRMSG = 'F-DRVLEN: UNITS= followed by unrecognised length units: ' &
                   // VALUE(1:20)
        END SELECT
      ELSE
        ERRMSG = 'F-DRVLEN: Only UNITS=(value) allowed in *TAN/*LEN section'
        FAIL = .TRUE.
        RETURN
      END IF
    ELSE   
      CALL TANFLD ( FIELD, LENVAL, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      IF ( LENVAL .LE. 0.0 ) THEN 
        FAIL = .TRUE.
        ERRMSG = 'F-DRVLEN: Specified Homog. Path Length .LE. 0, value =' &
                 // TRIM ( FIELD(1:20) )
      END IF
    END IF
    IF ( FAIL ) RETURN
    CALL WRTLOG ( ' '//FIELD, .TRUE. ) 
  END DO
!
! Check at least one tangent height supplied
  IF ( NTAN .EQ. 0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVLEN: No entries in *LEN section'
  END IF
!
  TAN%ITN = 1
  TAN%IAT = 1
  TAN%SFC = SFCFLG
  TAN%SEC = TAN%USR * UNITAN
  TAN(1)%CLC = .TRUE.
!
END SUBROUTINE DRVLEN
END MODULE DRVLEN_SUB

