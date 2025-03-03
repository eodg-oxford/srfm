MODULE TANUNI_SUB
CONTAINS
SUBROUTINE TANUNI ( PARAM, VALUE, FAIL, ERRMSG )
!
! VERSION
!   26MAY24 AD Checked.
!   25MAR19 AD Original. Checked.
!
! DESCRIPTION
!   Check UNITS specified in *TAN/*LEN section
!   Called by DRVTAN if PARAM=VALUE field found *TAN/*LEN section
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE TANCOM_DAT, ONLY: UNITAN, USRUNI ! Units for  path length
    USE FLGCOM_DAT, ONLY: HOMFLG ! T = use homogeneous path
!
! SUBROUTINES
    USE LOCASE_FNC ! Convert text string to lower case
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: PARAM  ! PARAM from PARAM=VALUE pair
    CHARACTER(*),  INTENT(IN)  :: VALUE  ! VALUE from PARAM=VALUE pair
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL :: GOTUNI = .FALSE. ! T=UNITS already set
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  UNITAN = 0.0   ! New value will be set unless fatal error occurs
!
  IF ( PARAM .EQ. 'UNITS' ) THEN
    IF ( GOTUNI ) THEN
      ERRMSG = 'F-TANUNI: Repeated setting of UNITS in *TAN/*LEN section'
    ELSE IF ( HOMFLG ) THEN
      GOTUNI = .TRUE.
      USRUNI = LOCASE ( VALUE ) 
      SELECT CASE ( USRUNI )
      CASE ( 'km' ) ; UNITAN = 1.0
      CASE ( 'm'  ) ; UNITAN = 1.0E-3
      CASE ( 'cm' ) ; UNITAN = 1.0E-5
      CASE ( 'mm' ) ; UNITAN = 1.0E-6
      CASE DEFAULT
        ERRMSG = 'F-TANUNI: UNITS= followed by unrecognised length units: ' &
                 // VALUE
      END SELECT
    ELSE
      ERRMSG = 'F-TANUNI: length units can only be set with HOM flag'
    END IF
  ELSE
    ERRMSG = 'F-TANUNI: Only UNITS=(value) allowed in *TAN/*LEN section' 
  END IF
!
  FAIL = UNITAN .EQ. 0.0
  IF ( .NOT. FAIL ) &
    CALL WRTLOG ( 'I-TANUNI: Defining *TAN/*LEN units as: ' // VALUE ) 
!
END SUBROUTINE TANUNI
END MODULE TANUNI_SUB

