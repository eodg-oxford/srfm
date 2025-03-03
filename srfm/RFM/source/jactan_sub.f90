MODULE JACTAN_SUB
CONTAINS
SUBROUTINE JACTAN ( JDX, TARGET, NFLD, FAIL, ERRMSG )
!
! VERSION
!   30NOV23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Pass list of tangent heights to Jacobian check
!   Called by DRVJAC if 'TAN' in *JAC section list of tangent heights
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE TANCOM_DAT ! Tangent path data
!
! SUBROUTINES
    USE JACALT_SUB ! Set Jacobian perturbation altitudes
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: JDX    ! Index of Jacobian target species
    CHARACTER(*),  INTENT(IN)  :: TARGET ! Target species for Jacobian
    INTEGER(I4),   INTENT(IN)  :: NFLD   ! Field# within record of 'TAN'
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: ITAN ! Counter for tangent heights
!
! EXECUTABLE CODE ------------------------------------------------------------
!
  FAIL = .TRUE.
! TAN only valid if limb-viewing mode
  IF ( .NOT. LIMTAN ) THEN
    ERRMSG = 'F-JACTAN: ''TAN'' entry in *JAC section ' & 
             // 'only valid for limb-viewing mode'
! and with tangent paths expressed as refracted (ie actual) atmos.altitudes
! ie the same coordinates as used in the normal *JAC section
  ELSE IF ( USRELE .OR. USRGEO ) THEN
    ERRMSG = 'F-JACTAN: ''TAN'' entry in *JAC section ' & 
             // 'only valid with *TAN section (refrac.alts)'
  ELSE IF ( NFLD .NE. 2 ) THEN
    ERRMSG = 'F-JACTAN: ''TAN'' entry must be 3rd field in record'
  ELSE
    FAIL = .FALSE.
  END IF
  IF ( FAIL ) RETURN
!
  DO ITAN = 1, NTAN
    CALL JACALT ( JDX, TARGET, .FALSE., TAN(ITAN)%USR, FAIL, ERRMSG, ITAN=ITAN )
    IF ( FAIL ) RETURN
    CALL WRTLOG ( TAN(ITAN)%STR, .TRUE. )
  END DO
! 
END SUBROUTINE JACTAN
END MODULE JACTAN_SUB
