MODULE CHKTAN_SUB
CONTAINS
SUBROUTINE CHKTAN ( FAIL, ERRMSG )
!
! VERSION
!   26FEB23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Check Limb-viewing tangent paths
!   Called by DRVCHK after all rfm.drv contents read.
!   Requires inputs from both *TAN and *OBS sections.
!
! VARIABLE KINDS
    USE KIND_DAT
! 
! GLOBAL DATA
    USE TANCOM_DAT ! Tangent heights
!
! SUBROUTINES
    USE CHKLIM_SUB ! Check Limb-viewing tangent paths
!
  IMPLICIT NONE 
!
! ARGUMENTS
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: ITAN ! Counter for tangent paths
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO ITAN = 1, NTAN
    CALL CHKLIM ( ITAN, USRELE, USRGEO, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
  END DO
!
END SUBROUTINE CHKTAN
END MODULE CHKTAN_SUB

