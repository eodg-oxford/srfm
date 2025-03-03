MODULE SUB25W_SUB
CONTAINS
SUBROUTINE SUB25W ( ABSORP )
!
! VERSION
!   09MAY24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Subtract absorption value at 25cm-1 from line centre 
!   Called by LINSHP for H2O lines if CTM flag enabled
!   Required for use with the CKD H2O continuum coefficients.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ADJCOM_DAT ! Path-adjusted line data
    USE PHYCON_DAT, ONLY: PI
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(INOUT) :: ABSORP(:) ! Absorption coefficient
!
! LOCAL VARIABLES
    REAL(R4) :: ABS25W ! Absorption evaluated at 25cm-1 from line centre
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! ABS25W is just the Lorentzian line shape evaluated at 25cm-1 from centre
  ABS25W = STRADJ * WIDADJ / PI / ( 625.0 + WIDADJ**2 )
!
! Assume CKD continuum allows for line contributions to absorption .LE. 0 
  ABSORP = ABSORP - ABS25W
!
END SUBROUTINE SUB25W
END MODULE SUB25W_SUB
