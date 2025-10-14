MODULE VVWCOR_SUB
CONTAINS
SUBROUTINE VVWCOR ( DWNO, K )
!
! VERSION
!   13MAY24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Apply Van Vleck-Weisskopf correction to line shape
!   Called by LINSHP if VVW flag set true.
!   Not called if VVW lineshape already used for calculation.
!   Uses path-adjusted line data in ADJCOM
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL VARIABLES
    USE ADJCOM_DAT ! Path-adjusted line data
    USE PHYCON_DAT, ONLY: PI
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(IN)    :: DWNO(:) ! Array of Wavenumbers [/cm]
    REAL(R4), INTENT(INOUT) :: K(:)    ! Absorption 
!
! LOCAL VARIABLES
    REAL(R4) :: TOP  ! Numerator of Lorentz expression
    REAL(R4) :: W2   ! Width^2
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  TOP = STRADJ * WIDADJ / PI                   
  W2 = WIDADJ**2
!
  K = SNGL ( ( DWNO / WNOADJ )**2 ) * &
             ( K +  TOP / ( SNGL ( ( DWNO + WNOADJ )**2 ) + W2 ) )
!
END SUBROUTINE VVWCOR
END MODULE VVWCOR_SUB
