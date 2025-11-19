MODULE LORSHP_SUB
CONTAINS
SUBROUTINE LORSHP ( DWNO, K )
!
! VERSION
!   07MAY24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Calculate Lorentz Line shape
!   Called by LINSHP if Lorentz lineshape selected.
!   Uses path-adjusted line data in ADJCOM
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
    REAL(R8), INTENT(IN)  :: DWNO(:) ! Array of Wavenumbers [cm-1]
    REAL(R4), INTENT(OUT) :: K(:)    ! Absorption 
!
! LOCAL VARIABLES
    REAL(R4) :: TOP ! Numerator of Lorentz expression
    REAL(R4) :: W2  ! Width^2
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  TOP = STRADJ * WIDADJ / PI                   
  W2 = WIDADJ**2
!
  K = TOP / ( SNGL ( ( DWNO - WNOADJ )**2 ) + W2 )
!
END SUBROUTINE LORSHP
END MODULE LORSHP_SUB                                                 
