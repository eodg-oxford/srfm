MODULE DOPSHP_SUB
CONTAINS
SUBROUTINE DOPSHP ( DWNO, K )
!
! VERSION
!   06MAY24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Calculate Doppler Line shape
!   Called by LINSHP if Doppler lineshape selected.
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
    REAL(R8), INTENT(IN)  :: DWNO(:) ! Array of Wavenumbers [/cm]
    REAL(R4), INTENT(OUT) :: K(:)    ! Absorption 
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: LN2 = LOG ( 2.0 ) 
    REAL(R4), PARAMETER :: SL2RPI = SQRT ( LN2 / PI ) 
!
! LOCAL VARIABLES
    REAL(R4) :: ARG ! exponential term multiplying (wn-wn0)^2
    REAL(R4) :: COF ! Scaling factor for exponential.
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  ARG = - LN2 / DOPADJ**2           
  COF = SL2RPI * STRADJ / DOPADJ      ! SQRT(ln 2/PI) = 0.46971864
!
  K = COF * EXP ( ARG * SNGL ( ( DWNO - WNOADJ )**2 ) )
!
END SUBROUTINE DOPSHP
END MODULE DOPSHP_SUB
