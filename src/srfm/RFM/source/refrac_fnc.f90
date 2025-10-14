MODULE REFRAC_FNC
CONTAINS
FUNCTION REFRAC ( PRE, TEM, WNO )
!
! VERSION
!   31JAN24 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Calculate refractivity
!   Called by ATMAUX
!
! REFERENCE
!   Tables of Physical and Chemical Constants (16th Ed), 
!     G. W. C. Kaye and T. H. Laby
!     Published by Longman, 1995.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE PHYCON_DAT, ONLY: ATMB ! [mb/atm]
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: PRE(:) ! Pressure [hPa]
    REAL(R4), INTENT(IN) :: TEM(:) ! Temperature [K]
    REAL(R8), INTENT(IN) :: WNO    ! Wavenumber [cm-1]
!
! FUNCTION TYPE
    REAL(R4) :: REFRAC(SIZE(PRE))  ! function returns array size matching pre
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: TREF = 288.16 ! Ref.Temp. for Edlen refrac. calc.
!
! LOCAL VARIABLES
    REAL(R4) :: SIGMA2 ! (1/Wavelength[/micron])^2
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! K&L formula expressed in sigma^2, where sigma=1/wavelength[microns] 
  SIGMA2 = SNGL(WNO*0.0001)**2
!
  REFRAC = 1.0E-8 * PRE / ATMB * TREF / TEM * &
   ( 8342.54 + 2406147.0 / ( 130.0 - SIGMA2 ) + 15998.0 / (  38.9 - SIGMA2 ) ) 
!
END FUNCTION REFRAC
END MODULE REFRAC_FNC
