MODULE COOWGT_FNC
CONTAINS
FUNCTION COOWGT ( HGT, DNS )
!
! VERSION
!   07JUL26 AD Checked
!   01AUG25 AD Rewritten, simplified and converted to function 
!   31JUL24 AD Checked.
!   01MAY17 AD F90 conversion of part of rfmflx.for. Checked.
!
! DESCRIPTION
!   Calculate cooling rate weights
!   Called once by COOLAY if COO flag enabled.
!
!   Cooling rates are calculated from dT/dt = dF/dz*(1/rho.cp)
!   where rho = atmospheric density
!   dF/dz is calculated from a quadratic fit to F(z) at 3 successive levels
!   in the atmospheric profile.
!   F is in nW/cm2..., z in km, rho in molec/cm3, cp in J/K/kmol=W.s/K/kmol, 
!   so dT/dt*(1e-9W/nW)*(1e-5km/cm)*(Na molec/kmol)*(864000s/day)= K/day

! For cooling rates, fit quadratic to atm.levels i,j,k and calculate gradient
! at level j. Coefficients do not depend on sequence i,j,k so can also be used
! for heights in any order
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE TANCOM_DAT ! Tangent path data
    USE PHYCON_DAT, ONLY: AVOG   ! Avogadro's constant [kmole-1]
    USE PHYADJ_DAT, ONLY: CPKMOL ! Molar heat cap of air [J/K/kmole]
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: HGT(:)  ! Profile Level altitudes [km]
    REAL(R4), INTENT(IN) :: DNS     ! Central level density [molec/cm3]
!
! FUNCTION TYPE
    REAL(R4) :: COOWGT(3) 
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: MAXWGT = 5 !  Max.no o/p levels affected by atm#
!                                           Quadratic fit implies a maximum=5
! LOCAL VARIABLES
    REAL(R4)    :: COOFAC      ! Factor for cooling rates
    REAL(R4)    :: XIJ,XJK,XKI ! Altitude separation [km] between atm.levels
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  XIJ = HGT(1) - HGT(2)
  XJK = HGT(2) - HGT(3)
  XKI = HGT(3) - HGT(1)
  COOFAC = (AVOG * 86400.0E-14) / DNS / CPKMOL
  COOWGT(1) = -COOFAC * XJK / XIJ / XKI
  COOWGT(2) =  COOFAC * ( 1.0/XJK - 1.0/XIJ ) 
  COOWGT(3) =  COOFAC * XIJ / XJK / XKI
!
END FUNCTION COOWGT
END MODULE COOWGT_FNC
