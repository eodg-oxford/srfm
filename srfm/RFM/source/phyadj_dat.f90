MODULE PHYADJ_DAT
!
! VERSION
!   12MAY24 AD Checked.
!   12MAR19 AD Add WNORFR. Checked.
!   01MAY17 AD F90 conversion of phycom.inc. Checked.
!
! DESCRIPTION
!   Adjustable physical constants used within the RFM
!   Modified by values read from *PHY section of driver table
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL VARIABLES
    REAL(R4) :: CPKMOL = 29012.0    ! Molar heat cap of air [J/K/kmole]
    REAL(R4) :: GRAVTY = 9.80665    ! Accn due to gravity [m/s^2]
    REAL(R4) :: TEMSPA = 2.7        ! Cosmic background temperature [K]
    REAL(R4) :: WGTAIR = 28.964     ! Molecular weight of air [kg/kmole]
!
    REAL(R8) :: RADCRV = 6367.421D0 ! Radius of curvature [km]
    REAL(R8) :: WNORFR = 1000.0D0   ! Wavenumber [cm-1] for refractivity calc.
!
END MODULE PHYADJ_DAT
