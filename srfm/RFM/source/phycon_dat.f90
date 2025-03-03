MODULE PHYCON_DAT
!
! VERSION
!   17JAN24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Mathematical and physical constants
!   See also PHYADJ_DAT for 'adjustable' constants.
!
! REFERENCE
!   NIST: http://physics.nist.gov/cuu/Units/index.html
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
!
! GLOBAL CONSTANTS
    REAL(R4), PARAMETER :: ATMB = 1013.25        ! [mb/atm]
    REAL(R4), PARAMETER :: AVOG = 6.02214199E26  ! Avogadro constant [kmole-1]
    REAL(R4), PARAMETER :: RGAS = 8314.472       ! Gas constant [J/K/kmole] 
    REAL(R4), PARAMETER :: PI = 3.141592654      ! pi
    REAL(R4), PARAMETER :: PREREF = 1.0          ! Line data ref. pres. [atm] 
    REAL(R4), PARAMETER :: TEMREF = 296.0        ! Line data reference temp [K]
    REAL(R4), PARAMETER :: DG2RAD = PI/180.0     ! deg to radians conv.factor
!
    REAL(R8), PARAMETER :: C1 = 1.191042723D-3   ! 1st r.c [nW/(cm2.sr.(cm-1)4)]
    REAL(R8), PARAMETER :: C2 = 1.4387752D0      ! 2nd rad.con.=hc/k [K/cm-1]
    REAL(R8), PARAMETER :: VLIGHT = 2.99792458D8 ! Speed of light [m/s]
    REAL(R8), PARAMETER :: GHZ2CM = 1.0D7/VLIGHT ! GHz to cm-1 conv.factor
!
END MODULE PHYCON_DAT
