MODULE COOWGT_SUB
CONTAINS
SUBROUTINE COOWGT 
!
! VERSION
!   31JUL24 AD Checked.
!   01MAY17 AD F90 conversion of part of rfmflx.for. Checked.
!
! DESCRIPTION
!   Calculate cooling rate weights
!   Called once by DRVCHK if COO flag enabled.
!
!   Cooling rates are calculated from dT/dt = dF/dz*(1/rho.cp)
!   where rho = atmospheric density
!   dF/dz is calculated from a quadratic fit to F(z) at 3 successive levels
!   in the atmospheric profile.
!   F is in nW/cm2..., z in km, rho in molec/cm3, cp in J/K/kmol=W.s/K/kmol, 
!   so dT/dt*(1e-9W/nW)*(1e-5km/cm)*(Na molec/kmol)*(864000s/day)= K/day
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE PHYCON_DAT, ONLY: AVOG   ! Avogadro's constant [kmole-1]
    USE PHYADJ_DAT, ONLY: CPKMOL ! Molar heat cap of air [J/K/kmole]
!
  IMPLICIT NONE
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: MAXWGT = 5 !  Max.no o/p levels affected by atm#
!                                           Quadratic fit implies a maximum=5
! LOCAL VARIABLES
    INTEGER(I4) :: IATM        ! Atmospheric level/layer counter
    INTEGER(I4) :: ITAN        ! Output level index (1:NTAN)
    INTEGER(I4) :: JATM        ! Atmospheric level/layer counter
    INTEGER(I4) :: KATM        ! Atmospheric level/layer counter
    INTEGER(I4) :: NWTI,NWTJ,NWTK ! No.weights for atm levels i,j,katm
    REAL(R4)    :: COOFAC      ! Factor for cooling rates
    REAL(R4)    :: XIJ,XJK,XKI ! Altitude separation [km] between atm.levels
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  ALLOCATE ( NCOATM(NATM), ICOATM(MAXWGT,NATM), WCOATM(MAXWGT,NATM) ) 
  NCOATM = 0
!
  DO JATM = 1, NATM
    ITAN = ITNATM(JATM) 
    IF ( ITAN .EQ. 0 ) CYCLE             ! not an output level
! For cooling rates, fit quadratic to atm.levels i,j,k and calculate gradient
! at level j. Coefficients do not depend on sequence i,j,k so can also be used
! if j is one of the end points (if JATM=1 or NATM).
    IATM = JATM-1                        ! sequence i,j,k
    KATM = JATM+1
    IF ( IATM .LT. 1 ) IATM = 3          ! sequence j,k,i 
    IF ( KATM .GT. NATM ) KATM = NATM-2  ! sequence k,i,j
    NWTI = NCOATM(IATM) + 1
    NWTJ = NCOATM(JATM) + 1
    NWTK = NCOATM(KATM) + 1
    IF ( NWTI .GT. MAXWGT .OR. NWTJ .GT. MAXWGT .OR. NWTK .GT. MAXWGT ) &
      STOP 'F-COOWGT: No. Weights > expected max array size'
    ICOATM(NWTI,IATM) = ITAN
    ICOATM(NWTJ,JATM) = ITAN
    ICOATM(NWTK,KATM) = ITAN
    XIJ = HGTATM(IATM) - HGTATM(JATM)
    XJK = HGTATM(JATM) - HGTATM(KATM)
    XKI = HGTATM(KATM) - HGTATM(IATM)
    COOFAC = (AVOG * 86400.0E-14) / DNSATM(JATM) / CPKMOL
    WCOATM(NWTI,IATM) = - COOFAC * XJK / XIJ / XKI
    WCOATM(NWTJ,JATM) = COOFAC * ( 1.0/XJK - 1.0/XIJ ) 
    WCOATM(NWTK,KATM) = COOFAC * XIJ / XJK / XKI
    NCOATM(IATM) = NWTI
    NCOATM(JATM) = NWTJ
    NCOATM(KATM) = NWTK
  END DO  
!
END SUBROUTINE COOWGT
END MODULE COOWGT_SUB
