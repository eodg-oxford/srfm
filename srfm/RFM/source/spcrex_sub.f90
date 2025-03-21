MODULE SPCREX_SUB
CONTAINS
SUBROUTINE SPCREX 
!
! VERSION
!   08DEC23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   Calculate the contribution of Rayleigh extinction
!   Called by RFMSPC for each wide mesh interval if REX flag enabled.
!
!   Eq(2) of Bodhaine et al gives scattering cross-section per molecule
!     sigma = [24 pi^3 (n^2-1)^2]/[lambda^4 rho^2 (n^2+2)^2] * K
!   where n is refractive index, lambda = wavelength, rho=number density,
!   K is the King Factor, or depolarization term.
!   Note that (n^2-1)^2 is effectively units of molecules
!   Writing n = N+1, where N is the refractivity <<1, this is approximately
!     sigma = [24 pi^3 (2N)^2]/[lambda^4 rho^2 3^2] * K
!           = [32 pi^3/3] ( N^2/rho^2 ) (1/lambda^4) * K
!
!   Refractivity depends primarily on density, so N/rho is virtually constant
!   with altitude, so it should be OK to use value at lower boundary of layer
!   rather than anything more detailed.
!
!   From function REFRAC 
!     Refractivity, RFR = 1E-8 * PRE / ATMB * TREF / TEM * spectral_factor
!   From subroutine ATMAUX
!     Density, DNS = PRE * AVOG * 1E-4 / RGAS / TEM
!   So RFR/DNS = 1E-4 * TREF * RGAS / ATMB / AVOG * spectral_factor
!              = 1E-4 * 288 * 8314 / 1013 / 6E26 = 3.93E-28 * spectral_factor
!
!   Bodhaine expression N/rho is implicitly for 1 molecule, so square has 
!   dimension [cm^6/molec] rather than [cm^6/molec^2]
!   Similarly for King factor, which depends on air mixture therefore varies
!   only slowly with altitude
!
! REFERENCES
!   Bodhaine, B. A., N. B. Wood, E. G. Dutton, and J. R. Slusser,
!     On Rayleigh Optical Depth Calculations, J. Atm. Ocean. Tech.,
!     16(11), 1854--1861, 1999.
!               
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FINCOM_DAT ! Finemesh data
    USE REXCOM_DAT ! Rayleigh extinction path data
    USE PHYCON_DAT, ONLY: AVOG, PI ! Avogadro's number, pi
!
  IMPLICIT NONE
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: PIFACT = 32.0 * PI**3 / 3.0 ! numerical factor
    REAL(R4), PARAMETER :: RFRDEN = 3.93E-28           ! [cm3/molec] Refrac/Dens
    REAL(R4), PARAMETER :: SIGMA = PIFACT * ( AVOG * RFRDEN ) * RFRDEN 
!
! LOCAL VARIABLES
    INTEGER(I4) :: ICLC         ! Index of corresponding calculated path
    INTEGER(I4) :: IREX         ! Counter for Rayleigh Extinction paths
    REAL(R4)    :: KINGDF(NFIN) ! King factor for depolarization
    REAL(R4)    :: RWNO2(NFIN)  ! Wavenumber^2 [cm-2]
    REAL(R4)    :: RWNO4(NFIN)  ! Wavenumber^4 [cm-4]
    REAL(R4)    :: WPM2(NFIN)   ! (Wavenumber [/micron])^2
    REAL(R4)    :: SPECF(NFIN)  ! Spectral factor for refractivity
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  RWNO2 = REAL ( WNOFIN**2 ) 
  RWNO4 = RWNO2**2
  WPM2 = RWNO2 * 1.0E-8    ! Conversion from /cm to /micron, squared.
  SPECF = ( 8342.54 + 2406147.0 / ( 130.0 - WPM2 ) + &
                      15998.0 / (  38.9 - WPM2 )      )**2
  DO IREX = 1, NREX
    ICLC = REX(IREX)%ICL
    KINGDF = REX(IREX)%A0 + RWNO2 * REX(IREX)%A2 + RWNO4 * REX(IREX)%A4
    ABSFIN(:,ICLC) = SIGMA * RWNO4 * SPECF * KINGDF * REX(IREX)%AMT 
  END DO
!
END SUBROUTINE SPCREX
END MODULE SPCREX_SUB


