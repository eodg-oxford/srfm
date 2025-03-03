MODULE REXPTH_SUB
CONTAINS
SUBROUTINE REXPTH 
!
! VERSION
!   22SEP24 AD Checked.
!   01MAY17 AD F90 original. Checked.
! 
! DESCRIPTION    
!   Set paths for Rayleigh Extinction
!   Called by RFMPTH if REX flag enabled.
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
    USE ATMCOM_DAT ! Atmospheric profile data
    USE CLCCOM_DAT ! Calculated path segments
    USE IDXCON_DAT ! RFM/HITRAN indices for specific molecules
    USE REXCOM_DAT ! Rayleigh extinction path data
!
! SUBROUTINES
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
    USE VAL1DI_GEN ! Interpolate value from 1D array
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: ICLC  ! Counter for calculated path segments
    INTEGER(I4) :: IGSAIR ! Index of air in GAS arrays
    INTEGER(I4) :: IREX  ! Counter for REX calc path segments
    INTEGER(I4) :: IVCO2  ! Index of CO2 in VMRATM (and GAS) arrays
    INTEGER(I4) :: IVN2  ! Index of N2 in VMRATM (and GAS) arrays
    INTEGER(I4) :: IVO2  ! Index of O2 in VMRATM (and GAS) arrays
    REAL(R4)    :: LNP ! Log(pressure) [p in mb]
    REAL(R4)    :: VMRA   = 0.009 * 1.0E6  ! VMR [ppmv] of Argon
    REAL(R4)    :: VMRCO2 = 0.0004 * 1.0E6 ! VMR [ppmv] of CO2
    REAL(R4)    :: VMRN2  = 0.78 * 1.0E6   ! VMR [ppmv] of N2
    REAL(R4)    :: VMRO2  = 0.21 * 1.0E6   ! VMR [ppmv] of O2
    REAL(R4)    :: VMRSUM ! Sum of contributing VMR values
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Routine should only be called if REXFLG enabled so 'air' should be loaded
  IGSAIR = IDXGAS ( IDXAIR )
  IF ( IGSAIR .EQ. 0 ) STOP 'F-REXPTH: Logical error'
!
! The following coud return 0 if these molecules are not used in RFM calc.
  IVN2  = IDXGAS ( IDXN2 )
  IVO2  = IDXGAS ( IDXO2 )
  IVCO2 = IDXGAS ( IDXCO2 )
!
! No. of Calc paths for 'air' requiring Rayleigh Extinction calculation
  NREX = COUNT ( CLC%IGS .EQ. IGSAIR ) 
  ALLOCATE ( REX(NREX) ) 
!
  IREX = 0
  DO ICLC = 1, NCLC
    IF ( CLC(ICLC)%IGS .NE. IGSAIR ) CYCLE
    IREX = IREX + 1
    REX(IREX)%ICL = ICLC
    REX(IREX)%AMT = CLC(ICLC)%AMT
    LNP = LOG ( CLC(ICLC)%PRE )
    IF ( IVN2 .NE. 0 ) &
      VMRN2 = VAL1DI ( LNPATM, LNP, VMRATM(:,IVN2), .NOT. LINVMR(IVN2) ) 
    IF ( IVO2 .NE. 0 ) &
      VMRO2 = VAL1DI ( LNPATM, LNP, VMRATM(:,IVO2), .NOT. LINVMR(IVO2) ) 
    IF ( IVCO2 .NE. 0 ) &
      VMRCO2 = VAL1DI ( LNPATM, LNP, VMRATM(:,IVCO2), .NOT. LINVMR(IVCO2) ) 
    VMRSUM = VMRN2 + VMRO2 + VMRCO2 + VMRA
    REX(IREX)%A0 = ( 1.034 * VMRN2 + 1.096 * VMRO2 + &
                  1.15 * VMRCO2 +  1.0 * VMRA      ) / VMRSUM
    REX(IREX)%A2 = ( 3.17E-12 * VMRN2 + 1.385E-11 * VMRO2 ) / VMRSUM
    REX(IREX)%A4 = 1.448E-20 * VMRO2 / VMRSUM
  END DO
!
END SUBROUTINE REXPTH
END MODULE REXPTH_SUB
