MODULE RADSFC_SUB
CONTAINS
SUBROUTINE RADSFC ( OPT, RAD, JDX )
!
! VERSION
!   08JAN24 AD Checked.
!   01JUL17 AD F90 Original. Checked.
!
! DESCRIPTION
!   Surface contribution to radiative transfer calculation
!   Called by SPCRAD for any surface-intersecting ray paths.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FINCOM_DAT ! Finemesh data
    USE SFCCOM_DAT ! Surface parameters
    USE JACCOM_DAT, ONLY: JDXSFE, JDXSFT ! Indices of surface perturbations
    USE PTBCON_DAT, ONLY: PTBSFE, PTBSFT ! Jacobian perturbation sizes
!
! SUBROUTINES
    USE INTERP_GEN ! Interpolate array
    USE PLANCK_FNC ! Planck Function
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(INOUT) :: OPT(:) ! Absorption along path
    REAL(R8), INTENT(INOUT) :: RAD(:) ! Radiance along path
    INTEGER(I4), INTENT(IN) :: JDX    ! Index of Jacobian
!
! LOCAL VARIABLES
    REAL(R8) :: BBF(NFIN) ! Surface Planck emission
    REAL(R8) :: EMS(NFIN) ! Surface emissivity [0:1]
!
! EXECUTABLE CODE ------------------------------------------------------------
!
  IF ( JDX .EQ. JDXSFT ) THEN
    BBF = PLANCK ( TEMSFC + PTBSFT, WNOFIN )
  ELSE
    BBF = PLANCK ( TEMSFC, WNOFIN ) 
  END IF
!
  IF ( RFLSFC ) THEN
    EMS = INTERP ( WNOSFC, WNOFIN, EMSSFC )
    IF ( JDX .EQ. JDXSFE ) EMS = EMS + PTBSFE
    RAD = RAD * ( 1.0D0 - EMS ) + BBF * EMS
    OPT = OPT - LOG ( MAX ( TINY(1.0D0), (1.0D0-EMS) ) ) 
  ELSE
    RAD = BBF
    OPT = 0.0D0
  END IF
!
END SUBROUTINE RADSFC
END MODULE RADSFC_SUB

