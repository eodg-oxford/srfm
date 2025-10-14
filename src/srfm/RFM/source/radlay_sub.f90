MODULE RADLAY_SUB
CONTAINS
SUBROUTINE RADLAY ( ITAN, IATM, IDIR, RAD, OPT )
!
! VERSION
!   30APR24 AD Checked.
!   05MAR19 AD Remove JTAN argument. Scale path AMT by TAN(ITAN)%SEC
!              Assume IDXPTH will always assign path. Checked.
!   01JUL17 AD F90 original. Checked.
!
! DESCRIPTION
!   Radiative transfer calculation through a layer
!   Called by SPCRAD
!   Integrates from the remote side of the atmosphere to the observer, so that
!   input radiance is reduced by layer transmittance before adding on the layer
!   contribution.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FINCOM_DAT ! Finemesh data
    USE FLGCOM_DAT ! Option flags
    USE PTHCOM_DAT ! Path segment data
    USE TANCOM_DAT ! Tangent path data
    USE ATMCOM_DAT, ONLY: NVMR   ! No.different absorbers
!
! SUBROUTINES
    USE IDXPTH_FNC ! Index in PTHCOM of tan/atm/gas/dir
    USE PLANCK_FNC ! Planck function
    USE SRCBFX_FNC ! Modified source function for BFX flag
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN)    :: ITAN   ! Index of tangent path
    INTEGER(I4), INTENT(IN)    :: IATM   ! Index of atmos.layer
    INTEGER(I4), INTENT(IN)    :: IDIR   ! +1 = upward path, -1 = downward path
    REAL(R8),    INTENT(INOUT) :: RAD(:) ! Radiance along path
    REAL(R8),    INTENT(OUT)   :: OPT(:) ! Layer optical depth
!
! LOCAL VARIABLES
    INTEGER(I4) :: ICLC      ! Index of calculated path 
    INTEGER(I4) :: IPTH      ! Path counter
    INTEGER(I4) :: IVMR      ! Absorber counter
    REAL(R4)    :: AMT       ! Absorber amount in path [kmol/cm^2]
    REAL(R4)    :: TEM       ! Path CG temperature [K]
    REAL(R8)    :: SRC(NFIN) ! Absorption-wgtd source fn. for layer
    REAL(R8)    :: TRA(NFIN) ! Transmittance of layer
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  SRC = 0.0D0
  OPT = 0.0D0
  DO IVMR = 1, NVMR
    IPTH = IDXPTH ( ITAN, IATM, IVMR, IDIR )
! Possibly a Jacobian tan ray without an assigned perturbed path
!    IF ( IPTH .EQ. 0 ) IPTH = IDXPTH ( TAN(ITAN)%ITN, IATM, IVMR, IDIR )  
    ICLC = PTH(IPTH)%ICL
    AMT  = PTH(IPTH)%AMT * TAN(ITAN)%SEC
    TEM  = PTH(IPTH)%TEM
    IF ( PTH(IPTH)%NTE ) THEN
      SRC = SRC + DBLE ( CNTFIN(:,ICLC) * AMT ) * PLANCK ( TEM, WNOFIN ) 
    ELSE
      SRC = SRC + DBLE ( ABSFIN(:,ICLC) * AMT ) * PLANCK ( TEM, WNOFIN ) 
    END IF
    OPT = OPT + DBLE ( ABSFIN(:,ICLC) * AMT )
  END DO
!
  SRC = SRC / SIGN ( MAX ( ABS(OPT), TINY(1.0D0) ), OPT )
  IF ( BFXFLG ) THEN
    IF ( GRAFLG ) THEN
      SRC = SRCBFX ( ITAN, IATM, IDIR, WNOFIN, OPT, SRC, PTH(IPTH)%PSI )
    ELSE
      SRC = SRCBFX ( ITAN, IATM, IDIR, WNOFIN, OPT, SRC )
    END IF
  END IF
!
  TRA = EXP ( -OPT )  
  RAD = RAD * TRA + SRC * ( 1.0D0 - TRA ) 
!
END SUBROUTINE RADLAY
END MODULE RADLAY_SUB
