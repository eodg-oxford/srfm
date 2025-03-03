MODULE FLXLAY_SUB
CONTAINS
SUBROUTINE FLXLAY ( IATM, IDIR, SRCPTB, XQAD, RAD, OPTLAY )
!
! VERSION
!   07NOV23 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Radiative flux calculation through a layer
!   Called by FLXATM
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FINCOM_DAT ! Finemesh data
    USE PTHCOM_DAT ! Path segment data
    USE ATMCOM_DAT, ONLY: NATM, NVMR ! No.different atmos. levels, absorbers
    USE FLGCOM_DAT, ONLY: BFXFLG ! T = Source Fn varies with layer optical depth
    USE PTBCON_DAT, ONLY: PTBSRC ! Fraction for Source Fn perturbation
!
! SUBROUTINES
    USE IDXPTH_FNC ! Index in PTHCOM of tan/atm/gas/dir
    USE PLANCK_FNC ! Planck function
    USE SRCBFX_FNC ! Modified source function for BFX flag
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN)    :: IATM      ! Index of atmos.layer
    INTEGER(I4), INTENT(IN)    :: IDIR      ! +1 = upward path, -1 = down path
    LOGICAL,     INTENT(IN)    :: SRCPTB    ! T=Perturb source function
    REAL(R8),    INTENT(IN)    :: XQAD(:)   ! 1/sec(theta) scale factors
    REAL(R8),    INTENT(INOUT) :: RAD(:,:)  ! Radiances along theta paths
    REAL(R8),    INTENT(OUT)   :: OPTLAY(:) ! Layer optical depth
!
! LOCAL VARIABLES
    INTEGER(I4) :: ICLC  ! Index of calculated path 
    INTEGER(I4) :: IPTH  ! Path counter
    INTEGER(I4) :: IQAD  ! Counter over sec theta quadrature angles
    INTEGER(I4) :: IVMR  ! Absorber counter
    INTEGER(I4) :: NQAD  ! No. different zenith angles for quadrature
    REAL(R4)    :: AMT   ! Absorber amount in path [kmol/cm^2]
    REAL(R4)    :: TEM   ! Path CG temperature [K]
    REAL(R8)    :: SRCLAY(SIZE(RAD,1)) ! Absorption-wgtd source fn. for layer
    REAL(R8)    :: TRALAY(SIZE(RAD,1)) ! Transmittance of layer
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NQAD = SIZE ( RAD, 2 )      
!
  SRCLAY = 0.0D0
  OPTLAY = 0.0D0
! 
! For top level of atmosphere, return with RAD unchanged and OPTLAY = 0
  IF ( IATM .EQ. NATM ) RETURN
!
  DO IVMR = 1, NVMR
    IPTH = IDXPTH ( 1, IATM, IVMR )
    ICLC = PTH(IPTH)%ICL
    AMT  = PTH(IPTH)%AMT
    TEM  = PTH(IPTH)%TEM
    SRCLAY = SRCLAY + DBLE ( ABSFIN(:,ICLC) * AMT ) * PLANCK ( TEM, WNOFIN ) 
    OPTLAY = OPTLAY + DBLE ( ABSFIN(:,ICLC) * AMT )
  END DO
!
  SRCLAY = SRCLAY / SIGN ( MAX ( ABS(OPTLAY), TINY(1.0D0) ), OPTLAY )
!
  IF ( SRCPTB ) SRCLAY = SRCLAY * DBLE ( 1.0 + PTBSRC ) 
!
  DO IQAD = 1, SIZE ( RAD, 2 ) 
    TRALAY = EXP ( -OPTLAY / XQAD(IQAD) )  
    IF ( BFXFLG ) &
      SRCLAY = SRCBFX ( 1, IATM, IDIR, WNOFIN, OPTLAY/XQAD(IQAD), SRCLAY )
    RAD(:,IQAD) = RAD(:,IQAD) * TRALAY + SRCLAY * ( 1.0D0 - TRALAY ) 
  END DO
!
END SUBROUTINE FLXLAY
END MODULE FLXLAY_SUB

 
