MODULE LUTPTH_SUB
CONTAINS
SUBROUTINE LUTPTH ( ILUT, IGAS, NLKP )
!
! VERSION
!   19JUN24 AD Checked.
!   04OCT17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   Set LUT path interpolation
!   Called by LUTTAB for each LUT file/spec.rang
!
! VARIABLE KINDS
    USE KIND_DAT
!                  
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE LUTCOM_DAT ! TAB LUT data
    USE ATMCOM_DAT, ONLY: LINVMR ! T= linear interp, F=log interp
!
! SUBROUTINES
    USE LUTWGT_SUB ! Calculate LUT p,T,q-axis interpolation factors for path
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: ILUT ! Index of LUT
    INTEGER(I4), INTENT(IN) :: IGAS ! Index of absorber
    INTEGER(I4), INTENT(IN) :: NLKP ! No. p,T,q interpolation points
!
! LOCAL VARIABLES
    LOGICAL     :: LININT ! T=linear interpolation in VMR for this absorber
    INTEGER(I4) :: ICLC   ! Counter for calculated paths
    INTEGER(I4) :: ILP    ! Counter for calc paths using this absorber
    INTEGER(I4) :: NLP    ! Tot. no. calculated paths for this absorber
    INTEGER(I4), ALLOCATABLE :: IDXLKP(:) ! Interpolation points
    REAL(R4),    ALLOCATABLE :: WGTLKP(:) ! Interpolation weights
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  NLP = COUNT ( CLC%IGS .EQ. IGAS ) 
  LUT(ILUT)%NLP = NLP
  ALLOCATE ( LUT(ILUT)%ICL(NLP) )
!
  LUT(ILUT)%NLK = NLKP
  ALLOCATE ( LUT(ILUT)%WGT(NLKP,NLP), LUT(ILUT)%IDX(NLKP,NLP) )
  ALLOCATE ( WGTLKP(NLKP), IDXLKP(NLKP) ) 
!
  LININT = LINVMR(IGAS) 
!
  ILP = 0
  DO ICLC = 1, NCLC
    IF ( CLC(ICLC)%IGS .NE. IGAS ) CYCLE
    ILP = ILP + 1
    CALL LUTWGT ( (ILP.EQ.NLP), LININT, IGAS, CLC(ICLC)%PRE, CLC(ICLC)%TEM, &
                  CLC(ICLC)%PPA, IDXLKP, WGTLKP ) 
    LUT(ILUT)%ICL(ILP)   = ICLC
    LUT(ILUT)%IDX(:,ILP) = IDXLKP
    LUT(ILUT)%WGT(:,ILP) = WGTLKP
  END DO
!
END SUBROUTINE LUTPTH
END MODULE LUTPTH_SUB
