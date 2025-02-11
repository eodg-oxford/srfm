MODULE SVDPTH_SUB
CONTAINS
SUBROUTINE SVDPTH ( ISVD, NLNP, LNP1, DLNP, NTEM, TEM1, DTEM )
!
! VERSION
!   03JAN24 AD Checked.
!   01JUL17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   Set SVD path interpolation
!   Called by REASVD.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE SVDCOM_DAT ! SVD-compressed LUT data
!
! SUBROUTINES
    USE SVDWGT_SUB ! Calculate SVD-LUT p,T-axis interpolation weights for path
!                  
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: ISVD ! Index of SVD-LUT
    INTEGER(I4), INTENT(IN) :: NLNP ! No. -lnp axis points
    REAL(R4),    INTENT(IN) :: LNP1 ! 1st -lnp axis point
    REAL(R4),    INTENT(IN) :: DLNP ! -lnp axis interval
    INTEGER(I4), INTENT(IN) :: NTEM ! No. Tem axis points
    REAL(R4),    INTENT(IN) :: TEM1 ! 1st Tem axis point
    REAL(R4),    INTENT(IN) :: DTEM ! Tem axis interval
!
! LOCAL VARIABLES
    INTEGER(I4) :: ICLC ! Counter for calculated paths
    INTEGER(I4) :: IGAS ! Index of absorber in GASCOM
    INTEGER(I4) :: ISP  ! Counter for calc paths for this absorber
    INTEGER(I4) :: NLK  ! No. p,T interpolation points
    INTEGER(I4) :: NSP  ! No. calc paths using this absorber
    INTEGER(I4), ALLOCATABLE :: IDXLKP(:) ! LUT (p,T) points for interpolation
    REAL(R4),    ALLOCATABLE :: WGTLKP(:) ! Interpolation weights
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  IGAS = SVD(ISVD)%IGS
  NSP = COUNT ( CLC%IGS .EQ. IGAS ) 
  SVD(ISVD)%NSP = NSP
  ALLOCATE ( SVD(ISVD)%ICL(NSP) ) 
!
  NLK = 1
  IF ( NLNP .GT. 1 ) NLK = NLK * 2
  IF ( NTEM .GT. 1 ) NLK = NLK * 2
  SVD(ISVD)%NLK = NLK
  ALLOCATE ( SVD(ISVD)%WGT(NLK,NSP), SVD(ISVD)%IDX(NLK,NSP) ) 
  ALLOCATE ( IDXLKP(NLK), WGTLKP(NLK) ) 
!
  ISP = 0
  DO ICLC = 1, NCLC
    IF ( CLC(ICLC)%IGS .NE. IGAS ) CYCLE
    ISP = ISP + 1
    CALL SVDWGT ( (ISP .EQ. NSP), IGAS, NLNP, LNP1, DLNP, NTEM, TEM1, DTEM, &
                  CLC(ICLC)%PRE, CLC(ICLC)%TEM, IDXLKP, WGTLKP )   
    SVD(ISVD)%ICL(ISP)   = ICLC
    SVD(ISVD)%IDX(:,ISP) = IDXLKP
    SVD(ISVD)%WGT(:,ISP) = WGTLKP
  END DO
!
END SUBROUTINE SVDPTH
END MODULE SVDPTH_SUB
