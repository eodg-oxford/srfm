MODULE SPCSVD_SUB
CONTAINS 
SUBROUTINE SPCSVD 
!
! VERSION
!   02DEC23 AD Checked.
!   01MAY17 AD F90 conversion of rfmsvd.for. Checked.
! 
! DESCRIPTION    
!   Calculate absorption using SVD-compressed Look-Up Tables
!   Called by RFMSPC for each wide mesh interval if SVD flag enabled
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE FINCOM_DAT ! Finemesh data
    USE SVDCOM_DAT ! SVD-compressed LUT data
!                  
  IMPLICIT NONE
!
! LOCAL VARIABLES
    LOGICAL     :: USELOG ! T=Tabulated function is ln(k)
    LOGICAL     :: USE4RT ! T=Tabulated function is SQRT(SQRT(k))
    INTEGER(I4) :: ICLC   ! Index of calc path in CLCCOM
    INTEGER(I4) :: IGAS   ! Index of absorber
    INTEGER(I4) :: ISP    ! Counter for calc paths for each SVD
    INTEGER(I4) :: ISVD   ! Counter for SVD-LUT data sets
    INTEGER(I4) :: IV1    ! First U matrix spectral point to use
    INTEGER(I4) :: IV2    ! Last U matrix spectral point to use
    INTEGER(I4) :: NL     ! No. singular vectors to use
    INTEGER(I4),  POINTER :: IDX(:)  ! Indices of interp points in p,T axis
    REAL(R4)              :: K(NFIN) ! Absorption coefficient [moles/m2]
    REAL(R4), ALLOCATABLE :: KW(:)   ! Interpolated K for each path
    REAL(R4),     POINTER :: WGT(:)  ! Weights for interpolation points
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  DO ISVD = 1, NSVD
    IGAS = SVD(ISVD)%IGS
    USELOG = SVD(ISVD)%TAB .EQ. 'LOG'
    USE4RT = SVD(ISVD)%TAB .EQ. '4RT'
    IV1 = IF1FIN - SVD(ISVD)%IGO 
    IV2 = IF2FIN - SVD(ISVD)%IGO
    NL = SVD(ISVD)%NL
    ALLOCATE ( KW(NL) )
!
    DO ISP = 1, SVD(ISVD)%NSP
      ICLC = SVD(ISVD)%ICL(ISP)
      IDX  => SVD(ISVD)%IDX(:,ISP)
      WGT  => SVD(ISVD)%WGT(:,ISP)
      KW = MATMUL ( SVD(ISVD)%K(1:NL,IDX), WGT ) 
      K  = MATMUL ( SVD(ISVD)%U(IV1:IV2,1:NL), KW ) 
      IF ( USELOG ) K = EXP ( K ) 
      IF ( USE4RT ) K = K**4
! Factor 1.0e7 to convert AMT [kmoles/cm^2]) to [moles/m^2] for compatibility
! with pretabulated K coefficients.
      ABSFIN(:,ICLC) = ABSFIN(:,ICLC) + K * CLC(ICLC)%AMT * 1.0E7
    END DO
    DEALLOCATE ( KW ) 
  END DO
!
END SUBROUTINE SPCSVD
END MODULE SPCSVD_SUB

