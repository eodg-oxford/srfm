MODULE SPCLUT_SUB
CONTAINS
SUBROUTINE SPCLUT ( FAIL, ERRMSG )
!
! VERSION
!   05JUN24 AD Checked.
!   01MAY17 AD F90 conversion of rfmlut.for. Checked.
! 
! DESCRIPTION    
!   Calculate the absorption using external Abs.Coeff. Tables.
!   Called by RFMSPC for each wide mesh interval if LUT flag enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE FINCOM_DAT ! Finemesh data
    USE LUTCOM_DAT ! TAB LUT data
!
! SUBROUTINES
    USE REALUT_SUB ! Interpolate data from LUT file for given wavenumber
!                  
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL     :: ZERABS ! T=zero absorption for all p,T,q points
    INTEGER(I4) :: ICLC   ! Index in CLCCOM of calculated path
    INTEGER(I4) :: IFIN   ! Fine mesh grid point counter 
    INTEGER(I4) :: ILP    ! Counter for paths using LUT
    INTEGER(I4) :: ILUT   ! Counter for LUT files
    REAL(R4)    :: K      ! Interp.Abs.Coeffs at high, low wavenumbers
    REAL(R8)    :: WNO    ! Waveno of current Fine Mesh grid point
    INTEGER(I4), ALLOCATABLE :: IDX(:) ! Indices of p,T,q pts for interp.
    REAL(R4),    ALLOCATABLE :: TAB(:) ! Spectrally interp. record from TAB file
    REAL(R4),    ALLOCATABLE :: WGT(:) ! Interpolation weights for p,T,q points
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Loop over all LUT for this spectral range
  DO ILUT = 1, NLUT
    ALLOCATE ( TAB(LUT(ILUT)%NX) )
    ALLOCATE ( IDX(LUT(ILUT)%NLK), WGT(LUT(ILUT)%NLK) )
!
! Loop over fine mesh grid points
    DO IFIN = 1, NFIN
      WNO = WNOFIN(IFIN)
!
      CALL REALUT ( ILUT, WNOFIN(IFIN), ZERABS, TAB, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      IF ( ZERABS ) CYCLE
!
      DO ILP = 1, LUT(ILUT)%NLP
        ICLC = LUT(ILUT)%ICL(ILP)
        IDX = LUT(ILUT)%IDX(:,ILP)
        WGT = LUT(ILUT)%WGT(:,ILP) 
        K = DOT_PRODUCT ( TAB(IDX), WGT ) 
! Factor 1.0e4 to convert AMT [kmoles/cm^2]) to [kmoles/m^2] for compatibility
! with pretabulated K coefficients.
        ABSFIN(IFIN,ICLC) = ABSFIN(IFIN,ICLC) + &
                            EXP(K) * CLC(ICLC)%AMT * 1.0E4
      END DO
    END DO
    DEALLOCATE ( TAB, IDX, WGT ) 
  END DO
!
  FAIL = .FALSE.
!
END SUBROUTINE SPCLUT
END MODULE SPCLUT_SUB
