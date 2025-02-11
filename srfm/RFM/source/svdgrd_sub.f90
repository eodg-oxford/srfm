MODULE SVDGRD_SUB
CONTAINS 
SUBROUTINE SVDGRD ( UT, NV, V1, DV, IGO, NG ) 
!
! VERSION
!   03JUL24 AD Checked.
!   01MAY17 AD F90 original. Checked.
! 
! DESCRIPTION    
!   Thin out SVD UT matrix to match (irregular) fine grid
!   Called by REASVD
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FULCOM_DAT ! Full grid data
!                  
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4),    INTENT(INOUT) :: UT(:,:) ! [NL,NV] Transpose of U matrix
    INTEGER(I4), INTENT(IN)    :: NV      ! Original no wno pts in UT
    REAL(R8),    INTENT(IN)    :: V1      ! Wno [cm-1] of first pt in UT
    REAL(R8),    INTENT(IN)    :: DV      ! Wno spacing [cm-1] of UT
    INTEGER(I4), INTENT(OUT)   :: IGO     ! Full grid index offset
    INTEGER(I4), INTENT(OUT)   :: NG      ! No. spectral pts in new UT
!
! LOCAL VARIABLES
    INTEGER(I4) :: IFUL ! Counter for full spectral grid points
    INTEGER(I4) :: IV   ! Counter for original spectral points
    REAL(R8)    :: WNO  ! Waveno [cm-1] of full grid point
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  NG = 0
  IGO = 0
  DO IFUL = 1, NFUL
    WNO = WNOFUL(IFUL)
    IV = 1 + NINT ( ( WNO - V1 ) / DV )      
    IF ( IV .LT. 1 ) THEN
      IGO = IFUL
    ELSE IF ( IV .GT. NV ) THEN
      STOP 'F-SVDGRD: Logical error#1' ! assume SVD spans complete spc.range
    ELSE
      NG = NG + 1
! Assume that WNOFUL never starts below V1 and not on finer grid
      IF ( IV .LT. NG ) STOP 'F-SVDGRD: Logical error#2'
      UT(:,NG) = UT(:,IV) 
    END IF
  END DO
!
END SUBROUTINE SVDGRD
END MODULE SVDGRD_SUB
