MODULE SPCILS_SUB
CONTAINS
SUBROUTINE SPCILS ( ISPC ) 
!
! VERSION
!   07MAR22 AD Bug#31 Average cooling rate spectra as well
!   30NOV18 AD Bug#12 Set IRRFUL = .FALSE. after convolution. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Convolve spectra with ILS
!   Called by RFMSPC for each spectral range if ILS or AVG flags enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FULCOM_DAT ! Full grid data
    USE SPCCOM_DAT ! Spectral range data
    USE FLGCOM_DAT, ONLY: AVGFLG, COOFLG ! spectrally average, cooling rate
    USE TANCOM_DAT, ONLY: MTAN ! No. tangent paths
!
! SUBROUTINES
    USE ILSGRD_SUB ! Calculate ILS fn for irregular grid
    USE ILSINT_SUB ! Interpolate ILS function
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: ISPC ! Index of current spectral range
!
! LOCAL VARIABLES
    LOGICAL     :: IRRGRD ! T=using irregular spectral grid, F=regular grid
    INTEGER(I4) :: I1,I2  ! Fine grid extent of ILS function
    INTEGER(I4) :: IFIN   ! Counter for fine grid points
    INTEGER(I4) :: IFNOUT ! Location of o/p grid points in fine grid
    INTEGER(I4) :: IOUT   ! Counter for output grid points
    INTEGER(I4) :: IP1ILS ! Lower bound of ILS function rel. to output point
    INTEGER(I4) :: IP2ILS ! Upper bound of ILS function rel. to output point
    INTEGER(I4) :: ITAN   ! Counter for output tangent paths
    INTEGER(I4) :: MFIN   ! No. fine grid points per output grid point
    INTEGER(I4) :: NOUT   ! No. output grid points
    REAL(R8), ALLOCATABLE :: FNCILS(:)   ! Interpolated ILS function
    REAL(R8), ALLOCATABLE :: FNCIRR(:)   ! ILS function for irregular grid
    REAL(R8), ALLOCATABLE :: COOILS(:,:) ! Convolved cooling rates
    REAL(R8), ALLOCATABLE :: RADILS(:,:) ! Convolved radiances
    REAL(R8), ALLOCATABLE :: TRAILS(:,:) ! Convolved transmittances
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  MFIN = NINT ( SPC(ISPC)%WNR / WNRFUL ) 
!
  IF ( AVGFLG ) THEN
    IP1ILS = -MFIN
    IP2ILS = MFIN
    ALLOCATE ( FNCILS(IP1ILS:IP2ILS) )
    DO IFIN = -MFIN, MFIN
      FNCILS(IFIN) = DBLE ( MFIN - ABS ( IFIN ) ) 
    END DO
    FNCILS = FNCILS / SUM ( FNCILS )
  ELSE
    CALL ILSINT ( SPC(ISPC)%ILS, WNRFUL, IP1ILS, IP2ILS, FNCILS ) 
  END IF
!
  NOUT = SPC(ISPC)%NPT
  ALLOCATE ( RADILS(NOUT,MTAN), TRAILS(NOUT,MTAN) ) 
  IF ( COOFLG ) ALLOCATE ( COOILS(NOUT,MTAN) )
!
! Allocate WGTGRD for maximum size assuming every fine grid pt within ILS is used
  IRRGRD = SPC(ISPC)%IGD .NE. 0 
  IF ( IRRGRD ) ALLOCATE ( FNCIRR(IP2ILS-IP1ILS+1) )
!
  IFNOUT = 1 + NINT ( ( SPC(ISPC)%WNL - SPC(ISPC)%WXL ) / WNRFUL )
  I1 = 1
  DO IOUT = 1, NOUT
    IF ( IRRGRD ) THEN   ! RADFUL is regular fine grid
      CALL ILSGRD ( IFNOUT, IP1ILS, IP2ILS, FNCILS, I1, I2, FNCIRR ) 
      DO ITAN = 1, MTAN 
        RADILS(IOUT,ITAN) = DOT_PRODUCT ( FNCIRR, RADFUL(I1:I2,ITAN) )
        TRAILS(IOUT,ITAN) = DOT_PRODUCT ( FNCIRR, TRAFUL(I1:I2,ITAN) )
        IF ( COOFLG ) &
          COOILS(IOUT,ITAN) = DOT_PRODUCT ( FNCIRR, COOFUL(I1:I2,ITAN) )
      END DO
    ELSE                      ! RADFUL is irregular grid
      I1 = IFNOUT + IP1ILS
      I2 = IFNOUT + IP2ILS
      DO ITAN = 1, MTAN 
        RADILS(IOUT,ITAN) = DOT_PRODUCT ( FNCILS, RADFUL(I1:I2,ITAN) )
        TRAILS(IOUT,ITAN) = DOT_PRODUCT ( FNCILS, TRAFUL(I1:I2,ITAN) )
        IF ( COOFLG ) & 
          COOILS(IOUT,ITAN) = DOT_PRODUCT ( FNCILS, COOFUL(I1:I2,ITAN) )
      END DO
    END IF
    IFNOUT = IFNOUT + MFIN
  END DO
!
  CALL MOVE_ALLOC ( RADILS, RADFUL )
  CALL MOVE_ALLOC ( TRAILS, TRAFUL )
  IF ( COOFLG ) CALL MOVE_ALLOC ( COOILS, COOFUL ) 
  NFUL = NOUT
  IRRFUL = .FALSE.  ! Full grid is always regular after any convolution
!
  DEALLOCATE ( WNOFUL )
  ALLOCATE ( WNOFUL(NFUL) )
  DO IOUT = 1, NOUT
    WNOFUL(IOUT) = SPC(ISPC)%WNL + ( IOUT - 1 ) * SPC(ISPC)%WNR
  END DO
!
END SUBROUTINE SPCILS
END MODULE SPCILS_SUB
