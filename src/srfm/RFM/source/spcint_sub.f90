 MODULE SPCINT_SUB
CONTAINS
SUBROUTINE SPCINT
!
! VERSION
!   07DEC23 AD Checked.
!   01MAY17 AD F90 conversion of rfmint.for. Checked.
!
! DESCRIPTION
!   Interpolate spectra to regular grid
!   Called by RFMSPC for each spectral range if GRD flag enabled 
!   but without AVG or ILS flags.
!   Assumes regular output fine grid required with GRD flag.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FULCOM_DAT ! Full grid data
    USE FLGCOM_DAT, ONLY:GRDFLG ! Option flags
!
! SUBROUTINES
    USE INTERP_GEN ! Interpolate array
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IFUL   ! Counter for regular grid points
    INTEGER(I4) :: ITAN   ! Counter for output tangent paths
    INTEGER(I4) :: MTAN   ! No. tangent paths in arrays
    INTEGER(I4) :: NIRR   ! No. irreg spectral grid points
    REAL(R8), ALLOCATABLE :: DATIRR(:,:) ! Original irreg gridded data
    REAL(R8), ALLOCATABLE :: WNOIRR(:)   ! Original irregular Wavenumber grid
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( .NOT. GRDFLG ) RETURN
! 
! Move irregular arrays to *IRR
  CALL MOVE_ALLOC ( WNOFUL, WNOIRR ) 
  NIRR = NFUL
!
! Set up regular wavenumber grid
  NFUL = 1 + NINT ( ( WNUFUL - WNLFUL ) / WNRFUL ) 
  ALLOCATE ( WNOFUL(NFUL) )
  WNOFUL = (/ ( WNLFUL + (IFUL-1)*WNRFUL, IFUL = 1, NFUL ) /)
!
  MTAN = SIZE ( OPTFUL, 2 ) 
!
  CALL MOVE_ALLOC ( OPTFUL, DATIRR ) 
  ALLOCATE ( OPTFUL(NFUL,MTAN) ) 
  DO ITAN = 1, MTAN
    OPTFUL(:,ITAN) = INTERP ( WNOIRR, WNOFUL, DATIRR(:,ITAN) )
  END DO
!
  CALL MOVE_ALLOC ( RADFUL, DATIRR ) 
  ALLOCATE ( RADFUL(NFUL,MTAN) ) 
  DO ITAN = 1, MTAN
    RADFUL(:,ITAN) = INTERP ( WNOIRR, WNOFUL, DATIRR(:,ITAN) )
  END DO
!
  CALL MOVE_ALLOC ( TRAFUL, DATIRR ) 
  ALLOCATE ( TRAFUL(NFUL,MTAN) ) 
  DO ITAN = 1, MTAN
    TRAFUL(:,ITAN) = INTERP ( WNOIRR, WNOFUL, DATIRR(:,ITAN) )
  END DO
!
  IF ( ALLOCATED ( COOFUL ) ) THEN  ! interpolate cooling rate spectra if set
    CALL MOVE_ALLOC ( COOFUL, DATIRR ) 
    ALLOCATE ( COOFUL(NFUL,MTAN) ) 
    DO ITAN = 1, MTAN
      COOFUL(:,ITAN) = INTERP ( WNOIRR, WNOFUL, DATIRR(:,ITAN) )
    END DO
  END IF
!
  IRRFUL = .FALSE.
!
END SUBROUTINE SPCINT
END MODULE SPCINT_SUB
