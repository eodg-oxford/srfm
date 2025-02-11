MODULE ILSINT_SUB
CONTAINS
SUBROUTINE ILSINT ( IILS, WNRFIN, IP1FIN, IP2FIN, ILSFIN ) 
!
! VERSION
!   19MAY24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Interpolate ILS function
!   Called by SPCILS for each spectral range
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ILSCOM_DAT ! Instrument Lineshape data
!
! SUBROUTINES
    USE INTERP_GEN ! Interpolate array
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN)  :: IILS   ! Index of reqd tabulated ILS function
    REAL(R8),    INTENT(IN)  :: WNRFIN ! Reqd.resln [cm-1] for interp. ILS
    INTEGER(I4), INTENT(OUT) :: IP1FIN ! Lower index of interp.ILS on fine grid
    INTEGER(I4), INTENT(OUT) :: IP2FIN ! Upper index of interp.ILS on fine grid
    REAL(R8), ALLOCATABLE, INTENT(OUT) :: ILSFIN(:) ! Interpolated ILS function
!
! LOCAL VARIABLES
    INTEGER(I4) :: IFIN    ! Counter for fine grid points
    INTEGER(I4) :: IPT     ! Counter for tabulated points
    INTEGER(I4) :: NPTFIN  ! No. points in interpolated ILS
    INTEGER(I4) :: NPTILS  ! No. points in Tabulated ILS
    REAL(R8), ALLOCATABLE :: WNOFIN(:) ! Waveno axis [cm-1] for interpolation
    REAL(R8), ALLOCATABLE :: WNOILS(:) ! Waveno axis [cm-1] of tabulated fn.
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NPTILS = ILS(IILS)%NPT
!
  IP1FIN = NINT ( ILS(IILS)%PT1 / WNRFIN )          ! PT1 is usually negative
  IP2FIN = NINT ( ILS(IILS)%PT2 / WNRFIN )
  NPTFIN = IP2FIN - IP1FIN + 1
  ALLOCATE ( ILSFIN(NPTFIN) )
!
! Assume interpolation only required if no.points is different.
  IF ( NPTFIN .NE. NPTILS ) THEN
    ALLOCATE ( WNOFIN(NPTFIN) ) 
    DO IFIN = 1, NPTFIN
      WNOFIN(IFIN) = ( IP1FIN + IFIN - 1 ) * WNRFIN
    END DO
    ALLOCATE ( WNOILS(NPTILS) ) 
    DO IPT = 1, NPTILS
      WNOILS(IPT) = ILS(IILS)%PT1 + (IPT-1) * ILS(IILS)%PTD
    END DO
    ILSFIN = INTERP ( WNOILS, WNOFIN, ILS(IILS)%FNC )
  ELSE
    ILSFIN = ILS(IILS)%FNC
  END IF
!
  ILSFIN = ILSFIN / SUM ( ILSFIN )
!
END SUBROUTINE ILSINT
END MODULE ILSINT_SUB

