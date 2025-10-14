MODULE QADCOM_DAT
!
! VERSION
!   11DEC23 AD Checked.
!   01JUN17 AD F90 original. Checked.
!
! DESCRIPTION
!   Gaussian quadrature data
!   Loaded in INIQAD, GAUQAD.
!   If VRTFLG is set, ie no hemispherical integration, then default values are 
!   altered.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL VARIABLES
    INTEGER(I4) :: NQAD = 4   ! Default No.pts for Gaussian quadrature
    REAL(R8)    :: RPIQAD     ! Reciprocal pi factor
    REAL(R8), ALLOCATABLE :: XQAD(:) ! cosine(angle) for quadrature
    REAL(R8), ALLOCATABLE :: WQAD(:) ! quadrature weights
!
END MODULE QADCOM_DAT
