MODULE FULCOM_DAT
!
! VERSION
!   04FEB24 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Full grid data
!   This is initially used for fine grid calculation over the full spectral
!   range but, with spectral convolution, is then reset to the output grid.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL VARIABLES
    LOGICAL     :: IRRFUL ! T = irregular spectral grid, F = regular
    INTEGER(I4) :: IFUL1  ! Index of lower limit of current widemesh interval
    INTEGER(I4) :: IFUL2  ! Index of upper limit of current widemesh interval
    INTEGER(I4) :: NFUL   ! Tot. No points in full output spectra
    REAL(R8)    :: WNLFUL ! Lower Wno [cm-1] in WNOFUL
    REAL(R8)    :: WNRFUL ! Spectral resolution [cm-1] of WNOFUL
    REAL(R8)    :: WNUFUL ! Upper Wno [cm-1] in WNOFUL
    REAL(R8), TARGET, ALLOCATABLE :: COOFUL(:,:) ! [NFUL,MTAN] Cooling rates
    REAL(R8), TARGET, ALLOCATABLE :: OPTFUL(:,:) ! [NFUL,MTAN] Optical depth
    REAL(R8), TARGET, ALLOCATABLE :: RADFUL(:,:) ! [NFUL,MTAN] Radiance
    REAL(R8), TARGET, ALLOCATABLE :: TRAFUL(:,:) ! [NFUL,MTAN] Transmittance
    REAL(R8), TARGET, ALLOCATABLE :: WNOFUL(:)   ! [NFUL] Wavenumber axis [cm-1]
!
END MODULE FULCOM_DAT
