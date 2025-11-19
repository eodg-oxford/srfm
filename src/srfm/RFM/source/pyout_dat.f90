MODULE PYOUT_DAT
!
! VERSION
!   18FEB25 AK Added to support in-memory capture of RFM outputs.
!
! DESCRIPTION
!   Data containers for capturing RFM outputs when running via the Python
!   interface.  These structures retain the atmospheric profile and spectral
!   optical-depth information so that Python callers can access the results
!   without relying on intermediate files.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: PROFILE_RESULT
    LOGICAL :: HAS_DATA = .FALSE.
    INTEGER(I4) :: N_LEVELS = 0
    REAL(R8), ALLOCATABLE :: ALTITUDE(:)
    REAL(R8), ALLOCATABLE :: PRESSURE(:)
    REAL(R8), ALLOCATABLE :: TEMPERATURE(:)
  END TYPE PROFILE_RESULT
!
  TYPE :: SPECTRAL_RESULT
    LOGICAL :: HAS_DATA = .FALSE.
    INTEGER(I4) :: N_LEVELS = 0
    INTEGER(I4) :: N_POINTS = 0
    REAL(R8), ALLOCATABLE :: ALTITUDE(:)
    REAL(R8), ALLOCATABLE :: PRESSURE(:)
    REAL(R8), ALLOCATABLE :: TEMPERATURE(:)
    REAL(R8), ALLOCATABLE :: WNO(:)
    REAL(R8), ALLOCATABLE :: CUMULATIVE(:,:) ! (N_LEVELS, N_POINTS)
  END TYPE SPECTRAL_RESULT
!
  LOGICAL :: CAPTURE_ENABLED = .FALSE.
  TYPE(PROFILE_RESULT) :: PROFILE
  TYPE(SPECTRAL_RESULT), TARGET, ALLOCATABLE :: SPECTRA(:)
!
END MODULE PYOUT_DAT
