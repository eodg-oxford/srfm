MODULE REJCOM_DAT
!
! VERSION
!   11JUL24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Minimum line strength limits
!   Set by DRVREJ if REJFLG enabled.
!   Used in REAHIT when reading in HITRAN lines.
!   Default value is 0.0, meaning all lines of each molecule are included.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL VARIABLES
    LOGICAL :: USEREJ = .FALSE.        ! T=min strength limits set, F=ignore
    REAL(R4), ALLOCATABLE :: STRREJ(:) ! [NGAS] Min.Line strength for inclusion
!
END MODULE REJCOM_DAT

