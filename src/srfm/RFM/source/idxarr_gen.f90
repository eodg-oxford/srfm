MODULE IDXARR_GEN
!
! VERSION
!   16JUL26 AD Checked.
!   01AUG25 AD Add LOGICAL
!   09JUN24 AD Checked.
!   10AUG21 AD Original. Checked.
!
! DESCRIPTION
!   Index of value within an array, else 0
!   General purpose module.
!   Basically a simplified version of the intrinsic function FINDLOC
!   introduced in Fortran 2008.
!
! VARIABLE KINDS
    USE KIND_DAT
!
INTERFACE IDXARR
  MODULE PROCEDURE IDXARR_C, IDXARR_I, IDXARR_L
END INTERFACE

CONTAINS

PURE FUNCTION IDXARR_C ( XARR, X )
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*), INTENT(IN) :: XARR(:) ! Array to be searched
    CHARACTER(*), INTENT(IN) :: X       ! Value to be located
!
! FUNCTION TYPE
    INTEGER(I4) :: IDXARR_C  ! Function returns integer
!
! LOCAL VARIABLES
    INTEGER(I4) :: I    ! array index
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( ANY ( XARR .EQ. X ) ) THEN
    DO I = 1, SIZE(XARR)
      IF ( XARR(I) .EQ. X ) THEN
        IDXARR_C = I
        EXIT
      END IF
    END DO
!    STOP 'F-IDXARR_C: Logical error'
  ELSE
    IDXARR_C = 0
  END IF
!
END FUNCTION IDXARR_C

PURE FUNCTION IDXARR_I ( XARR, X )
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: XARR(:) ! Array to be searched
    INTEGER(I4), INTENT(IN) :: X       ! Value to be located
!
! FUNCTION TYPE
    INTEGER(I4) :: IDXARR_I  ! Function returns integer
!
! LOCAL VARIABLES
    INTEGER(I4) :: I    ! array index
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( ANY ( XARR .EQ. X ) ) THEN
    DO I = 1, SIZE(XARR)
      IF ( XARR(I) .EQ. X ) THEN 
        IDXARR_I = I
        EXIT
      END IF
    END DO
!    STOP 'F-IDXARR_I: Logical error'
  ELSE
    IDXARR_I = 0
  END IF
!
END FUNCTION IDXARR_I

PURE FUNCTION IDXARR_L ( XARR, X )
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL, INTENT(IN) :: XARR(:) ! Array to be searched
    LOGICAL, INTENT(IN) :: X       ! Value to be located
!
! FUNCTION TYPE
    INTEGER(I4) :: IDXARR_L  ! Function returns integer
!
! LOCAL VARIABLES
    INTEGER(I4) :: I    ! array index
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( ANY ( XARR .EQV. X ) ) THEN
    DO I = 1, SIZE(XARR)
      IF ( XARR(I) .EQV. X ) THEN 
        IDXARR_L = I
        EXIT
      END IF
    END DO
!    STOP 'F-IDXARR_L: Logical error'
  ELSE
    IDXARR_L = 0
  END IF
!
END FUNCTION IDXARR_L

END MODULE IDXARR_GEN
