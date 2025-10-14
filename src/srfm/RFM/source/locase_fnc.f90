MODULE LOCASE_FNC
CONTAINS
PURE FUNCTION LOCASE ( INPSTR )
!
! VERSION
!   09FEB24 AD Checked.
!   01MAY17 AD F90 conversion from F77 subroutine. Checked.
!
! DESCRIPTION
!   Convert text string to lower case
!   General purpose module
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
!
! ARGUMENT
    CHARACTER(*), INTENT(IN) :: INPSTR ! Input string (mixed case)
!
! FUNCTION TYPE
    CHARACTER ( LEN ( INPSTR ) ) :: LOCASE
!
! LOCAL CONSTANTS
    CHARACTER(26), PARAMETER :: LOWER = 'abcdefghijklmnopqrstuvwxyz' 
    CHARACTER(26), PARAMETER :: UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!
! LOCAL VARIABLES
    INTEGER(I4) :: I, J
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  DO I = 1, LEN(INPSTR)
    J = INDEX ( UPPER, INPSTR(I:I) )
    IF ( J .NE. 0 ) THEN
      LOCASE(I:I) = LOWER(J:J)
    ELSE
      LOCASE(I:I) = INPSTR(I:I)
    END IF
  END DO
!
END FUNCTION LOCASE
END MODULE LOCASE_FNC


