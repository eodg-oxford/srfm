MODULE UPCASE_FNC
CONTAINS
PURE FUNCTION UPCASE ( INPSTR )
!
! VERSION
!   15DEC23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Convert text string to upper case
!   General purpose module.
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
    CHARACTER(LEN(INPSTR))   :: UPCASE ! Variable length output string
!
! LOCAL CONSTANTS
    CHARACTER(26), PARAMETER :: UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
    CHARACTER(26), PARAMETER :: LOWER = 'abcdefghijklmnopqrstuvwxyz'
!
! LOCAL VARIABLES
    INTEGER(I4) :: I, J
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  DO I = 1, LEN(INPSTR)
    J = INDEX ( LOWER, INPSTR(I:I) )
    IF ( J .NE. 0 ) THEN
      UPCASE(I:I) = UPPER(J:J)
    ELSE
      UPCASE(I:I) = INPSTR(I:I)
    END IF
  END DO
!
END FUNCTION UPCASE
END MODULE UPCASE_FNC


