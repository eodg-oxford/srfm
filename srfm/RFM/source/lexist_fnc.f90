MODULE LEXIST_FNC
CONTAINS
LOGICAL FUNCTION LEXIST ( FILNAM )
!
! VERSION
!   12DEC23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Check if file exists
!   General purpose module.
!   Just a wrapper around the Fortran INQUIRE statement allowing the
!   test to be incorporated into a line of code as a logical variable.
!
  IMPLICIT NONE 
!
! ARGUMENTS
    CHARACTER(*), INTENT(IN) :: FILNAM ! File name to be tested
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  INQUIRE ( FILE=FILNAM, EXIST=LEXIST ) 
!
END FUNCTION LEXIST
END MODULE LEXIST_FNC

