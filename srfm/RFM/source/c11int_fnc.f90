MODULE C11INT_FNC
CONTAINS
CHARACTER(11) PURE FUNCTION C11INT ( INUM ) 
!
! VERSION
!   25JAN24 Checked.
!   01MAY17 Original. Checked.
!
! DESCRIPTION
!   Write integer as left-adjusted string
!   General purpose module
!   C*11 accommodates any 32bit (I4) integer.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: INUM ! Number to be written
!
! LOCAL VARIABLES
    CHARACTER(11) :: C11 ! Internal version of eventual output 
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  WRITE ( C11, '(I11)' ) INUM
  C11INT = ADJUSTL ( C11 ) 
!
END FUNCTION C11INT
END MODULE C11INT_FNC

