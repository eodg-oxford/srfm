MODULE LENREC_FNC
CONTAINS
INTEGER(I4) FUNCTION LENREC ( LUN )
!
! VERSION
!   21AUG24 AD Checked.
!   11AUG23 AD Original.
!
! DESCRIPTION
!   Determine length of next record in file
!   General purpose module.
!   After reading, returns pointer to start of record.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: LUN ! Logical Unit Number
!
! LOCAL VARIABLES
    INTEGER(I4)  :: L  ! Local counter for characters in record
    CHARACTER(1) :: C  ! Dummy character for reading
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  L = 0
  DO
    READ ( LUN, '(A)', ADVANCE='NO', EOR=100 ) C
    L = L + 1
  END DO
!
100 CONTINUE
  BACKSPACE ( LUN ) 
  LENREC = L
!
END FUNCTION LENREC
END MODULE LENREC_FNC

