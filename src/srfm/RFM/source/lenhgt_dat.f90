MODULE LENHGT_DAT
!
! VERSION
!   18JUL26 AD Checked.
!   01AUG25 AD Original. Replaces LENTAN in TANCOM.
!
! DESCRIPTION
!   Max length of height component of RFM output filenames
!   6 allows for sign + 5 digits
! 
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
!
! GLOBAL CONSTANTS
    INTEGER(I4), PARAMETER :: LENHGT = 6
!
END MODULE LENHGT_DAT

