MODULE RUN_ID_DAT
!
! VERSION
!   01DEC23 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   RFM Run ID string
!   Optional string added by user appended to output filenames to distinguish 
!   separate runs. Set in main program. Used in MAKNAM_SUB
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL VARIABLES
    CHARACTER(90) :: RUN_ID = '' ! Optional ID appended to filenames
!
END MODULE RUN_ID_DAT
