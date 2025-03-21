MODULE RFMLUN_DAT
!
! VERSION
!   25AUG24 AD Checked.
!   11AUG23 AD Remove LUNHIT, and renumber LUNTMP from 4 to 3
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Logical Unit Numbers of RFM files
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL CONSTANTS
    INTEGER(I4), PARAMETER :: LUNLOG = 1  ! LUN for log file
    INTEGER(I4), PARAMETER :: LUNDRV = 2  ! LUN for driver file
    INTEGER(I4), PARAMETER :: LUNTMP = 3  ! LUN for temporarily open files
!
! GLOBAL VARIABLES
    INTEGER(I4) :: LUNNXT = 10 ! Next free lun
!
END MODULE RFMLUN_DAT
