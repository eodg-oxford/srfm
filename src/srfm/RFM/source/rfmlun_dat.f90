MODULE RFMLUN_DAT
!
! VERSION
!   25AUG24 AD Checked.
!   11AUG23 AD Remove LUNHIT, and renumber LUNTMP from 4 to 3
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Logical Unit Numbers of RFM files.
!   RFMLUN_RESET closes transient units and rewinds LUNNXT so the shared
!   library can be re-entered safely.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
  SAVE
  PUBLIC :: RFMLUN_RESET
!
! GLOBAL CONSTANTS
    INTEGER(I4), PARAMETER :: LUNLOG = 1  ! LUN for log file
    INTEGER(I4), PARAMETER :: LUNDRV = 2  ! LUN for driver file
    INTEGER(I4), PARAMETER :: LUNTMP = 3  ! LUN for temporarily open files
!
! GLOBAL VARIABLES
    INTEGER(I4) :: LUNNXT = 10 ! Next free lun
!
CONTAINS

  SUBROUTINE RFMLUN_RESET()
    LOGICAL :: IS_OPEN

    INQUIRE ( UNIT=LUNTMP, OPENED=IS_OPEN )
    IF ( IS_OPEN ) CLOSE ( UNIT=LUNTMP )
    LUNNXT = 10
  END SUBROUTINE RFMLUN_RESET

END MODULE RFMLUN_DAT
