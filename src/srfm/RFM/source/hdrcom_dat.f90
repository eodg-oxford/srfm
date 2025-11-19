MODULE HDRCOM_DAT
!
! VERSION
!   17DEC23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Output header data.
!   VIDHDR is set in RFM and TXTHDR set in DRVHDR.
!   HDRCOM_RESET clears the cached header between repeated runs.
!
  IMPLICIT NONE
  SAVE
  PUBLIC :: HDRCOM_RESET
!
! GLOBAL VARIABLES
    CHARACTER(80) :: TXTHDR ! Text header from driver table
    CHARACTER(11) :: VIDHDR ! RFM Version identifier 
!
CONTAINS

  SUBROUTINE HDRCOM_RESET()
    TXTHDR = ''
  END SUBROUTINE HDRCOM_RESET

END MODULE HDRCOM_DAT
