MODULE HDRCOM_DAT
!
! VERSION
!   17DEC23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Output header data
!   VIDHDR set in RFM 
!   TXTHDR set in DRVHDR 
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL VARIABLES
    CHARACTER(80) :: TXTHDR ! Text header from driver table
    CHARACTER(11) :: VIDHDR ! RFM Version identifier 
!
END MODULE HDRCOM_DAT
