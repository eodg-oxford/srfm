MODULE OBSCOM_DAT
!
! VERSION
!   06OCT24 AD Checked.
!   05AUG19 AD Change ALTOBS to HGTOBS. Add PREOBS. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Observer location data
!   HGTOBS set in DRVOBS if OBS flag T.
!   PREOBS set in DRVOBS if OBS flag T.
!   IATOBS set in OBSCHK if OBS flag T.
!   PSIOBS set in DRVOBS if OBS+GRA flags T.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL VARIABLES
    INTEGER(I4) :: IATOBS = 0 ! Atm level# of observer (1:NATM)
    REAL(R4)    :: HGTOBS     ! Observer Altitude [km]
    REAL(R4)    :: PREOBS     ! Observer Pressure [hPa]
    REAL(R4)    :: PSIOBS     ! Observer horizontal angle [deg]
!
END MODULE OBSCOM_DAT
