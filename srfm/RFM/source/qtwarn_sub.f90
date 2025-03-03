MODULE QTWARN_SUB
CONTAINS
SUBROUTINE QTWARN ( IDXMOL, IDXISO, TEM, TEM1, TEM2 )
!
! VERSION
!   02SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Warning messages from QTFCT
!   Called by QTFCT
!   Issues warning messages if TIPS data missing or out of temperature range.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
    USE NAMMOL_FNC ! Return molecule name + i[#iso]
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
  SAVE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IDXMOL ! HITRAN gas ID
    INTEGER(I4), INTENT(IN) :: IDXISO ! HITRAN iso. ID, or -1 = no TIPS for mol.
    REAL(R4),    INTENT(IN) :: TEM    ! Path temperature [K] (if IDXMOL = -1 )
    REAL(R4), OPTIONAL, INTENT(IN) :: TEM1 ! Min Temp [K] in TIPS tabulation
    REAL(R4), OPTIONAL, INTENT(IN) :: TEM2 ! Max Temp [K] in TIPS tabulation
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: MAXWRN = 10 ! Max no of different warning messages
!
! LOCAL VARIABLES
    LOGICAL       :: WARNHI = .FALSE. ! T=T>Tmax warning already issued
    LOGICAL       :: WARNLO = .FALSE. ! T=T<Tmin warning already issued
    INTEGER(I4)   :: IDMWRN(MAXWRN)   ! List of IDXMOL already sent as warnings
    INTEGER(I4)   :: ISOWRN(MAXWRN)   ! List of ISO already issued as warnings
    INTEGER(I4)   :: NWRN = 0         ! No. of diff. warning messages so far
    CHARACTER(80) :: LOGMSG           ! Text for warning message
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! First check for warning messages due to low or high temperature limits
  IF ( PRESENT ( TEM1 ) ) THEN
    IF ( .NOT. WARNLO ) THEN
      LOGMSG = 'W-QTFCT: Extrapolating TIPS below Tmin, T=' // C9REAL(TEM) 
      WARNLO = .TRUE.
      CALL WRTLOG ( LOGMSG )
    END IF
    RETURN
  ELSE IF ( PRESENT ( TEM2 ) ) THEN
    IF ( .NOT. WARNHI ) THEN
      LOGMSG = 'W-QTFCT: Extrapolating TIPS above Tmax, T=' // C9REAL(TEM)
      WARNHI = .TRUE.
      CALL WRTLOG ( LOGMSG ) 
    END IF
    RETURN
  END IF
!
! For other warnings, check if maximum number of warnings already issues
  IF ( NWRN .EQ. MAXWRN ) RETURN
!
! Check if warning has been issued for this particular molec,iso combination
  IF ( ANY ( IDMWRN .EQ. IDXMOL .AND. ISOWRN .EQ. IDXISO ) ) RETURN
!
! Construct warning message
  IF ( IDXISO .EQ. -1 ) THEN         
    LOGMSG = 'W-QTFCT: No TIPS data for ' // NAMMOL(IDXMOL) // &
             ' - assuming T^(3/2) dependence'
  ELSE 
    LOGMSG = 'W-QTFCT: No TIPS data for ' // NAMMOL(IDXMOL,IDXISO) // &
             ' - assuming ~ Isotope#1'
  END IF
  CALL WRTLOG ( LOGMSG )
!
! Save details to avoid repetition
  NWRN = NWRN + 1
  IDMWRN(NWRN) = IDXMOL
  ISOWRN(NWRN) = IDXISO
  IF ( NWRN .EQ. MAXWRN ) &
    CALL WRTLOG ( 'W-QTWARN: Further such warnings suppressed' )
!
END SUBROUTINE QTWARN
END MODULE QTWARN_SUB
