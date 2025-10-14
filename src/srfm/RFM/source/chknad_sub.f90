MODULE CHKNAD_SUB
CONTAINS
SUBROUTINE CHKNAD ( TANTST, USRELE, FAIL, ERRMSG )
!
! VERSION
!   21JAN24 Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Check specification of Nadir/Zenith viewing paths
!   Called by TANCHK for each path with NAD or ZEN flags enabled. 
!   For zenith or nadir viewing TANTST can either be elevation angle (restricted
!   to +/- 90 degrees) or airmass/sec(zen.angle) ), restricted to .GE. 1.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT, ONLY: NADFLG, ZENFLG ! Flags for nadir/zenith viewing
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
!
  IMPLICIT NONE 
!
! ARGUMENTS
    REAL(R4),      INTENT(IN)  :: TANTST ! Value of sec(theta) or elev. to test
    LOGICAL,       INTENT(IN)  :: USRELE ! T=TANTST is elev. angle, F=sec
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    REAL(R4),     PARAMETER :: ELEMIN = 0.1   ! Min elev.ang [deg] 
    CHARACTER(3), PARAMETER :: ELESTR = '0.1' ! 0.1 equivalent to ~600 airmasses
!
! LOCAL VARIABLES
    CHARACTER(9) :: TANSTR ! TANTST written as string for error message.
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  FAIL = .FALSE.
  TANSTR = C9REAL ( TANTST ) 
!
  IF ( USRELE ) THEN            ! elevation angle
    IF ( ZENFLG .AND. ( TANTST .GT. 90.0 .OR. TANTST .LT. ELEMIN ) ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-CHKNAD: Specified Elev. Ang. outside ZEN range +' // &
               ELESTR // ':+90deg, value=' // TANSTR
    ELSE IF ( NADFLG .AND. ( TANTST .LT. -90.0 .OR. TANTST .GT. -ELEMIN ) ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-CHKNAD: Specified Elev. Ang. outside NAD range -' // &
               ELESTR // ':-90deg, value=' // TANSTR
    END IF
!
  ELSE                            ! airmass 
    IF ( TANTST .LT. 1.0 ) THEN  
      FAIL = .TRUE.
      ERRMSG = 'F-CHKNAD: Specified NADir/ZENith view sec(theta) < 1, value =' &
               // TANSTR
    END IF
  END IF
!
END SUBROUTINE CHKNAD
END MODULE CHKNAD_SUB
