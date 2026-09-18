MODULE ADDLEV_SUB
  USE KIND_DAT
  IMPLICIT NONE
  PRIVATE

  INTEGER(I4) :: IATPRV_STATE = 0
  REAL(R4)    :: HGTPRV_STATE = -1.0_R4

  PUBLIC :: ADDLEV
  PUBLIC :: ADDLEV_RESET
CONTAINS
SUBROUTINE ADDLEV ( FIELD, HGTBOA, HGTTOA, HGTLEV, FAIL, ERRMSG )
!
! VERSION
!   03JUL26 AD Checked.
!   01AUG25 AD Renamed/rewritten. Was LEVCHK
!   05OCT24 AD Checked.
!   05AUG19 AD Add USEHGT argument for ATMLEV. Check SETHGT. Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Check and add output altitude level
!   Called by DRVLEV for each field in *LEV section.
!
! GLOBAL DATA
    USE LENHGT_DAT ! Max length of height component of RFM output filenames
    USE LEVCOM_DAT ! Intermediate output levels 
!
! SUBROUTINES
    USE ATMLEV_SUB ! Find/insert atmospheric level for given altitude/pressure
    USE HGTSTR_FNC ! Convert altitude to C*5 string
    USE UPCASE_FNC ! Convert text string to upper case
!
    IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: FIELD  ! Altitude string to be tested
    REAL(R4),      INTENT(IN)  :: HGTBOA ! Altitude [km] of base of atmosphere
    REAL(R4),      INTENT(IN)  :: HGTTOA ! Altitude [km] of top of atmosphere
    REAL(R4),      INTENT(OUT) :: HGTLEV ! Altitude [km] read from FIELD
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)       :: IATM        ! Index of atm level corresp to HGTLEV
    INTEGER(I4)       :: IOS         ! Saved value of IOSTAT for error message
    CHARACTER(LENHGT) :: LEVSTR      ! Altitude converted to string
    TYPE(LEVTYP), ALLOCATABLE :: LEVSAV(:) ! Saved version of LEV during realloc
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .TRUE.
  IF ( NLEV .EQ. 0 ) THEN
    IATPRV_STATE = 0
    HGTPRV_STATE = HGTBOA - 1.0_R4
  END IF
!
  SELECT CASE ( UPCASE ( FIELD )  )
  CASE ( 'TOA' ) 
    HGTLEV = HGTTOA
    TOALEV = .TRUE.
  CASE ( 'BOA' )
    HGTLEV = HGTBOA
    BOALEV = .TRUE.
  CASE DEFAULT
    READ ( FIELD, *, IOSTAT=IOS ) HGTLEV 
    IF ( IOS .NE. 0 ) THEN
      ERRMSG = 'F-ADDLEV: Unreadable value in *LEV section:' // LEVSTR
      RETURN
    ELSE IF ( HGTLEV .GT. HGTTOA ) THEN
      ERRMSG = 'F-ADDLEV: Output level is above top of atmosphere'
      RETURN
    ELSE IF ( HGTLEV .LT. HGTBOA ) THEN
      ERRMSG = 'F-ADDLEV: Output level is below base of atmosphere'
      RETURN   
    END IF
  END SELECT

  IF ( HGTLEV .LT. HGTPRV_STATE ) THEN
    ERRMSG = 'F-ADDLEV: List of output levels not increasing monotonically'
    RETURN
  END IF

  IF ( NLEV .EQ. 0 ) THEN
    NLEV = 1
    ALLOCATE ( LEV(NLEV) )
  ELSE
    CALL MOVE_ALLOC ( LEV, LEVSAV ) 
    NLEV = NLEV + 1
    ALLOCATE ( LEV(NLEV) )
    LEV(1:NLEV-1) = LEVSAV
  END IF

  LEV(NLEV)%HGT = HGTLEV
! Convert value to string that will form part of output filenames and check
! that this is unique
  LEVSTR = HGTSTR ( HGTLEV ) 
  IF ( ANY ( LEV(1:NLEV-1)%STR .EQ. LEVSTR ) ) THEN
    ERRMSG = 'F-ADDLEV: Height not unique in part of output filename: ' &
             // LEVSTR
    RETURN
  END IF
  LEV(NLEV)%STR = LEVSTR
!
! Find level in atmosphere, inserting extra level if required
  CALL ATMLEV ( HGTLEV, .TRUE., IATM )    ! T=LEV specified as altitude
  IF ( IATM .LE. IATPRV_STATE ) STOP 'F-ADDLEV: Logical error'
!
  LEV(NLEV)%IAT = IATM
!
  IATPRV_STATE = IATM
  HGTPRV_STATE = HGTLEV
  FAIL = .FALSE.
!
END SUBROUTINE ADDLEV

SUBROUTINE ADDLEV_RESET()
  IATPRV_STATE = 0
  HGTPRV_STATE = -1.0_R4
END SUBROUTINE ADDLEV_RESET

END MODULE ADDLEV_SUB
