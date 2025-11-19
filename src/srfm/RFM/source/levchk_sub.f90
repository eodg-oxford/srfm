MODULE LEVCHK_SUB
  USE KIND_DAT
  REAL(R4), SAVE :: LEVPRV_STATE = 0.0_R4
CONTAINS
SUBROUTINE LEVCHK ( LEVSTR, FAIL, ERRMSG )
!
! VERSION
!   05OCT24 AD Checked.
!   05AUG19 AD Add USEHGT argument for ATMLEV. Check SETHGT. Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Check if string is legal altitude level
!   Called by DRVLEV for each field in *LEV section.
!   This loads altitude levels into LEVCOM arrays, using NLEV as a counter.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LEVCOM_DAT ! Intermediate output levels
    USE ATMCOM_DAT, ONLY: SETHGT, HGTSFC, HGTTOA ! Profile altitude [km] limits
!
! SUBROUTINES
    USE ATMLEV_SUB ! Find/insert atmospheric level for given altitude
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: LEVSTR ! Altitude string to be tested
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)    :: IATM             ! Index of atm level corresp to LEVTST
    INTEGER(I4)    :: IATPRV = 0       ! Index of previous atm. profile level
    INTEGER(I4)    :: IOS              ! Saved value of IOSTAT for error message
    REAL(R4)       :: LEVTST           ! Altitude value to be tested
    TYPE(LEVTYP), ALLOCATABLE :: LEVSAV(:) ! Saved version of LEV during realloc
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( NLEV .EQ. 0 ) THEN
    LEVPRV_STATE = HGTSFC - 1.0
    IATPRV = 0
  END IF
!
  READ ( LEVSTR, *, IOSTAT=IOS ) LEVTST 
  FAIL = .TRUE.
  IF ( IOS .NE. 0 ) THEN
    ERRMSG = 'F-LEVCHK: Unreadable value in *LEV section:' // LEVSTR
  ELSE IF ( ABS ( LEVTST ) .GE. 99999.5 ) THEN
    ERRMSG = 'F-LEVCHK: Cannot handle values .GE. 99999.5'
! Check that altitude profile has been specified
  ELSE IF ( .NOT. SETHGT ) THEN
    ERRMSG = 'F-LEVCHK: *HGT profile has not been supplied'
! Check that value lies within atmosphere
  ELSE IF ( LEVTST .GT. HGTTOA ) THEN
    ERRMSG = 'F-LEVCHK: Output level is above top of atmosphere'
  ELSE IF ( LEVTST .LT. HGTSFC ) THEN
    ERRMSG = 'F-LEVCHK: Output level is below base of atmosphere'
  ELSE IF ( LEVTST .LT. LEVPRV_STATE ) THEN
    ERRMSG = 'F-LEVCHK: List of output levels not increasing monotonically'
  ELSE
    FAIL = .FALSE.
  END IF
  IF ( FAIL ) RETURN
!
! Find level in atmosphere, inserting extra level if required
  CALL ATMLEV ( LEVTST, .TRUE., IATM )    ! T=LEV specified as altitude
  IF ( IATM .LE. IATPRV ) STOP 'F-LEVCHK: Logical error'
!
  IF ( ALLOCATED ( LEV ) ) CALL MOVE_ALLOC ( LEV, LEVSAV ) 
  NLEV = NLEV + 1
  ALLOCATE ( LEV(NLEV) )
  IF ( ALLOCATED ( LEVSAV ) ) LEV(1:NLEV-1) = LEVSAV
  LEV(NLEV)%IAT = IATM
  LEV(NLEV)%HGT = LEVTST
!
  IATPRV = IATM
  LEVPRV_STATE = LEVTST
!
END SUBROUTINE LEVCHK

SUBROUTINE LEVCHK_RESET()
  LEVPRV_STATE = 0.0_R4
END SUBROUTINE LEVCHK_RESET

END MODULE LEVCHK_SUB
