MODULE ATMPAR_SUB
CONTAINS
SUBROUTINE ATMPAR ( PARAM, VALUE, FAIL, ERRMSG )
!
! VERSION
!   01NOV19 AD Checked.
!   28JUN18 AD Bug#9 - ensure ATM arrays initialised with HOM flag
!   21JUN17 AD Original.
!
! DESCRIPTION
!   Assign atmospheric profile to single value
!   Called by DRVATM if PARAM=VALUE listed in *ATM section of driver table.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT, ONLY: NATM, LENATM ! Max length of .atm profile label
    USE FLGCOM_DAT, ONLY: HOMFLG       ! T = use homogeneous path
!
! SUBROUTINES
    USE ATMINI_SUB ! Initialise atmospheric profile data in ATMCOM
    USE ATMPRF_SUB ! Load profile into ATMCOM
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: PARAM  ! Profile type
    CHARACTER(*),  INTENT(IN)  :: VALUE  ! Profile value
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL           :: USEPAR ! T=parameter is required, F=not required
    INTEGER(I4)       :: IOS    ! Saved value of IOSTAT for error message
    REAL(R4)          :: RVAL   ! VALUE converted to real number
    CHARACTER(LENATM) :: LABEL  ! Label identifying profile
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  LABEL = PARAM
  IF ( LABEL(1:1) .NE. '*') LABEL = '*'//LABEL(1:LENATM-1)
!
  READ ( VALUE, *, IOSTAT=IOS ) RVAL
  IF ( IOS .NE. 0 ) THEN
    WRITE ( ERRMSG, * ) 'F-ATMPAR: Error reading value from ' // PARAM // &
      '=' // VALUE // '. IOS=', IOS
    FAIL = .TRUE.
    RETURN
  END IF
!
! If homogeneous path calculation it is possible that this is the first entry
! in the *ATM section so ensure ATM arrays are defined with HGT=0.0
  IF ( HOMFLG .AND. NATM .EQ. 0 ) CALL ATMINI ( (/ 0.0 /), .FALSE. ) 
!
  CALL ATMPRF ( LABEL, (/ 1.0 /), (/ RVAL /), USEPAR, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN
  IF ( .NOT. USEPAR ) CALL WRTLOG ( 'W-ATMPAR: Parameter ' // TRIM(LABEL) // &
    ' ignored - not required' ) 
!
END SUBROUTINE ATMPAR
END MODULE ATMPAR_SUB

