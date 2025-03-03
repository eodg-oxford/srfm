MODULE REAQAL_SUB
CONTAINS
PURE SUBROUTINE REAQAL ( QALSTR, IQAL, FAIL, ERRMSG, TWOQAL, JQAL )
!
! VERSION
!   19MAR23 AD Checked.
!   07NOV17 AD F90 rewritten. Checked.
!
! DESCRIPTION
!   Decode Qualifier strings '(IQAL)' or '(IQAL:JQAL)'
!   General purpose module.
!   Returns TWOQAL=TRUE if two qualifiers detected
!   Also removes '(...)' from QALSTR.
!   Any qualifier containing '*' is returned as IQAL or JQAL=0
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(INOUT) :: QALSTR ! String to be decoded
    INTEGER(I4),   INTENT(OUT)   :: IQAL   ! First/only value
    LOGICAL,       INTENT(OUT)   :: FAIL   ! Set TRUE if fatal error is detected
    CHARACTER(80), INTENT(OUT)   :: ERRMSG ! Error message written if FAIL
    LOGICAL, OPTIONAL, & 
                   INTENT(OUT)   :: TWOQAL ! Set TRUE if single qualifier found
    INTEGER(I4), OPTIONAL, &
                   INTENT(OUT)   :: JQAL   ! Second qualifier value
!
! LOCAL VARIABLES
    INTEGER(I4) :: IOS ! Saved value of IOSTAT
    INTEGER(I4) :: IPT ! Location of '('
    INTEGER(I4) :: JPT ! Location of ')'
    INTEGER(I4) :: KPT ! Location of ':'
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! Default values if wildcard '*' appears in qualifier string
  IQAL = 0
  IF ( PRESENT ( TWOQAL ) ) JQAL = 0
!
  IPT = INDEX ( QALSTR, '(' ) 
  JPT = INDEX ( QALSTR, ')' ) 
  FAIL = .TRUE.
  IF ( IPT .NE. 1 ) THEN
    ERRMSG = 'F-REAQAL: Qualifier string does not start with ''('''
  ELSE IF ( JPT .LT. IPT ) THEN
    ERRMSG = 'F-REAQAL: Qualifier string does not end with '')'''
  ELSE IF ( JPT .EQ. IPT+1 ) THEN
    ERRMSG = 'F-REAQAL: Qualifier string has no contents'
  ELSE
    FAIL = .FALSE.
  END IF
  IF ( FAIL ) RETURN
!
! Established QALSTR of form '(...)' & at least one character between brackets
!
! Check for ':' character indicating two qualifier values
  KPT = INDEX ( QALSTR(1:JPT), ':' ) 
  IF ( PRESENT ( TWOQAL ) ) THEN
    TWOQAL = KPT .GT. 0
  ELSE IF ( KPT .GT. 0 ) THEN 
    FAIL = .TRUE.
    ERRMSG = 'F-REAQAL: Qualifier string contains unexpected '':'''
    RETURN
  END IF
!
  FAIL = .TRUE.
!
! Single number enclosed by brackets
  IF ( KPT .EQ. 0 ) THEN
    IF ( QALSTR(IPT+1:JPT-1) .NE. '*' ) THEN  ! IQAL set = 0 previously
      READ ( QALSTR(IPT+1:JPT-1), *, IOSTAT=IOS ) IQAL
      IF ( IOS .NE. 0 ) THEN  
        WRITE ( ERRMSG, * ) &
          'F-REAQAL: Error reading contents of Qualifier string. IOSTAT=', IOS
        RETURN
      END IF
    END IF
  ELSE
! Pair of numbers enclosed by brackets, separated by colon 
    IF ( KPT .EQ. IPT+1 ) THEN   ! Read first value
      ERRMSG = 'F-REAQAL: Qualifier string 1st field is empty'
      RETURN
    ELSE IF ( QALSTR(IPT+1:KPT-1) .NE. '*' ) THEN   ! IQAL set = 0 previously
      READ ( QALSTR(IPT+1:KPT-1), *, IOSTAT=IOS ) IQAL
      IF ( IOS .NE. 0 ) THEN  
        WRITE ( ERRMSG, * ) &
         'F-REAQAL: Qualifier string error reading 1st field. IOSTAT=', IOS
        RETURN
      END IF
    END IF
    IF ( KPT .EQ. JPT-1 ) THEN     ! Read second value
      ERRMSG = 'F-REAQAL: Qualifier string 2nd field is empty'
      RETURN
    ELSE IF ( QALSTR(KPT+1:JPT-1) .NE. '*' ) THEN  ! JQAL set = 0 previously
      READ ( QALSTR(KPT+1:JPT-1), *, IOSTAT=IOS ) JQAL
      IF ( IOS .NE. 0 ) THEN
        WRITE ( ERRMSG, * ) & 
         'F-REAQAL: Qualifier string error reading 2nd field. IOSTAT=', IOS
        RETURN
      END IF
    END IF
  END IF
!
  QALSTR = QALSTR(JPT+1:)
  FAIL = .FALSE.
!
END SUBROUTINE REAQAL
END MODULE REAQAL_SUB
