MODULE ATMQAL_SUB
CONTAINS
SUBROUTINE ATMQAL ( FILATM, NQAL, QALPRF, FAIL, ERRMSG )
!
! VERSION
!   16AUG19 AD Original. Adapted from MORSE module. Checked.
!
! DESCRIPTION
!   Extract qualifiers from .atm filename
!   Called by ATMFIL for each file in *ATM section
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE UPCASE_FNC ! Convert text string to upper case
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(INOUT) :: FILATM ! .atm filename + qualifiers
    INTEGER(I4),   INTENT(OUT)   :: NQAL   ! No. of qualifiers (-ve=exclude)
    CHARACTER(*), ALLOCATABLE, &
                   INTENT(OUT)   :: QALPRF(:) ! list of profile labels
    LOGICAL,       INTENT(OUT)   :: FAIL   ! Set TRUE if a fatal error occurs
    CHARACTER(80), INTENT(OUT)   :: ERRMSG ! Error message written if FAIL=T
!
! LOCAL VARIABLES
    LOGICAL     :: ANYQAL ! T=qualifier(s) appended to filename
    LOGICAL     :: LEXCLD ! T=list of excluded profiles
    INTEGER(I4) :: I, J   ! Indices of brackets in FILATM
    INTEGER(I4) :: IQAL   ! Counter for profiles in qualifier string
    CHARACTER(LEN(FILATM)) :: QALSTR ! Part of FILATM containing qualifiers
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
  ERRMSG = ''
  ANYQAL = INDEX ( FILATM, '(' ) .GT. 0
  NQAL = 0
  IF ( .NOT. ANYQAL ) RETURN
!
  I = INDEX ( FILATM, '(' ) 
  J = INDEX ( FILATM, ')' ) 
  QALSTR = UPCASE ( ADJUSTL ( FILATM(I+1:J-1) ) )
  FILATM = FILATM(1:I-1)
!
! For list of excluded profiles, first character is '-' 
  LEXCLD = QALSTR(1:1) .EQ. '-' 
  IF ( LEXCLD ) QALSTR = ADJUSTL ( QALSTR(2:) )
!
  NQAL = 1
  DO I = 1, LEN ( QALSTR ) 
    IF ( QALSTR(I:I) .EQ. ';' ) NQAL = NQAL + 1
    IF ( QALSTR(I:I) .EQ. '*' ) QALSTR(I:I) = ' '  ! will be added later
  END DO
!
  ALLOCATE ( QALPRF(NQAL) )
  QALSTR = ADJUSTL ( QALSTR ) 
  IQAL = 1
  DO 
    I = INDEX ( QALSTR, ';' ) 
    IF ( I .EQ. 0 ) THEN    ! Last profile in qualifier list
      QALPRF(IQAL) = '*' // QALSTR 
      EXIT
    END IF
    QALPRF(IQAL) = '*' // QALSTR(1:I-1) 
    QALSTR = ADJUSTL ( QALSTR(I+1:) )
    IQAL = IQAL + 1
  END DO
!
  IF ( IQAL .NE. NQAL ) STOP 'F-ATMQAL: Logical error'
  IF ( LEXCLD ) NQAL = -NQAL
!
END SUBROUTINE ATMQAL
END MODULE ATMQAL_SUB
