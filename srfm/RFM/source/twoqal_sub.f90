MODULE TWOQAL_SUB
CONTAINS
SUBROUTINE TWOQAL ( QALSTR, WILD, ANYQAL, QAL, FAIL, ERRMSG )
!
! VERSION
!   27AUG24 AD Checked.
!   11AUG23 AD Original. 
!
! DESCRIPTION
!   Extract pair of numbers from '( : )' string
!   General purpose module.
!   If WILD is TRUE, then a wildcard character '*' will maintain the input
!   value of QAL(i). 
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*), INTENT(INOUT) :: QALSTR ! .hit filename + qualifiers
    LOGICAL,         INTENT(IN) :: WILD   ! T=accept wild-card '*' as default
    LOGICAL,        INTENT(OUT) :: ANYQAL ! T=identified pair values
    REAL(R8),     INTENT(INOUT) :: QAL(2) ! Qualifier values
    LOGICAL,        INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error occurs
    CHARACTER(80),  INTENT(OUT) :: ERRMSG ! Error message written if FAIL=T
!
! LOCAL VARIABLES
    INTEGER(I4)            :: I      ! Index of '(' in QALSTR
    INTEGER(I4)            :: J      ! Index of ')' in QALSTR
    INTEGER(I4)            :: K      ! Index of ':' in QALSTR
    INTEGER(I4)            :: IOS    ! Value of IOSTAT
    INTEGER(I4)            :: IQAL   ! Counter for qualifiers (1,2)
    CHARACTER(LEN(QALSTR)) :: SUBSTR ! Part of QALSTR containing qualifiers
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
  ANYQAL = .FALSE.
!
! Look for sequence '(' ... ':' ... ')' as indicating pair of qualifiers
  K = INDEX ( QALSTR, ':' ) 
  IF ( K .EQ. 0 ) RETURN
  I = INDEX ( QALSTR(1:K), '(', .TRUE. ) ! backward search
  IF ( I .EQ. 0 ) RETURN
  J = K + INDEX ( QALSTR(K+1:), ')' ) 
  IF ( J .EQ. K ) RETURN
!
  DO IQAL = 1, 2
    IF ( IQAL .EQ. 1 ) THEN  ! Read first qualifier value
      SUBSTR = QALSTR(I+1:K-1)
    ELSE                     ! Read second qualifier value
      SUBSTR = QALSTR(K+1:J-1)
    END IF
! Check for wildcard character,
    IF ( TRIM(SUBSTR) .EQ. '*' ) THEN
      IF ( .NOT. WILD ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-TWOQAL: Wildcard ''*'' not allowed in qualifier string'
        RETURN
      END IF
! Else read numerical value
    ELSE 
      READ ( SUBSTR, *, IOSTAT=IOS ) QAL(IQAL)
      IF ( IOS .NE. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-TWOQAL: Cannot get number from qualifier string: '//SUBSTR
        RETURN
      END IF     
    END IF
  END DO
!
  ANYQAL = .TRUE.
! Remove '( : )' from QALSTR
  QALSTR = QALSTR(1:I-1) // QALSTR(J+1:) 
 
END SUBROUTINE TWOQAL
END MODULE TWOQAL_SUB
