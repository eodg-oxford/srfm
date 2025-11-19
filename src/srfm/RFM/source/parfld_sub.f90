MODULE PARFLD_SUB
CONTAINS
PURE SUBROUTINE PARFLD ( RECORD, GOTPAR, PARAM, VALUE )
!
! VERSION
!   22SEP24 AD Checked.
!   02NOV17 AD F90 version. Checked.
!
! DESCRIPTION
!   Extract Parameter=Value string from record
!   General purpose module.
!   'PARAM' is converted to upper case.
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
    CHARACTER(*), INTENT(IN)  :: RECORD ! Field containing param=value
    LOGICAL,      INTENT(OUT) :: GOTPAR ! T=extracted PARAM and VALUE
    CHARACTER(*), INTENT(OUT) :: PARAM  ! PARAM extracted from record
    CHARACTER(*), INTENT(OUT) :: VALUE  ! VALUE extracted from record
!
! LOCAL VARIABLES
    INTEGER(I4)            :: IPT ! Location of '=' sign in RECORD
    INTEGER(I4)            :: LPT ! Loc. of last non-blank character in RECORD
    CHARACTER(LEN(RECORD)) :: LREC ! Left adjusted version of RECORD
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  GOTPAR = .FALSE.
  PARAM = ''
  VALUE = ''
!
  LREC = ADJUSTL ( RECORD ) 
  LPT = LEN_TRIM ( LREC )
  IPT = INDEX ( LREC, '=' )
  IF ( IPT .LE. 1 .OR. IPT .EQ. LPT ) RETURN
!
  GOTPAR = .TRUE.
  PARAM = UPCASE ( LREC(1:IPT-1) )
  VALUE = ADJUSTL ( LREC(IPT+1:LPT) ) 
!
END SUBROUTINE PARFLD
END MODULE PARFLD_SUB
