MODULE HITQAL_SUB
CONTAINS
SUBROUTINE HITQAL ( NAMHIT, IDMQAL, LEXCLD, FAIL, ERRMSG )
!
! VERSION
!   14AUG24 AD Checked.
!   11AUG23 AD Original. 
!
! DESCRIPTION
!   Extract qualifiers from .hit filename
!   Called by HITFIL for each file in *HIT section
!   Set appropriate element of IDMQAL = T if molecule name or molec# found
!   Sets LEXCLD=T if first item in brackets is a minus sign.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT, ONLY: LENGAS ! Maximum length of molecule nam
!
! SUBROUTINES
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE MOLIDX_SUB ! Give molecule name for HITRAN/RFM index, or vice-versa
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*), INTENT(INOUT) :: NAMHIT    ! .hit filename + qualifiers
    LOGICAL,        INTENT(OUT) :: IDMQAL(:) ! T=molec.listed, F=not listed 
    LOGICAL,        INTENT(OUT) :: LEXCLD    ! T=list of exclusions, F=incl.
    LOGICAL,        INTENT(OUT) :: FAIL      ! Set TRUE if a fatal error occurs
    CHARACTER(80),  INTENT(OUT) :: ERRMSG    ! Error message written if FAIL=T
!
! LOCAL VARIABLES
    INTEGER(I4)            :: I, J   ! Indices of brackets in NAMHIT
    INTEGER(I4)            :: IDXMOL ! HITRAN molecule#
    INTEGER(I4)            :: IOS    ! Value of IOSTAT
    CHARACTER(LENGAS)      :: NAMMOL ! Name of molecule
    CHARACTER(LEN(NAMHIT)) :: QALSTR ! Part of NAMHIT containing qualifiers
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
  IDMQAL = .FALSE.
!
  I = INDEX ( NAMHIT, '(' ) 
  IF ( I .EQ. 0 ) RETURN             ! No qualifier string present
!
  QALSTR = ADJUSTL ( NAMHIT(I+1:) )    ! extract '...)', exclude first '('
  NAMHIT = NAMHIT(1:I-1)             ! reduce NAMHIT to just filename part
!
  J = LEN_TRIM ( QALSTR ) 
  IF ( QALSTR(J:J) .NE. ')' ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-HITQAL: ''(...'' string after filename does not end with '')'''
    RETURN
  END IF
  QALSTR(J:J) = ';'  ! replace ')' with ';'
!
! For list of excluded molecules, first character is '-' 
  LEXCLD = ( QALSTR(1:1) .EQ. '-' ) 
  IF ( LEXCLD ) QALSTR = QALSTR(2:)  ! remove '-' and shift left
!
  I = 1
  DO J = 2, LEN_TRIM ( QALSTR )    ! should finish with ';'
    IF ( QALSTR(J:J) .EQ. ';' ) THEN
      NAMMOL = TRIM ( QALSTR(I:J-1) )
      IDXMOL = 0
      CALL MOLIDX ( IDXMOL, NAMMOL ) ! Try reading as formula eg 'h2o'
      IF ( IDXMOL .EQ. 0 ) THEN      ! Try reading as a number
        READ ( QALSTR(I:J-1), *, IOSTAT=IOS ) IDXMOL
        IF ( IOS .NE. 0 ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-HITQAL: Failed to identify a molecule from qualifier: ' &
                   // QALSTR(I:J-1) 
          RETURN
        END IF
      END IF
      IF ( IDXMOL .LE. 0 .OR. IDXMOL .GT. SIZE(IDMQAL) ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-HITQAL: Molec# from qualifier not in valid range,=' &
                 // C11INT(IDXMOL) 
        RETURN
      END IF
      IDMQAL(IDXMOL) = .TRUE.
      I = J+1
    END IF 
  END DO   
!
END SUBROUTINE HITQAL
END MODULE HITQAL_SUB
