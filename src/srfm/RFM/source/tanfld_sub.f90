MODULE TANFLD_SUB
CONTAINS
SUBROUTINE TANFLD ( FIELD, VALUE, FAIL, ERRMSG )
!
! VERSION
!   27JUL26 AD Checked.
!   01AUG25 AD Original. Extracted from TANCHK
!
! DESCRIPTION
!   Check if string is valid *TAN entry and insert in TANCOM
!   General purpose.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE TANCOM_DAT ! Tangent path data
    USE LENHGT_DAT ! Max length of height component of RFM output filenames
!
! SUBROUTINES
    USE HGTSTR_FNC ! Convert altitude to C*5 string
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: FIELD  ! Tangent string to be tested 
    REAL(R4),      INTENT(OUT) :: VALUE  ! Value read from FIELD
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)       :: IOS    ! Saved value of IOSTAT
    INTEGER(I4)       :: ITAN   ! Index of inserted value in TANCOM
    CHARACTER(LENHGT) :: TANSTR ! Tan.value written as a string 
    TYPE(TANTYP), ALLOCATABLE :: TANSAV(:) ! Saved TAN during reallocation
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! To be identified as a real number, the string FIELD must be readable without
! an error
  READ ( FIELD, *, IOSTAT=IOS ) VALUE
  IF ( IOS .NE. 0 ) THEN
    ERRMSG = 'F-TANFLD: Unreadable value in *TAN section: ' // FIELD
    FAIL = .TRUE.
    RETURN
  ELSE IF ( ABS ( VALUE ) .GE. 99999.5 ) THEN
    ERRMSG = 'F-TANFLD: Cannot handle values .GE. 99999.5'
    FAIL = .TRUE.
    RETURN
  END IF
!
! Convert value to string that will form part of output filenames and check
! that this is unique
  TANSTR = HGTSTR ( VALUE ) 
!
! Check that this is distinguishable from other USRTAN values, 
  IF ( ALLOCATED ( TAN ) ) THEN
    IF ( ANY ( TAN%STR .EQ. TANSTR ) ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-TANFLD: Repeated field=' // TANSTR
      RETURN
    END IF
    ITAN = COUNT ( TAN%USR .LT. VALUE ) + 1  ! insertion index
    CALL MOVE_ALLOC ( TAN, TANSAV ) 
    NTAN = NTAN + 1
    ALLOCATE ( TAN(NTAN) ) 
    IF ( ITAN .GT. 1 ) TAN(1:ITAN-1) = TANSAV(1:ITAN-1)
    IF ( ITAN .LT. NTAN ) TAN(ITAN+1:NTAN) = TANSAV(ITAN:)
  ELSE
    NTAN = 1
    ALLOCATE ( TAN(NTAN) ) 
    ITAN = 1
  END IF
!
  TAN(ITAN)%USR = VALUE
  TAN(ITAN)%STR = TANSTR
  TAN(ITAN)%CLC = .FALSE.
  TAN(ITAN)%IAT = 0   
  TAN(ITAN)%HGT = 0.0
  TAN(ITAN)%JDX = 0
  TAN(ITAN)%SKY = .FALSE.
  TAN(ITAN)%ISK = 0
  TAN(ITAN)%SFC = .FALSE.
  TAN(ITAN)%ITN = 0 
  TAN(ITAN)%SEC = 0.0
!
  MTAN = NTAN
!
END SUBROUTINE TANFLD
END MODULE TANFLD_SUB

