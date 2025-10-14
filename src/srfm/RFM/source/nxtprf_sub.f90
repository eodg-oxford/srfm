MODULE NXTPRF_SUB
CONTAINS
SUBROUTINE NXTPRF ( LUN, LABEL, PRF, FAIL, ERRMSG )
!
! VERSION
!   16FEB24 AD Bug#44 - test for NaN
!   28NOV23 AD Checked.
!   01MAY17 AD Original. Checked.
!
! DESCRIPTION
!   Get next profile from .atm file
!   Called by ATMFIL.
!   Assumes that PRF has correct profile size on input.
!   Returns profile label converted to upper case
!   PRF filled unless LABEL='*END' is returned
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
    INTEGER(I4),   INTENT(IN)  :: LUN    ! LUN for .atm format file
    CHARACTER(*),  INTENT(OUT) :: LABEL  ! Upper case label of profile (eg *TEM)
    REAL(R4),      INTENT(OUT) :: PRF(:) ! Profile of values
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: I      ! Position of last character in *label
    INTEGER(I4)   :: IOS    ! Saved value of IOSTAT
    CHARACTER(80) :: RECORD ! Record read from .atm file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Skip any comment or empty records before next *label record
  DO
    READ ( LUN, '(A)', IOSTAT=IOS, END=800 ) RECORD
    IF ( IOS .NE. 0 ) THEN
      FAIL = .TRUE.
      WRITE ( ERRMSG, * ) &
        'F-NXTPRF: Error reading .atm file record, IOSTAT=', IOS
      RETURN
    END IF
    IF ( RECORD(1:1) .NE. '!' .AND. RECORD .NE. '' ) EXIT
  END DO
!
  IF ( RECORD(1:1) .NE. '*' ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-NXTPRF: unexpected record in .atm file:'//RECORD(1:39)
    RETURN
  END IF
!
  I = INDEX ( RECORD//' ', ' ' ) - 1
  IF ( I .GT. LEN(LABEL) ) THEN
    FAIL = .TRUE.       
    WRITE ( ERRMSG, * ) 'F-NXTPRF: Label longer than expected:'//RECORD(1:43)
    RETURN
  END IF
  LABEL = UPCASE ( RECORD(1:I) ) 
! 
  IF ( LABEL .NE. '*END' ) THEN
    READ ( LUN, *, IOSTAT=IOS ) PRF
    IF ( IOS .NE. 0 ) THEN
      FAIL = .TRUE.
      WRITE ( ERRMSG, * ) 'F-NXTPRF: Error reading profile following label:' &
        // TRIM(LABEL) // '. IOSTAT=', IOS
      RETURN
    ELSE IF ( ANY ( ISNAN ( PRF ) ) ) THEN     ! Bug#44
      FAIL = .TRUE.
      WRITE ( ERRMSG, * ) 'F-NXTPRF: ' // TRIM(LABEL) // &
        ' profile contains a NaN value'
      RETURN
    END IF
  END IF
  RETURN
!
! Allow for any improperly terminated files
800 CONTINUE
  LABEL = '*END'
!
END SUBROUTINE NXTPRF
END MODULE NXTPRF_SUB

