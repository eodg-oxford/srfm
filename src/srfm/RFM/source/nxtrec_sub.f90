MODULE NXTREC_SUB
CONTAINS
SUBROUTINE NXTREC ( LUN, RECORD, ENDSEC, FAIL, ERRMSG )
!
! VERSION
!   14DEC23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Load next record from input file
!   General purpose module.
! 
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUN    ! Logical Unit Number
    CHARACTER(*),  INTENT(OUT) :: RECORD ! Last record read from file
    LOGICAL,       INTENT(OUT) :: ENDSEC ! T=reached end of section
    LOGICAL,       INTENT(OUT) :: FAIL   ! T=Fatal error detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IEND ! Pointer to last part before end-of-line comments 
    INTEGER(I4) :: IOS  ! Saved value of IOSTAT for error message
    INTEGER(I4) :: IPT  ! Pointer to each character in record
    INTEGER(I4) :: LREC ! Length of input RECORD string
    CHARACTER(LEN(RECORD)+10) :: RECEXT ! Extended version of RECORD
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
  ENDSEC = .FALSE.
  LREC = LEN ( RECORD )
!
  DO
    READ ( LUN, '(A)', IOSTAT=IOS ) RECEXT
    IF ( IOS .GT. 0 ) THEN                   ! fatal read error
      FAIL = .TRUE.
      WRITE ( ERRMSG, * ) 'F-NXTREC: Error reading file. IOSTAT=', IOS
      EXIT
    ELSE IF ( IOS .LT. 0 ) THEN              ! end-of-file reached
      ENDSEC = .TRUE.
      EXIT
    ELSE IF ( RECEXT(1:1) .EQ. '*' ) THEN    ! end-of-section reached
      BACKSPACE ( LUN, IOSTAT=IOS ) 
      ENDSEC = .TRUE.
      EXIT
    ELSE IF ( RECEXT .EQ. '' ) THEN          ! empty record, or all blanks
      CYCLE
    ELSE IF ( INDEX ( ADJUSTL(RECEXT), '!' ) .EQ. 1 ) THEN ! comment record
      CYCLE      
    ELSE IF ( LREC .EQ. 0 ) THEN             ! position before next data record
      BACKSPACE ( LUN, IOSTAT=IOS )
      EXIT
    ELSE                                     ! return next data record
      IEND = INDEX ( RECEXT, '!' ) - 1
      IF ( IEND .EQ. -1 .OR. IEND .GT. LREC ) THEN
        DO IPT = LREC+1, LEN(RECEXT)         ! Check for extra chars
          IF ( RECEXT(IPT:IPT) .NE. ' ' .AND. &
               RECEXT(IPT:IPT) .NE. '!'       ) THEN
            FAIL = .TRUE.
            ERRMSG = 'F-NXTREC: Extra characters beyond expected '// &
                     'end of record:' // RECEXT(IPT:)
            EXIT
          END IF
        END DO
        IEND = LREC
      END IF
      RECORD = RECEXT(1:IEND)  
      EXIT
    END IF
  END DO
!
END SUBROUTINE NXTREC
END MODULE NXTREC_SUB
