MODULE OPNFIL_SUB
CONTAINS
SUBROUTINE OPNFIL ( LUNFIL, NAMFIL, FAIL, ERRMSG )
!
! VERSION
!   10DEC23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Open input file
!   General purpose module.
!   Open ASCII file for READ, and skip/Log any initial comments.
!   First record written to Log file if it begins with '!'
!   All subsequent records beginning with '!' are skipped
!   Sets file pointer ready to read 1st record not beginning with '!'.
!
! VARAIBLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
!
! SUBROUTINES
    USE LEXIST_FNC ! Check if file exists
    USE NXTREC_SUB ! Load next record from input file
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE 
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNFIL ! LUN for opening/reading file
    CHARACTER(*),  INTENT(IN)  :: NAMFIL ! File name
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL           :: ENDSEC ! Set TRUE if end-of-section/file reached 
    INTEGER(I4)       :: IOS    ! Saved value of IOSTAT for error messages
    CHARACTER(LENREC) :: RECORD ! Text record read from file
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  CALL WRTLOG ( 'I-OPNFIL: Opening file: '//NAMFIL )
!
  IF ( .NOT. LEXIST ( NAMFIL ) ) THEN
    FAIL = .TRUE.
    IF ( LEN ( NAMFIL ) .GT. 65 ) THEN
      ERRMSG = 'F-OPNFIL: file not found: ' // NAMFIL(1:62) // '...'
    ELSE
      ERRMSG = 'F-OPNFIL: file not found: ' // TRIM ( NAMFIL ) 
    END IF
    RETURN
  END IF
!
  OPEN ( UNIT=LUNFIL, FILE=NAMFIL, STATUS='OLD', ACTION='READ', IOSTAT=IOS )
  IF ( IOS .NE. 0 ) THEN
    FAIL = .TRUE.
    WRITE ( ERRMSG, * ) 'F-OPNFIL: Failed to open file. IOSTAT=', IOS
    RETURN
  END IF 
!
! Print first record of file to Log file (should be a comment anyway)
  READ ( LUNFIL, '(A)', IOSTAT=IOS ) RECORD
  IF ( IOS .NE. 0 ) THEN
    FAIL = .TRUE.
    WRITE ( ERRMSG, * ) 'F-OPNFIL: Failed to read file. IOSTAT=', IOS
    RETURN
  END IF   
  CALL WRTLOG ( '  '//RECORD )
  BACKSPACE ( LUNFIL ) 
!
! Set pointer to first non-comment, non-blank record
  RECORD = ''
  CALL NXTREC ( LUNFIL, RECORD, ENDSEC, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
  IF ( .NOT. ENDSEC ) BACKSPACE ( LUNFIL ) 
!
END SUBROUTINE OPNFIL
END MODULE OPNFIL_SUB
