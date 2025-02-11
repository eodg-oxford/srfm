MODULE INIHDB_SUB
CONTAINS
SUBROUTINE INIHDB ( LUNHIT, WNOREQ, WNOFMT, FAIL, ERRMSG ) 
!
! VERSION
!   19AUG24 AD Checked.
!   11AUG23 AD Pass WNOFMT as argument.
!   30MAY23 AD Original. 
!
! DESCRIPTION
!   Initialise pointer in HITRAN line database file
!   Called by INIHIT.
!   This sets the pointer to read the next record with wno .ge. wnoreq
!   In the interests of stepping quickly through the file, this just reads
!   the wavenumber, so only a suitable format statement is required, and
!   no further information on the record structure.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
!
! ARGUMENTS      
    INTEGER(I4),   INTENT(IN)  :: LUNHIT ! Logical unit number of HITRAN file
    REAL(R8),      INTENT(IN)  :: WNOREQ ! Required starting Wavenumber [cm-1]
    CHARACTER(*),  INTENT(IN)  :: WNOFMT ! Format string for reading Wavenumber
    LOGICAL,       INTENT(OUT) :: FAIL   ! TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message  if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IOS  ! Saved value of IOSTAT for error messages
    REAL(R8)    :: WNO  ! Wavenumber [cm-1] read from database file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
!
! Read the wavenumber of the next record in the file
  READ ( LUNHIT, WNOFMT, IOSTAT=IOS ) WNO
  IF ( IOS .GT. 0 ) GOTO 900
!
! If currently beyond required wavenumber or reached EOF, return to start
  IF ( WNO .GT. WNOREQ .OR. IOS .LT. 0 ) THEN
    REWIND ( LUNHIT, ERR=900, IOSTAT=IOS )
    READ ( LUNHIT, *, ERR=900, IOSTAT=IOS )   ! skip header record
    WNO = -1.0D0
  END IF
!
! Step forward until WNOREQ reached or exceeded
  DO WHILE ( WNO .LT. WNOREQ ) 
    READ ( LUNHIT, WNOFMT, IOSTAT=IOS, ERR=900 ) WNO
  END DO
!
! Reposition so last record will be read again 
  BACKSPACE ( LUNHIT, ERR=900, IOSTAT=IOS ) 
!
900 CONTINUE
  FAIL = IOS .GT. 0
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-INIHDB: Failed to read record in HITRAN database file. IOSTAT=', IOS
!
END SUBROUTINE INIHDB
END MODULE INIHDB_SUB

