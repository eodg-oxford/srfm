MODULE INIPAR_SUB
CONTAINS
SUBROUTINE INIPAR ( LUNHIT, WNOREQ, FAIL, ERRMSG ) 
!
! VERSION
!   21AUG24 AD Checked.
!   11AUG23 AD Original. Simplified from HITREC
!
! DESCRIPTION
!   Initialise pointer in HITRAN line data .par file
!   Called by INIHIT.   
!   This sets the pointer to read the next record with wno .ge. wnoreq
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNHIT ! Logical Unit Number
    REAL(R8),      INTENT(IN)  :: WNOREQ ! Required wavenumber
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IOS ! Saved value of IOSTAT for error messages
    REAL(R8)    :: WNO ! Wavenumber [cm-1]
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
!
! Ensure WNO is defined if currently at end-of-file
  WNO = WNOREQ * 2.0
!
! Only require wavenumber from HITRAN .par file records
  READ ( LUNHIT, '(3X,F12.6)', ERR=900, END=700, IOSTAT=IOS ) WNO
!
700 CONTINUE
  IF ( WNO .GT. WNOREQ ) THEN
    REWIND ( LUNHIT, ERR=900, IOSTAT=IOS )
    WNO = -1.0D0
  END IF
!
  DO WHILE ( WNO .LT. WNOREQ )
    READ ( LUNHIT, '(3X,F12.6)', ERR=900, END=800, IOSTAT=IOS ) WNO
  END DO
!
  BACKSPACE ( LUNHIT, ERR=900, IOSTAT=IOS )
  RETURN
!
! Not expected to reach end-of-file searching forwards for WNOREQ
800 CONTINUE
  STOP 'F-INIPAR: Logical error'
!
900 CONTINUE
  FAIL = .TRUE.
  WRITE ( ERRMSG, * ) &
    'F-INIPAR: Failed to read record in HITRAN File. IOSTAT=', IOS
!
END SUBROUTINE INIPAR
END MODULE INIPAR_SUB
