MODULE PARFIL_SUB
CONTAINS
SUBROUTINE PARFIL ( LUNHIT, NAMHIT, WNOREQ, IDMREQ, USEFIL, IDMFIL, FAIL, ERRMSG )
!
! VERSION
!   22AUG24 AD Checked.
!   11AUG23 AD Changed arguments
!   19MAY23 AD Original. Based on part of OPNHIT.
!
! DESCRIPTION
!   Check HITRAN line data .par file
!   Called by HITFIL once for each .par file file
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNHIT    ! LUN for open HITRAN file
    CHARACTER(*),  INTENT(IN)  :: NAMHIT    ! Name of HITRAN file
    REAL(R8),      INTENT(IN)  :: WNOREQ(2) ! Min/Max reqd wno [cm-1]
    LOGICAL,       INTENT(IN)  :: IDMREQ(:) ! T=required molecule
    LOGICAL,       INTENT(OUT) :: USEFIL    ! T=file contains usable contents
    LOGICAL,       INTENT(OUT) :: IDMFIL(:) ! State of found molecules
    LOGICAL,       INTENT(OUT) :: FAIL      ! TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG    ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL     :: CHKIDM ! T=check molecules in file
    LOGICAL     :: CHKWNO ! T=check wavenumber range of file
    INTEGER(I4) :: IDM    ! HITRAN index of molecule
    INTEGER(I4) :: IOS    ! Saved value of IOSTAT for error messages
    REAL(R8)    :: WNO    ! Wavenumber [cm-1] of record 
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL   = .FALSE.
  CHKWNO = WNOREQ(2) .GT. 0.0D0
  CHKIDM = ANY ( IDMREQ ) 
  IF ( CHKIDM ) IDMFIL = .FALSE.
!
  OPEN ( UNIT=LUNHIT, FILE=NAMHIT, STATUS='OLD', ACTION='READ', &
         IOSTAT=IOS, ERR=900 )
!
! Read first record
  READ ( LUNHIT, '(I2,X,F12.6)', IOSTAT=IOS, END=800, ERR=900 ) IDM, WNO
  IF ( CHKWNO ) THEN 
    IF ( WNO .GT. WNOREQ(1) ) THEN
      CALL WRTLOG ( 'W-PARFIL: File ignored. Starts above reqd min wavenumber' )
      RETURN
    END IF
  END IF
  BACKSPACE ( LUNHIT ) 
! 
  DO WHILE ( CHKWNO .OR. CHKIDM ) 
    READ ( LUNHIT, '(I2,X,F12.6)', IOSTAT=IOS, END=800, ERR=900 ) IDM, WNO
    CHKWNO = WNO .LT. WNOREQ(2)
    IF ( CHKIDM ) THEN
      IDMFIL(IDM) = IDMREQ(IDM)
      CHKIDM = ANY ( IDMREQ .NEQV. IDMFIL )
    END IF
  END DO
  USEFIL = .TRUE.
  RETURN
!
800 CONTINUE
  IF ( WNO .LT. WNOREQ(2) ) THEN
    CALL WRTLOG ( 'W-PARFIL: File ignored. Ends below reqd max wavenumber')
  ELSE IF ( .NOT. ANY ( IDMFIL ) ) THEN
    CALL WRTLOG ( 'W-PARFIL: File ignored, contains no required molecules' )
  ELSE
    USEFIL = .TRUE.
  END IF
!
900 CONTINUE
  FAIL = IOS .GT. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-PARFIL: Read failure on HITRAN file. IOSTAT=', IOS
!
END SUBROUTINE PARFIL
END MODULE PARFIL_SUB
