MODULE HDBFIL_SUB
CONTAINS
SUBROUTINE HDBFIL ( LUNHIT, NAMHIT, WNOREQ, IDMREQ, USEFIL, IDMFIL, HDB, FAIL, ERRMSG )
!
! VERSION
!   09AUG24 AD Checked.
!   11AUG23 AD Return data via arguments
!   30MAY23 AD Original.
!
! DESCRIPTION
!   Check HITRAN line database file
!   Called by HITFIL
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE HDBCOM_DAT ! HDB data structure
!
! SUBROUTINES
    USE HDBCHK_SUB ! Check HITRAN database file header record
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNHIT    ! LUN for open HITRAN file
    CHARACTER(*),  INTENT(IN)  :: NAMHIT    ! Name of HITRAN file
    REAL(R8),      INTENT(IN)  :: WNOREQ(2) ! Lowest  reqd wno [cm-1]
    LOGICAL,       INTENT(IN)  :: IDMREQ(:) ! Flags for HITRAN molecules
    LOGICAL,       INTENT(OUT) :: USEFIL    ! T=file contains usable contents
    LOGICAL,     INTENT(INOUT) :: IDMFIL(:) ! T=use molec# from file
    TYPE(HDBTYP),  INTENT(OUT) :: HDB       ! Database file auxiliary data
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
  CALL HDBCHK ( LUNHIT, HDB, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN
!
  REWIND ( LUNHIT, ERR=900, IOSTAT=IOS )
  READ ( LUNHIT, *, IOSTAT=IOS, ERR=900 ) ! Skip header record  
  READ ( LUNHIT, HDB%WNOFMT, IOSTAT=IOS, ERR=900 ) WNO ! Wno of 1st rec
  IF ( CHKWNO ) THEN
    IF ( WNO .GT. WNOREQ(1) ) THEN
      CALL WRTLOG ( 'W-PARFIL: File ignored. Starts above reqd min wavenumber' )
      RETURN
    END IF
  END IF
  BACKSPACE ( LUNHIT, IOSTAT=IOS, ERR=900 ) 
!
! Read through file until end or all molecules found
  DO WHILE ( CHKIDM ) 
    READ ( LUNHIT, HDB%IDMFMT, IOSTAT=IOS, ERR=900, END=700 ) IDM
    IF ( IDM .LE. 0 .OR. IDM .GE. SIZE(IDMREQ) ) CYCLE
    IDMFIL(IDM) = IDMREQ(IDM)
    CHKIDM = .NOT. ( ALL ( IDMREQ .EQV. IDMFIL ) )    ! not All molecules found
  END DO
!
700 CONTINUE
  BACKSPACE ( LUNHIT, IOSTAT=IOS, ERR=900 ) ! Ready to read WNO if reqd
  IF ( IOS .LT. 0 ) &                       ! If EOF, need to go back another
    BACKSPACE ( LUNHIT, IOSTAT=IOS, ERR=900 ) ! Ready to read WNO if reqd
!
! Read through until end or required wavenumber range spanned
  DO WHILE ( CHKWNO ) 
    READ ( LUNHIT, HDB%WNOFMT, IOSTAT=IOS, ERR=900, END=800 ) WNO
    CHKWNO = WNO .LT. WNOREQ(2) 
  END DO
!
800 CONTINUE
!
  IF ( WNO .LT. WNOREQ(2) ) THEN
    CALL WRTLOG ( 'W-HDBFIL: File ignored. Ends below reqd max wavenumber')
    USEFIL = .FALSE.
  ELSE IF ( .NOT. ANY ( IDMFIL ) ) THEN
    CALL WRTLOG ( 'W-HDBFIL: File ignored, contains no required molecules' )
    USEFIL = .FALSE.
  ELSE
    USEFIL = .TRUE.
  END IF
  RETURN
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-HDBFIL: Error reading HITRAN database file. IOSTAT=', IOS
!
END SUBROUTINE HDBFIL
END MODULE HDBFIL_SUB
