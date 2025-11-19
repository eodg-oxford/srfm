MODULE BINFIL_SUB
CONTAINS
SUBROUTINE BINFIL ( LUNHIT, NAMHIT, WNOREQ, IDMREQ, USEFIL, IDMFIL, FAIL, ERRMSG )
!
! VERSION
!   05AUG24 AD Checked.
!   11AUG23 AD Return file values via arguments
!   19MAY23 AD Simplified from original OPNHIT
!
! DESCRIPTION
!   Check HITRAN binary data file 
!   Called by HITFIL once for each .bin data file
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE SETIFP_SUB ! Set forward pointers for HITRAN binary file
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNHIT    ! LUN for open HITRAN file
    CHARACTER(*),  INTENT(IN)  :: NAMHIT    ! Name of HITRAN file
    REAL(R8),      INTENT(IN)  :: WNOREQ(2) ! Lowest  reqd wno [cm-1]
    LOGICAL,       INTENT(IN)  :: IDMREQ(:) ! Highest reqd wno [cm-1]
    LOGICAL,       INTENT(OUT) :: USEFIL    ! T=file contains usable contents
    LOGICAL,       INTENT(OUT) :: IDMFIL(:) ! T=use molec# from file
    LOGICAL,       INTENT(OUT) :: FAIL      ! TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG    ! Error message if FAIL is TRUE
!
! GLOBAL CONSTANTS
    INTEGER(I4), PARAMETER :: PARLEN(22) = 0 ! Dummy array to find RECL
!
! LOCAL VARIABLES
    LOGICAL       :: CHKIDM ! T=check molecules in file
    LOGICAL       :: CHKWNO ! T=check wavenumber range of file
    INTEGER(I4)   :: IDUM   ! Dummy integer(i4) for read
    INTEGER(I4)   :: IFMT   ! Binary file format identifier
    INTEGER(I4)   :: IOS    ! Saved value of IOSTAT for error messages
    INTEGER(I4)   :: IREC1  ! Rec#1 of first line param record in file
    INTEGER(I4)   :: IREC2  ! Rec#2 of last line param record in file
    INTEGER(I4)   :: LSTAT  ! Used to identify HITRAN Fwd.Pointer block
    INTEGER(I4)   :: RECLEN ! RECL parameter for opening file (22 or 88)
    REAL(R8)      :: WNOLOW ! Lowest wavenumber [cm-1] record in file
    REAL(R8)      :: WNOUPP ! Highest wavenumber [cm-1] record in file
    CHARACTER(48) :: LABEL  ! Label read from HITRAN file
    INTEGER(I4)   :: IFP(SIZE(IDMFIL)) ! Forward pointers
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
  CHKWNO = WNOREQ(2) .GT. 0.0D0
  CHKIDM = ANY ( IDMREQ ) 
!
  INQUIRE ( IOLENGTH=RECLEN ) PARLEN
  CALL WRTLOG ( 'I-OPNBIN: Opening binary file with RECL=' // C11INT(RECLEN) )
  OPEN ( UNIT=LUNHIT, FILE=NAMHIT, STATUS='OLD', ACTION='READ', &
         ACCESS='DIRECT', RECL=RECLEN, IOSTAT=IOS, ERR=900 )
!
! Read header and extract first/last record#
  READ ( LUNHIT, REC=1, IOSTAT=IOS, ERR=900 ) &
    LSTAT, IFMT, IREC1, IREC2, LABEL
  CALL WRTLOG ( 'I-BINFIL: HITRAN File Label=' // LABEL )
!
! Find wavenumber range by reading first,last records
  READ ( LUNHIT, REC=IREC1, IOSTAT=IOS, ERR=900 ) LSTAT, IDUM, IDUM, WNOLOW
  READ ( LUNHIT, REC=IREC2, IOSTAT=IOS, ERR=900 ) LSTAT, IDUM, IDUM, WNOUPP 
!
  IF ( CHKWNO ) THEN
    IF ( WNOLOW .GT. WNOREQ(1) .OR. WNOUPP .LT. WNOREQ(2) ) THEN
      CALL WRTLOG ( 'W-BINFIL: File ignored - ' // &
                    'does not span required wavenumber range' )
      USEFIL = .FALSE.
      RETURN
    END IF
  END IF
!
! Check which molecules are in the file by setting forward pointers from 1st rec
  IF ( CHKIDM ) THEN
    CALL SETIFP ( LUNHIT, IREC1, IFP, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IDMFIL = IDMREQ .AND. IFP .LT. IREC2 ! Flag reqd molecules in the file
    USEFIL = ANY ( IDMFIL ) 
  ELSE
    USEFIL = .TRUE.
  END IF
!
900 CONTINUE
  FAIL = IOS .GT. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) & 
    'F-BINFIL: Read failure on HITRAN binary file. IOSTAT=', IOS
!
END SUBROUTINE BINFIL
END MODULE BINFIL_SUB
