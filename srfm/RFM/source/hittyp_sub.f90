MODULE HITTYP_SUB
CONTAINS
SUBROUTINE HITTYP ( NAMHIT, TYPHIT, FAIL, ERRMSG )
!
! VERSION
!   15AUG24 AD Checked.
!   11AUG23 AD Remove LUNHIT argument, open with LUNTMP
!              Allow for qualifiers appended to filename
!   30MAY23 AD Simplified to three types. Add LUNHIT, TYPHIT arguments.
!   25APR20 AD Correction to reading data from .par file. Checked.
!   29JAN20 AD Original.
!
! DESCRIPTION
!   Identify type of HITRAN data file
!   Called by DRVHIT
!   This is for the old-format *HIT section where just the filename
!   is specified instead of PARAM=VALUE structure.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE LEXIST_FNC ! Check if file exists
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMHIT ! Name of HITRAN file
    CHARACTER(6),  INTENT(OUT) :: TYPHIT ! 'BINFIL', 'PARFIL' or 'HDBFIL'
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL CONSTANTS
    INTEGER(I4),  PARAMETER :: I2DUM(2) = 0 ! Dummy array to find RECL
!
! LOCAL VARIABLES
    INTEGER(I4) :: IBRACK ! Position of '(' character (if any) in NAMHIT
    INTEGER(I4) :: IFMT   ! Binary file format identifier
    INTEGER(I4) :: IOS    ! Saved value of IOSTAT for error messages
    INTEGER(I4) :: ISO    ! Isotope ID 
    INTEGER(I4) :: L      ! Length of NAMHIT excluding '(' character onwards
    INTEGER(I4) :: LSTAT  ! Used to identify HITRAN Fwd.Pointer block
    INTEGER(I4) :: MOL    ! Molecule ID
    INTEGER(I4) :: RECLEN ! RECL parameter for opening file (2 or 8)
    REAL(R8)    :: WNO    ! Wavenumber [cm-1]
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
!
  IBRACK = INDEX ( NAMHIT, '(' ) 
  IF ( IBRACK .GT. 0 ) THEN
    L = IBRACK - 1
  ELSE 
    L = LEN_TRIM ( NAMHIT ) 
  END IF
!
! Send message to LOG file saying which file is about to be opened
  CALL WRTLOG ( 'I-HITTYP: Opening HITRAN File: ' // NAMHIT(1:L) )
!
  IF ( .NOT. LEXIST ( NAMHIT(1:L) ) ) THEN 
    FAIL = .TRUE.
    ERRMSG = 'F-HITTYP: file not found'
    RETURN
  END IF
!
! Initially just require first two I4 integers from first record
  INQUIRE ( IOLENGTH=RECLEN ) I2DUM
!
! Start by assuming it is a direct-access (binary) file
  OPEN ( UNIT=LUNTMP, FILE=NAMHIT(1:L), STATUS='OLD', ACTION='READ', &
         ACCESS='DIRECT', RECL=RECLEN, IOSTAT=IOS, ERR=900 )
!
! Read header and extract format
  READ ( LUNTMP, REC=1, IOSTAT=IOS, ERR=900 ) LSTAT, IFMT
  CLOSE ( LUNTMP, IOSTAT=IOS, ERR=900 ) 
! 
! old version of HITBIN had IFMT=56  
  IF ( IFMT .EQ. 1 .OR. IFMT .EQ. 56 ) THEN 
    TYPHIT = 'BINFIL'
    CALL WRTLOG ( 'I-HITTYP: Identified as HITRAN binary file' )
!
! Next try opening as ASCII file
  ELSE 
    OPEN ( UNIT=LUNTMP, FILE=NAMHIT(1:L), STATUS='OLD', ACTION='READ', &
           IOSTAT=IOS, ERR=900 )
! 
    READ ( LUNTMP, '(I2,Z1,F12.6)', IOSTAT=IOS ) MOL, ISO, WNO
    CLOSE ( LUNTMP ) 
    IF ( IOS .EQ. 0 .AND. MOL .GT. 0 .AND. WNO .GT. 0.0D0 ) THEN
      TYPHIT = 'PARFIL'
      CALL WRTLOG ( 'I-HITTYP: Identified as HITRAN 160-character .par file' )
! else assume HDB file
    ELSE
      TYPHIT = 'HDBFIL'
      CALL WRTLOG ( 'I-HITTYP: Assuming a HITRAN database file' )
    END IF
  END IF
  RETURN
! 
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-HITTYP: Read failure on HITRAN file. IOSTAT=', IOS
!
END SUBROUTINE HITTYP
END MODULE HITTYP_SUB
