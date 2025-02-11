MODULE CHKNAM_SUB
CONTAINS
SUBROUTINE CHKNAM ( FAIL, ERRMSG )
!
! VERSION
!   25MAR23 AD Checked.
!   08NOV17 AD F90 original. Checked.
!
! DESCRIPTION
!   Check RFM output filename templates
!   Called by DRVCHK for any driver file sections defining output filenames.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE NAMCOM_DAT ! RFM output filenames
    USE GASCOM_DAT, ONLY: NGAS   ! No. different molecules
    USE GRACOM_DAT, ONLY: NPSI   ! No. different horiz.profile locations
    USE SPCCOM_DAT, ONLY: NSPC   ! No. different spectral ranges
    USE TANCOM_DAT, ONLY: NTAN   ! No. output tangent paths
!
! SUBROUTINES
    USE LOCASE_FNC ! Convert text string to lower case
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    INTEGER, PARAMETER :: NTYP = 11
    CHARACTER(3), PARAMETER :: TYPLST(NTYP) = & ! Different file types
      (/ 'ABS', 'BBT', 'COO', 'OPT', 'PRF', 'PTH', 'RAD', 'RJT', 'TAB', &
         'TRA', 'WID' /)
!
! LOCAL VARIABLES
    LOGICAL           :: LMULTI ! T = multiple output files for type required
    INTEGER(I4)       :: IPT    ! Pointer to '*' in supplied filename
    INTEGER(I4)       :: ITYP   ! Counter for different output file types
    CHARACTER(LENNAM) :: LOWNAM ! Lower case version of FILNAM
    CHARACTER(LENNAM), POINTER :: FILNAM ! Output file template
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  DO ITYP = 1, NTYP
    SELECT CASE ( TYPLST(ITYP) )
    CASE ( 'ABS' ) 
      IF ( .NOT. ABSFLG ) CYCLE
      FILNAM => ABSNAM
      LMULTI = NTAN .GT. 1 .OR. NSPC .GT. 1 .OR. JACFLG .OR. LEVFLG
    CASE ( 'BBT' ) 
      IF ( .NOT. BBTFLG ) CYCLE
      FILNAM => BBTNAM
      LMULTI = NTAN .GT. 1 .OR. NSPC .GT. 1 .OR. JACFLG .OR. LEVFLG
    CASE ( 'COO' ) 
      IF ( .NOT. COOFLG ) CYCLE
      FILNAM => COONAM
      LMULTI = NTAN .GT. 1 .OR. NSPC .GT. 1 .OR. JACFLG .OR. LEVFLG
    CASE ( 'OPT' ) 
      IF ( .NOT. OPTFLG ) CYCLE
      FILNAM => OPTNAM
      LMULTI = NTAN .GT. 1 .OR. NSPC .GT. 1 .OR. JACFLG .OR. LEVFLG
    CASE ( 'PRF' ) 
      IF ( .NOT. GRAFLG ) CYCLE
      FILNAM => PRFNAM
      LMULTI = NPSI .GT. 1 
    CASE ( 'PTH' ) 
      IF ( .NOT. PTHFLG ) CYCLE
      FILNAM => PTHNAM
      LMULTI = NTAN .GT. 1 .OR. NPSI .GT. 1
    CASE ( 'RAD' ) 
      IF ( .NOT. RADFLG ) CYCLE
      FILNAM => RADNAM
      LMULTI = NTAN .GT. 1 .OR. NSPC .GT. 1 .OR. JACFLG .OR. LEVFLG
    CASE ( 'RJT' ) 
      IF ( .NOT. RJTFLG ) CYCLE
      FILNAM => RJTNAM
      LMULTI = NTAN .GT. 1 .OR. NSPC .GT. 1 .OR. JACFLG .OR. LEVFLG
    CASE ( 'TAB' ) 
      IF ( .NOT. TABFLG ) CYCLE
      FILNAM => TABNAM
      LMULTI = NGAS .GT. 1 .OR. NSPC .GT. 1
    CASE ( 'TRA' ) 
      IF ( .NOT. TRAFLG ) CYCLE
      FILNAM => TRANAM
      LMULTI = NTAN .GT. 1 .OR. NSPC .GT. 1 .OR. JACFLG .OR. LEVFLG
    CASE ( 'WID' ) 
      IF ( .NOT. WIDFLG ) CYCLE
      FILNAM => WIDNAM
      LMULTI = NGAS .GT. 1 .OR. NSPC .GT. 1
    CASE DEFAULT
      STOP 'F-CHKNAM: Logical error'
    END SELECT
!
! If multiple output files required, check filename has a '*' symbol
    IF ( LMULTI .AND. INDEX ( FILNAM, '*') .EQ. 0 ) THEN
      FAIL = .TRUE. 
      ERRMSG = 'F-CHKNAM: ''*'' character required in ' // TYPLST(ITYP) // &
                ' filename template'
      RETURN
    END IF
!
    LOWNAM = LOCASE ( FILNAM ) 
    IF ( BINFLG ) THEN
      IPT = INDEX ( LOWNAM, '.asc' ) 
      IF ( IPT .NE. 0 ) FILNAM(IPT:IPT+3) = '.bin'
    ELSE
      IPT = INDEX ( LOWNAM, '.bin' ) 
      IF ( IPT .NE. 0 ) FILNAM(IPT:IPT+3) = '.asc'
    END IF
!
  END DO
!
  NULLIFY ( FILNAM )
!
END SUBROUTINE CHKNAM
END MODULE CHKNAM_SUB
