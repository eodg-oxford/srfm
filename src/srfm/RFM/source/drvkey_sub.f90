MODULE DRVKEY_SUB
!
! Persistent bookkeeping for driver section keys.
!
  USE KIND_DAT
!
  IMPLICIT NONE
  PRIVATE
!
  INTEGER(I4),  PARAMETER :: MAXKEY = 30 ! >Max no different *KEY values
  INTEGER(I4),  PARAMETER :: NREQ   = 6  ! No. of reqd *KEY values
  CHARACTER(4), PARAMETER :: REQLST(NREQ) = &
      (/ '*HDR', '*FLG', '*SPC', '*GAS', '*ATM', '*TAN' /)
!
  INTEGER(I4) :: KEY_COUNT = 0
  CHARACTER(4) :: KEY_HISTORY(MAXKEY) = ''
!
  PUBLIC :: DRVKEY
  PUBLIC :: DRVKEY_RESET
!
CONTAINS
SUBROUTINE DRVKEY ( LUNDRV, KEY, FAIL, ERRMSG ) 
!
! VERSION
!   04JUL24 AD Change *LEV to *HGT when used for *TAN with flux calcs.
!   24JUN19 AD Remove CIA flag. Checked.
!   19DEC17 AD F90 version. Checked.
!
! DESCRIPTION
!   Check section key from driver table.
!   Called by RFMDRV for each driver table section.
!   DRVKEY_RESET clears the cached sequence so the validator can be reused.
!   First 6 keys are mandatory and checked for correct sequence, others are
!   checked against flags. KEY set to key value (eg '*HDR') if useful, or
!   to 'skip' if section not required. 
!   On final call, with LUNDRV = -1, list of keys found compared with list of
!   required keys.
!
! VARIABLE KINDS
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE LENREC_DAT ! Max length of input text record
    USE TANCOM_DAT, ONLY: USRELE, USRGEO ! User-specified elev. or geom.tan.hts
!
! SUBROUTINES
    USE UPCASE_FNC ! Convert text string to upper case
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for driver file, or -1 for check
    CHARACTER(4),  INTENT(OUT) :: KEY    ! Next key (eg '*HDR') 
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
! LOCAL VARIABLES
    LOGICAL           :: USEKEY   ! T=Use this section, F=ignore
    INTEGER(I4)       :: IOS      ! Saved value of IOSTAT
    CHARACTER(LENREC) :: RECORD   ! Text record read from driver file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( LUNDRV .EQ. -1 ) THEN   ! just check all expected sections found
    FAIL = .TRUE.
    IF ( FINFLG .AND. .NOT. ANY ( KEY_HISTORY .EQ. '*FIN' ) ) THEN
      ERRMSG = 'F-DRVKEY: FIN flag requires *FIN section in driver file'
    ELSE IF ( FOVFLG .AND. .NOT. ANY ( KEY_HISTORY .EQ. '*FOV' ) ) THEN
      ERRMSG = 'F-DRVKEY: FOV flag requires *FOV section in driver file'
    ELSE IF ( GRDFLG .AND. .NOT. ANY ( KEY_HISTORY .EQ. '*GRD' .OR. &
                                       KEY_HISTORY .EQ. '*LUT' ) ) THEN
      ERRMSG = 'F-DRVKEY: FOV flag requires *FOV section in driver file'
    ELSE IF ( ILSFLG .AND. .NOT. ANY ( KEY_HISTORY .EQ. '*ILS' ) ) THEN
      ERRMSG = 'F-DRVKEY: ILS flag requires *ILS section in driver file'
    ELSE IF ( JACFLG .AND. .NOT. ANY ( KEY_HISTORY .EQ. '*JAC' ) ) THEN
      ERRMSG = 'F-DRVKEY: JAC flag requires *JAC section in driver file'
    ELSE IF ( LEVFLG .AND. .NOT. ANY ( KEY_HISTORY .EQ. '*LEV' ) ) THEN
      ERRMSG = 'F-DRVKEY: LEV flag requires *LEV section in driver file'
    ELSE IF ( LUTFLG .AND. .NOT. ANY ( KEY_HISTORY .EQ. '*LUT' ) ) THEN
      ERRMSG = 'F-DRVKEY: LUT flag requires *LUT section in driver file'
    ELSE IF ( OBSFLG .AND. .NOT. ANY ( KEY_HISTORY .EQ. '*OBS' ) ) THEN
      ERRMSG = 'F-DRVKEY: OBS flag requires *OBS section in driver file'
    ELSE IF ( REJFLG .AND. .NOT. ANY ( KEY_HISTORY .EQ. '*REJ' ) ) THEN
      ERRMSG = 'F-DRVKEY: REJ flag requires *REJ section in driver file'
    ELSE IF ( SHPFLG .AND. .NOT. ANY ( KEY_HISTORY .EQ. '*SHP' ) ) THEN
      ERRMSG = 'F-DRVKEY: SHP flag requires *SHP section in driver file'
    ELSE
      FAIL = .FALSE.
    END IF
    CALL DRVKEY_RESET()
    RETURN
  END IF    
!
  READ ( LUNDRV, '(A80)', IOSTAT=IOS ) RECORD
  IF ( IOS .NE. 0 ) THEN
    FAIL = .TRUE.
    WRITE ( ERRMSG, * ) &
      'F-DRVKEY: Error reading next driver table section *KEY. IOSTAT=', IOS
    RETURN
  END IF
!
  KEY_COUNT = KEY_COUNT + 1
  KEY = UPCASE ( RECORD(1:4) )
!
! Convert valid alternatives to *TAN key
  IF ( KEY_COUNT .EQ. 6 ) THEN 
    IF ( TABFLG ) THEN
      IF ( KEY .EQ. '*DIM' ) KEY = '*TAN'          ! temporarily
    ELSE IF ( HOMFLG ) THEN
      IF ( KEY .EQ. '*LEN' ) KEY = '*TAN' 
    ELSE IF ( FLXFLG ) THEN
      IF ( KEY .EQ. '*LEV' .OR. &                  ! old value, or
           KEY .EQ. '*HGT'      ) KEY = '*TAN'     ! new value
    ELSE IF ( NADFLG .OR. ZENFLG ) THEN            ! & not FLXFLG 
      IF ( KEY .EQ. '*SEC' ) KEY = '*TAN'
      IF ( KEY .EQ. '*ELE' ) THEN
        USRELE = .TRUE.
        KEY = '*TAN'
      END IF
    ELSE                                           ! limb-viewing assumed
      IF ( KEY .EQ. '*ELE' ) THEN
        USRELE = .TRUE.
        KEY = '*TAN'
      ELSE IF ( KEY .EQ. '*GEO' ) THEN
        USRGEO = .TRUE.
        KEY = '*TAN'
      END IF
    END IF
    IF ( USRELE .AND. .NOT. OBSFLG ) THEN 
      FAIL = .TRUE.
      ERRMSG = 'F-DRVKEY: *ELE section only valid with OBS flag'
      RETURN
    END IF
  END IF   
!
  IF ( KEY_COUNT .LE. NREQ ) THEN         ! mandatory sections
    IF ( KEY .NE. REQLST(KEY_COUNT) ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-DRVKEY: Expected ' // REQLST(KEY_COUNT) // &
               ' but found ' // KEY // ' section in driver file.'
      RETURN
    END IF
  END IF      
!
! If creating .tab output files, use *DIM for section header
  IF ( KEY_COUNT .EQ. 6 .AND. TABFLG ) KEY = '*DIM'
!
! Check for repeated key
  IF ( ANY ( KEY_HISTORY .EQ. KEY ) ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVKEY: Duplicated section header: '//KEY
    RETURN
  END IF
  IF ( KEY_COUNT .GT. MAXKEY ) STOP 'F-DRVKEY: logical error'
  KEY_HISTORY(KEY_COUNT) = KEY 
!
! Set KEY to 'skip' for any sections which are not required
  SELECT CASE ( KEY ) 
  CASE ( '*ABS' ) ; USEKEY = ABSFLG
  CASE ( '*BBT' ) ; USEKEY = BBTFLG
  CASE ( '*COO' ) ; USEKEY = COOFLG
  CASE ( '*FIN' ) ; USEKEY = FINFLG
  CASE ( '*FOV' ) ; USEKEY = FOVFLG
  CASE ( '*GRD' ) ; USEKEY = GRDFLG
  CASE ( '*ILS' ) ; USEKEY = ILSFLG
  CASE ( '*JAC' ) ; USEKEY = JACFLG
  CASE ( '*LEV' ) ; USEKEY = LEVFLG
  CASE ( '*LUT' ) ; USEKEY = LUTFLG
  CASE ( '*NTE' ) ; USEKEY = NTEFLG
  CASE ( '*OBS' ) ; USEKEY = OBSFLG
  CASE ( '*OPT' ) ; USEKEY = OPTFLG
  CASE ( '*PRF' ) ; USEKEY = PRFFLG
  CASE ( '*PTH' ) ; USEKEY = PTHFLG
  CASE ( '*RAD' ) ; USEKEY = RADFLG
  CASE ( '*REJ' ) ; USEKEY = REJFLG
  CASE ( '*RJT' ) ; USEKEY = RJTFLG
  CASE ( '*SFC' ) ; USEKEY = SFCFLG
  CASE ( '*SHP' ) ; USEKEY = SHPFLG
  CASE ( '*SVD' ) ; USEKEY = SVDFLG
  CASE ( '*TAB' ) ; USEKEY = TABFLG
  CASE ( '*TRA' ) ; USEKEY = TRAFLG
  CASE ( '*WID' ) ; USEKEY = WIDFLG
  CASE DEFAULT ;    USEKEY = .TRUE.
  END SELECT
  IF ( .NOT. USEKEY ) KEY = 'skip'
!
! Normal exit - print section header to Log file
  CALL WRTLOG ( RECORD )
  FAIL = .FALSE.
!
END SUBROUTINE DRVKEY

SUBROUTINE DRVKEY_RESET()
  KEY_COUNT   = 0
  KEY_HISTORY = ''
END SUBROUTINE DRVKEY_RESET

END MODULE DRVKEY_SUB
