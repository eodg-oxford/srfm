MODULE DRVSPC_SUB
CONTAINS
SUBROUTINE DRVSPC ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   22NOV23 AD Checked.
!   01MAY17 AD F90 conversion of inpspc.for. Checked.
!
! DESCRIPTION
!   Read RFM driver table *SPC section
!   Called once by RFMDRV.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE LENREC_DAT ! Max length of input text record
!
! SUBROUTINES
    USE NXTREC_SUB ! Load next record from input file
    USE SPCCHK_SUB ! Check spectral range and resolution data
    USE SPCFIL_SUB ! Load spectral range/resln from file
    USE SPCLAB_SUB ! Add Spc.range label to SPCCOM
    USE SPCTXT_SUB ! Read Spectral Range & Resln from text string
    USE TXTFLD_SUB ! Identify start and end points of text field in record
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for Driver File
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL  :: ANYSPC = .FALSE. ! T = at least one record in section
    LOGICAL           :: ENDSEC  ! Set TRUE when last record in section is found
    INTEGER(I4)       :: IEND    ! Pointer to end of field in text record
    INTEGER(I4)       :: ISTA    ! Pointer to start of field in text record
    INTEGER(I4)       :: NFLD    ! Counter for number of fields in record
    CHARACTER(LENREC) :: RECORD  ! Text record containing spec.range data
    CHARACTER(80)     :: WRNMSG  ! Warning message sent to log file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO                      ! Continue while more records in section
! Load next driver table record containing data
    CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( ENDSEC ) EXIT
    ANYSPC = .TRUE.
!
! Find how many fields in driver table record
    ISTA = 1
    NFLD = 0
    DO WHILE ( ISTA .NE. 0 )
      NFLD = NFLD + 1
      CALL TXTFLD ( RECORD, NFLD+1, ISTA, IEND )
    END DO
!
    SELECT CASE ( NFLD )    
    CASE ( 1 )                   ! 1 field: assume Spc.Range filename
      CALL SPCLAB ( '', FAIL, ERRMSG )     ! enter null label
      IF ( FAIL ) RETURN
      CALL TXTFLD ( RECORD, 1, ISTA, IEND ) ! Extract .spc filename
      CALL SPCFIL ( RECORD(ISTA:IEND), FAIL, ERRMSG )  
    CASE ( 2 )                   ! 2 fields: assume Label + Spc.range file
      CALL TXTFLD ( RECORD, 1, ISTA, IEND ) ! Extract Spc range label
      CALL SPCLAB ( RECORD(ISTA:IEND), FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      CALL TXTFLD ( RECORD, 2, ISTA, IEND ) ! Extract .spc filename
      CALL SPCFIL ( RECORD(ISTA:IEND), FAIL, ERRMSG )
    CASE ( 3 )                   ! 3 fields: assume unlabelled Spc.range
      CALL SPCLAB ( '', FAIL, ERRMSG )      ! enter null label
      IF ( FAIL ) RETURN
      CALL SPCTXT ( RECORD, FAIL, ERRMSG ) 
    CASE ( 4 )                   ! 4 fields: assume labelled Spc.range
      CALL TXTFLD ( RECORD, 1, ISTA, IEND )
      CALL SPCLAB ( RECORD(ISTA:IEND), FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      CALL SPCTXT ( RECORD(IEND+1:), FAIL, ERRMSG )
    CASE DEFAULT         ! 0 ot 5 or more fields: undefined structure
      ERRMSG = 'F-INPSPC: Unable to parse record: ' // RECORD(1:46)
      FAIL = .TRUE.
    END SELECT
    IF ( FAIL ) RETURN
!
    CALL SPCCHK ( FAIL, ERRMSG )
    IF ( FAIL ) RETURN
!
  END DO
!
! If no spectral range specified, note this in Log file (Warning status if
! spectral output also requested, ie ABS, BBT, OPT, RAD, RJT, TRA or WID flags).
  IF ( .NOT. ANYSPC ) THEN
    IF ( ABSFLG .OR. BBTFLG .OR. OPTFLG .OR. &
         RADFLG .OR. RJTFLG .OR. TRAFLG .OR. WIDFLG ) THEN
      WRNMSG = 'W-SPCCHK: No spectral calculations specified'
    ELSE
      WRNMSG = 'I-SPCCHK: No spectral calculations specified'
    END IF
    CALL WRTLOG ( WRNMSG )
    IF ( FAIL ) RETURN
  END IF
!
END SUBROUTINE DRVSPC
END MODULE DRVSPC_SUB
