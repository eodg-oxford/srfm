MODULE DRVOUT_SUB
CONTAINS
SUBROUTINE DRVOUT ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   12OCT20 AD Truncate TYPFIL from C4 to C3 for log file message. Checked.
!   08NOV17 AD Original. Checked.
!
! DESCRIPTION
!   Read RFM driver table *OUT section
!   Called by RFMDRV once if *OUT section found in driver table.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE NAMCOM_DAT ! RFM output filenames
!
! SUBROUTINES
    USE NXTREC_SUB ! Load next record from input file
    USE PARFLD_SUB ! Extract Parameter=Value string from record
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for Driver File
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: MAXPAR = 12 ! Max no. different PARAM types
!
! LOCAL VARIABLES
    LOGICAL           :: ENDSEC   ! T=end of *OUT section reached
    LOGICAL           :: GOTPAR   ! T=successfully extracted PARAM = VALUE 
    INTEGER(I4)       :: NPAR = 0 ! No of parameters stored so far
    CHARACTER(LENREC) :: FILNAM   ! Field read from driver table
    CHARACTER(LENREC) :: RECORD   ! Record read from driver table
    CHARACTER(6)      :: TYPFIL   ! Type of output file
    CHARACTER(6)      :: PARLST(MAXPAR) = '' ! List of PARAM values so far
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  DO             ! Continue while no more records in section
    CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( ENDSEC ) EXIT
    CALL PARFLD ( RECORD, GOTPAR, TYPFIL, FILNAM )
    IF ( .NOT. GOTPAR ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-DRVOUT: Record not of expected form PARAM = VALUE'
      RETURN
    END IF
    IF ( ANY ( PARLST .EQ. TYPFIL ) ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-DRVOUT: Repeated PARAM field: ' // TYPFIL
      RETURN
    ELSE
      NPAR = NPAR + 1
      PARLST(NPAR) = TYPFIL
    END IF      
    SELECT CASE ( TYPFIL )
      CASE ( 'OUTDIR' ) ; DIRNAM = FILNAM
      CASE ( 'ABSFIL' ) ; ABSNAM = FILNAM
      CASE ( 'BBTFIL' ) ; BBTNAM = FILNAM
      CASE ( 'COOFIL' ) ; COONAM = FILNAM
      CASE ( 'OPTFIL' ) ; OPTNAM = FILNAM
      CASE ( 'PRFFIL' ) ; PRFNAM = FILNAM
      CASE ( 'PTHFIL' ) ; PTHNAM = FILNAM
      CASE ( 'RADFIL' ) ; RADNAM = FILNAM
      CASE ( 'RJTFIL' ) ; RJTNAM = FILNAM
      CASE ( 'TABFIL' ) ; TABNAM = FILNAM
      CASE ( 'TRAFIL' ) ; TRANAM = FILNAM
      CASE ( 'WIDFIL' ) ; WIDNAM = FILNAM
      CASE DEFAULT
        FAIL = .TRUE.
        ERRMSG = 'F-DRVOUT: Unrecognised output file type: ' // TYPFIL
        RETURN
    END SELECT
    IF ( TYPFIL .EQ. 'OUTDIR' ) THEN
      CALL WRTLOG ( 'I-DRVOUT: Output directory: ' // FILNAM )
    ELSE
      CALL WRTLOG ( 'I-DRVOUT: User-supplied ' // TYPFIL(1:3) // ' filename: ' &
                  // FILNAM )
    END IF
  END DO
!
END SUBROUTINE DRVOUT
END MODULE DRVOUT_SUB
