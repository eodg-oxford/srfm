MODULE NXTFFL_SUB
CONTAINS
SUBROUTINE NXTFFL ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG, ONEREC )
!
! VERSION
!   18FEB23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Load next field from rfm.drv, expanding filenames
!   General purpose module.
!   This version of NXTFLD will first attempt to open any field as a filename
!   and, if successful, return subsequent fields from that file
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE LEXIST_FNC ! Check if file exists
    USE NXTFL2_SUB ! Load next field from secondary file
    USE NXTREC_SUB ! Load next record from input file
    USE OPNFIL_SUB ! Open input file
    USE TXTFLD_SUB ! Identify start and end points of text field in record
!
  IMPLICIT NONE
  SAVE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for Driver Table
    CHARACTER(*),  INTENT(OUT) :: FIELD  ! Field extracted from record
    INTEGER(I4),   INTENT(OUT) :: LENGTH ! Length of extracted field
    LOGICAL,       INTENT(OUT) :: FAIL   ! T=Fatal error detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
    LOGICAL, OPTIONAL, &
                   INTENT(IN)  :: ONEREC ! T=limit input to one record 
!
! LOCAL VARIABLES
    LOGICAL :: ENDSEC           ! T = end of section encountered
    LOGICAL :: FILINP = .FALSE. ! T = currently taking input from LUNTMP
    LOGICAL :: NEWREC = .TRUE.  ! T = new record to be read
    INTEGER(I4)       :: IEND   ! Pointer to end of field in record
    INTEGER(I4)       :: IFLD   ! Field counter within each record
    INTEGER(I4)       :: ISTA   ! Pointer to start of field in RECORD
    CHARACTER(LENREC) :: RECORD ! Record from Driver Table
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  FAIL = .FALSE.
  ENDSEC = .FALSE.
!
! Read next record from file if required
  DO 
    IF ( NEWREC ) THEN 
      CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      IF ( ENDSEC ) THEN                              ! Reached end of section
        LENGTH = 0
        RETURN
      END IF
      NEWREC = .FALSE.
      IFLD = 0
    END IF
!
! Read next field
    IF ( FILINP ) THEN             ! Currently reading data from file
      CALL NXTFL2 ( LUNTMP, FIELD, LENGTH, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      IF ( LENGTH .NE. 0 ) RETURN  ! Exit with next field from file
      CLOSE ( LUNTMP )             ! Reached end-of-file, go back to lundrv
      FILINP = .FALSE.
    END IF
!
! Read next field from LUNDRV
    IFLD = IFLD + 1
    CALL TXTFLD ( RECORD, IFLD, ISTA, IEND )
!
! If no more fields, either return (ONEREC) or read next record
    IF ( ISTA .EQ. 0 ) THEN        ! no (more) fields in record - get next rec
      NEWREC = .TRUE.              ! new record from LUNDRV required
      LENGTH = 0                   ! flag for no current FIELD loaded
      IF ( PRESENT ( ONEREC ) ) THEN
        IF ( ONEREC ) RETURN
      END IF
      CYCLE
    END IF
!
! Identified next field in LUNDRV record
    FIELD = RECORD(ISTA:IEND)
!
! If this is a valid filename, open it, set FILINP=TRUE and back to top to
! read first field from the file
    IF ( LEXIST ( FIELD ) ) THEN
      CALL OPNFIL ( LUNTMP, FIELD, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      FILINP = .TRUE. 
      CYCLE
    END IF    
!
    EXIT
  END DO
!
! If not a filename, return FIELD and its LENGTH
  LENGTH = IEND - ISTA + 1
!
END SUBROUTINE NXTFFL
END MODULE NXTFFL_SUB
