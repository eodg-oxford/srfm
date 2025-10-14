MODULE NXTFLD_SUB
CONTAINS
SUBROUTINE NXTFLD ( LUN, FIELD, LENGTH, FAIL, ERRMSG )
!
! VERSION
!   22FEB23 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Load next field from section of driver file
!   General purpose module.
!   If NEWREC = TRUE
!     Read next record
!     Read first field of record
!   Else
!     Read next field of current record
!   End if
!   If next record begins with section marker ('*'), backspaces so record
!   can be read again and sets LENGTH=0 for output, and NEWREC=TRUE for next
!   call.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
!
! SUBROUTINES
    USE NXTREC_SUB ! Load next record from input file
    USE TXTFLD_SUB ! Identify start and end points of text field in record
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUN    ! Logical Unit Number
    CHARACTER(*),  INTENT(OUT) :: FIELD  ! Field extracted from record
    INTEGER(I4),   INTENT(OUT) :: LENGTH ! Length of extracted field
    LOGICAL,       INTENT(OUT) :: FAIL   ! T=Fatal error detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL     :: ENDSEC             ! TRUE if end of section encountered
    LOGICAL     :: NEWREC = .TRUE.    ! TRUE if new record to be read
    INTEGER(I4) :: IEND               ! Pointer to end of field in record
    INTEGER(I4) :: IFLD = 0           ! Field counter within each record
    INTEGER(I4) :: ISTA               ! Pointer to start of field in RECORD
    CHARACTER(LENREC), SAVE :: RECORD ! Record from Driver Table
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  FAIL = .FALSE.
!
  DO
    IF ( NEWREC ) THEN ! Read next record from file if required
      CALL NXTREC ( LUN, RECORD, ENDSEC, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      IF ( ENDSEC ) THEN                              ! Reached end of section
        LENGTH = 0
        EXIT
      END IF
      NEWREC = .FALSE.
      IFLD = 0
    END IF
!
    IFLD = IFLD + 1    ! Read next field
    CALL TXTFLD ( RECORD, IFLD, ISTA, IEND )
    IF ( ISTA .GT. 0 ) THEN
      FIELD = RECORD(ISTA:IEND)
      LENGTH = IEND - ISTA + 1
      EXIT
    END IF
    NEWREC = .TRUE.
!
  END DO
!
END SUBROUTINE NXTFLD
END MODULE NXTFLD_SUB

