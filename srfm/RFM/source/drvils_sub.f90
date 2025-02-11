MODULE DRVILS_SUB
CONTAINS
SUBROUTINE DRVILS ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   11NOV23 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Read RFM driver table *ILS section
!   Called once by RFMDRV.
!   Expects to read the names of one or more instrument lineshape files.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
!
! SUBROUTINES
    USE ILSFIL_SUB ! Read ILS data file
    USE ILSSPC_SUB ! Assign ILS Functions to each spectral range
    USE NXTFLD_SUB ! Load next field from section of driver file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for Driver File
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL           :: ANYFIL = .FALSE. ! T=at least one filename in section
    INTEGER(I4)       :: LENGTH  ! Length of name of ILS file
    CHARACTER(LENREC) :: FILILS  ! Name of ILS file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO                   ! Loop over each ils filename in section
    CALL NXTFLD ( LUNDRV, FILILS, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
!
    ANYFIL = .TRUE.
    CALL ILSFIL ( FILILS(1:LENGTH), FAIL, ERRMSG )   
    IF ( FAIL ) RETURN
  END DO
!
  IF ( .NOT. ANYFIL ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVILS: no filename in *ILS section of driver file'
    RETURN
  END IF
!
! If no more files, assign ILS files for each spectral range
  CALL ILSSPC ( FAIL, ERRMSG )
!
END SUBROUTINE DRVILS
END MODULE DRVILS_SUB

