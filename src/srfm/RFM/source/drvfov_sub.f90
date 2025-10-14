MODULE DRVFOV_SUB
CONTAINS
SUBROUTINE DRVFOV ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   10FEB24 AD Checked.
!   01MAY17 AD F90 conversion of inpfov.for. Tested.
!
! DESCRIPTION
!   Read RFM driver table *FOV section
!   Called once by RFMDRV if FOV flag enabled.
!   Expects to find a single field in section: the name of the FOV file.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
!
! SUBROUTINES
    USE ENDCHK_SUB ! Check end of Driver Table section has been reached
    USE FOVFIL_SUB ! Read FOV data file
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
    INTEGER(I4)       :: LENGTH ! Length of name of FOV file
    CHARACTER(LENREC) :: FILFOV ! Name of FOV file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Read name of .fov file and open
  CALL NXTFLD ( LUNDRV, FILFOV, LENGTH, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
  IF ( LENGTH .EQ. 0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVFOV: no filename in *FOV section of driver file'
    RETURN
  END IF
!
! Read in data
  CALL FOVFIL ( FILFOV(1:LENGTH), FAIL, ERRMSG )   
  IF ( FAIL ) RETURN
!
! Check no more fields in this section
  CALL ENDCHK ( LUNDRV, '*FOV', FAIL, ERRMSG )
!
END SUBROUTINE DRVFOV
END MODULE DRVFOV_SUB
