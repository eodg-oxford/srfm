MODULE DRVCIA_SUB
CONTAINS
SUBROUTINE DRVCIA ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   23JUL24 AD Checked.
!   01MAY17 AD F90 conversion of inpcia.for. Tested.
!
! DESCRIPTION
!   Read RFM driver table *CIA section
!   Called by RFMDRV.
!   Read Collision-Induced Abs. filenames from RFM driver table.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE CIACOM_DAT, ONLY: NCIA   ! No. CIA datasets loaded
!
! SUBROUTINES
    USE CIADEF_SUB ! Use default CIA filename to find any missing files
    USE CIAFIL_SUB ! Read Collision-Induced Abs. data file
    USE NXTFLD_SUB ! Load next field from section of driver file
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
    INTEGER(I4)       :: LENGTH      ! Length of field read from driver file
    CHARACTER(LENREC) :: NAMDEF = '' ! CIA filename template
    CHARACTER(LENREC) :: NAMCIA      ! Name of CIA file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO              !    Read next field in driver table
    CALL NXTFLD ( LUNDRV, NAMCIA, LENGTH, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
!
    IF ( INDEX ( NAMCIA(1:LENGTH), '*' ) .EQ. 0 ) THEN ! explicit filename
      CALL CIAFIL ( NAMCIA(1:LENGTH), FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    ELSE IF ( NAMDEF .EQ. '' ) THEN                    ! filename template
      NAMDEF = NAMCIA
    ELSE 
      ERRMSG = 'F-DRVCIA: *CIA section contains more than one filename template'
      FAIL = .TRUE.
      RETURN
    END IF
  END DO
!
  IF ( NAMDEF .NE. '' ) THEN
    CALL CIADEF ( NAMDEF, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
  END IF
!
  IF ( NCIA .EQ. 0 ) &
    CALL WRTLOG ( 'W-DRVCIA: No useful filenames in *CIA section' )
!
END SUBROUTINE DRVCIA
END MODULE DRVCIA_SUB
