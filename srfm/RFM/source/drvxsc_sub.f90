MODULE DRVXSC_SUB
CONTAINS
SUBROUTINE DRVXSC ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   23NOV23 AD Checked.
!   01MAY17 AD F90 conversion of inpxsc.for. Tested.
!
! DESCRIPTION
!   Read RFM driver table *XSC section
!   Called by RFMDRV once
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE XSCCOM_DAT, ONLY: NXSC   ! No. of X/S datasets loaded
!
! SUBROUTINES
    USE NXTFLD_SUB ! Load next field from section of driver file
    USE XSCFIL_SUB ! Check if .xsc file required
    USE XSCDEF_SUB ! Use default .xsc filename to find any missing files
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for Driver File
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)       :: LENGTH      ! Length of field read from driver file
    CHARACTER(LENREC) :: NAMXSC      ! Filename read from *XSC section
    CHARACTER(LENREC) :: NAMDEF = '' ! Default filename template
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO                                  ! Read next field in driver table section
    CALL NXTFLD ( LUNDRV, NAMXSC, LENGTH, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
!
    IF ( INDEX ( NAMXSC(1:LENGTH), '*' ) .EQ. 0 ) THEN  ! explicit filename
      CALL XSCFIL ( NAMXSC(1:LENGTH), FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN 
    ELSE IF ( NAMDEF .EQ. '' ) THEN                     ! filename template
      NAMDEF = NAMXSC
    ELSE
      ERRMSG = 'F-DRVXSC: *XSC section has more than one filename template'
      FAIL = .TRUE.
      RETURN
    END IF
  END DO
!
  IF ( NAMDEF .NE. '' ) THEN
    CALL XSCDEF ( NAMDEF, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
  END IF
!
  IF ( NXSC .EQ. 0 ) &
    CALL WRTLOG ( 'W-DRVXSC: No useful filenames in *XSC section' )
!
END SUBROUTINE DRVXSC
END MODULE DRVXSC_SUB

