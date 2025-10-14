MODULE DRVGRD_SUB
CONTAINS
SUBROUTINE DRVGRD ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   20NOV23 AD Checked.
!   01MAY17 AD F90 conversion of inpgrd.for. Checked.
!
! DESCRIPTION
!   Read RFM driver table *GRD section
!   Called by RFMDRV once if *GRD key listed in RFM Driver table and GRDFLG=T
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE GFLCOM_DAT, ONLY: NGFL   ! No. .grd files stored in GFLCOM
!
! SUBROUTINES
    USE GRDDEF_SUB ! Use default .grd filename to find any missing files
    USE GRDFIL_SUB ! Open .grd file and check contents
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
    INTEGER(I4)       :: MGFL        ! Value of NGFL on input
    CHARACTER(LENREC) :: NAMDEF = '' ! GRD filename template
    CHARACTER(LENREC) :: NAMGRD      ! Name of GRD file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Irreg .spc files could already be stored in GFLCOM from *SPC section
  MGFL = NGFL
!
  DO                                 ! Read next field in driver table section
    CALL NXTFLD ( LUNDRV, NAMGRD, LENGTH, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
!
    IF ( INDEX ( NAMGRD(1:LENGTH), '*' ) .EQ. 0 ) THEN  ! explicit filename
      CALL GRDFIL ( NAMGRD(1:LENGTH), FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    ELSE IF ( NAMDEF .EQ. '' ) THEN                     ! filename template
      NAMDEF = NAMGRD
    ELSE
      ERRMSG = 'F-DRVGRD: *GRD section has more than one filename template'
      FAIL = .TRUE.
      RETURN
    END IF
  END DO
!
  IF ( NAMDEF .NE. '' ) THEN
    CALL GRDDEF ( NAMDEF, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
  END IF
!
  IF ( NGFL .EQ. MGFL ) &
    CALL WRTLOG ( 'W-DRVGRD: No useful filenames in *GRD section' ) 
!
END SUBROUTINE DRVGRD
END MODULE DRVGRD_SUB

