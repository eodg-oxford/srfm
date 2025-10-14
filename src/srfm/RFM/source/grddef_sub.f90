MODULE GRDDEF_SUB
CONTAINS
SUBROUTINE GRDDEF ( NAMDEF, FAIL, ERRMSG )
!
! VERSION
!   13NOV23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Use default .grd filename to find any missing files
!   Called once by DRVGRD if filename template found in *GRD section of 
!   driver table.
!
! VARIABLE KINDS
    USE KIND_DAT
! 
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE SPCCOM_DAT ! Spectral range data
!
! SUBROUTINES
    USE GRDFIL_SUB ! Open .grd file and check contents
    USE LEXIST_FNC ! Check if file exists
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMDEF ! Default name of LUT file
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)       :: IPT    ! Location of '*' character in NAMDEF
    INTEGER(I4)       :: ISPC   ! Spectral range counter
    CHARACTER(LENREC) :: NAMGRD ! Name of GRD file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL WRTLOG ( 'I-GRDDEF: Checking default .grd files' )
!
  IPT = INDEX ( NAMDEF, '*' )
  IF ( IPT .EQ. 0 ) STOP 'F-GRDDEF: Logical error'
!      
  DO ISPC = 1, NSPC                           ! Loop over spectral ranges
    IF ( SPC(ISPC)%IGD .NE. 0 ) CYCLE         ! .grd file already assigned
    NAMGRD = NAMDEF(1:IPT-1) // TRIM ( SPC(ISPC)%LAB ) // NAMDEF(IPT+1:)
    CALL WRTLOG ( 'I-GRDDEF: looking for file:' // NAMGRD )
    IF ( LEXIST ( NAMGRD ) ) THEN
      CALL GRDFIL ( NAMGRD, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    ELSE
      CALL WRTLOG ( '... file not found' )
    END IF
  END DO
!
END SUBROUTINE GRDDEF
END MODULE GRDDEF_SUB
