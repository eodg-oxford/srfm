 MODULE DRVLUT_SUB
CONTAINS
SUBROUTINE DRVLUT ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   06JUN24 AD Checked.
!   01MAY17 AD F90 conversion of inplut.for. Checked.
!
! DESCRIPTION
!   Read RFM driver table *LUT section
!   Called by RFMDRV if *LUT key listed in RFM Driver table and LUTFLG=T
!   Read list of Look-Up tables of absorption coefficients.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE LFLCOM_DAT ! Look-Up Table files
    USE GASCOM_DAT, ONLY: NGAS ! No. different absorbers
    USE SPCCOM_DAT, ONLY: NSPC ! No. spectral ranges
!
! SUBROUTINES
    USE LUTFIL_SUB ! Check if LUT file is required
    USE LUTDEF_SUB ! Use default LUT filename to find any missing files
    USE NXTFLD_SUB ! Load next field from section of driver file
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
    CHARACTER(LENREC) :: NAMLUT      ! Filename read from *LUT section
    CHARACTER(LENREC) :: NAMDEF = '' ! Default filename template
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Allocate matrix indicating which LUT files are used with each molec,range
  ALLOCATE ( IDXLFL(NGAS,NSPC) )
  IDXLFL = 0
!
  DO
    CALL NXTFLD ( LUNDRV, NAMLUT, LENGTH, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
!
    IF ( INDEX ( NAMLUT(1:LENGTH), '*' ) .EQ. 0 ) THEN ! explicit filename
      CALL LUTFIL ( NAMLUT(1:LENGTH), FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    ELSE IF ( NAMDEF .EQ. '' ) THEN                    ! filename template
      NAMDEF = NAMLUT
    ELSE
      ERRMSG = 'F-DRVLUT: *LUT section has more than one filename template'
      FAIL = .TRUE.
    END IF
  END DO
!
  IF ( NAMDEF .NE. '' ) THEN
    CALL LUTDEF ( NAMDEF, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
  END IF
!
 IF ( NLFL .EQ. 0 ) &
    CALL WRTLOG ( 'W-DRVLUT: No useful filenames in *LUT section' )
!
! Normal exit
!
END SUBROUTINE DRVLUT
END MODULE DRVLUT_SUB
