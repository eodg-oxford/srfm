MODULE DRVSVD_SUB
CONTAINS
SUBROUTINE DRVSVD ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   20JUN24 AD Checked.
!   20DEC17 AD Adapted from DRVLUT. Checked.
!
! DESCRIPTION
!   Read RFM driver table *SVD section
!   Called by RFMDRV if *SVD key listed in RFM Driver table and SVDFLG=T
!   Read list of Look-Up tables of absorption coefficients.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE SFLCOM_DAT ! SVD-compressed LUT files
    USE GASCOM_DAT, ONLY: NGAS ! No. different absorbers
    USE SPCCOM_DAT, ONLY: NSPC ! No. spectral ranges
!
! SUBROUTINES
    USE NXTFLD_SUB ! Load next field from section of driver file
    USE SVDDEF_SUB ! Use default SVD filename to find any missing files
    USE SVDFIL_SUB ! Check if SVD file is required
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
    CHARACTER(LENREC) :: NAMDEF = '' ! Default filename template
    CHARACTER(LENREC) :: NAMSVD      ! Filename read from *SVD section
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Allocate matrix indicating which SVD files are used with each molec,range
  ALLOCATE ( IDXSFL(NGAS,NSPC) )
  IDXSFL = 0
!
  DO
    CALL NXTFLD ( LUNDRV, NAMSVD, LENGTH, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
!
    IF ( INDEX ( NAMSVD(1:LENGTH), '*' ) .EQ. 0 ) THEN ! explicit filename
      CALL SVDFIL ( NAMSVD(1:LENGTH), FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    ELSE IF ( NAMDEF .EQ. '' ) THEN                    ! filename template
      NAMDEF = NAMSVD
    ELSE
      ERRMSG = 'F-DRVSVD: *SVD section has more than one filename template'
      FAIL = .TRUE.
    END IF
  END DO
!
  IF ( NAMDEF .NE. '' ) THEN
    CALL SVDDEF ( NAMDEF, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
  END IF
!
 IF ( NSFL .EQ. 0 ) &
    CALL WRTLOG ( 'W-DRVSVD: No useful filenames in *SVD section' )
!
! Normal exit
!
END SUBROUTINE DRVSVD
END MODULE DRVSVD_SUB
