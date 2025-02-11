MODULE SVDDEF_SUB
CONTAINS
SUBROUTINE SVDDEF ( NAMDEF, FAIL, ERRMSG )
!
! VERSION
!   12JUN24 AD Checked.
!   20DEC17 AD Adapted from LUTDEF. Checked.
!
! DESCRIPTION
!   Use default SVD filename to find any missing files
!   Called once by DRVSVD if filename template found in *SVD section of 
!   driver table.
!
! VARIABLE KINDS
    USE KIND_DAT
! 
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE LENREC_DAT ! Max length of input text record
    USE SPCCOM_DAT ! Spectral range data
    USE SFLCOM_DAT, ONLY: IDXSFL ! Indices of assigned SFL files
!
! SUBROUTINES
    USE IDGOLD_FNC ! Convert new RFM index for .xsc data to old value
    USE LEXIST_FNC ! Check if file exists
    USE SVDFIL_SUB ! Check if SVD file is required
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
    INTEGER(I4)       :: IGAS   ! Absorber counter
    INTEGER(I4)       :: IPT    ! Location of '*' character in NAMDEF
    INTEGER(I4)       :: ISPC   ! Counter for spectral ranges
    CHARACTER(2)      :: IDSTR  ! String containing molecule index
    CHARACTER(LENREC) :: NAMSVD ! Name of LUT file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL WRTLOG ( 'I-SVDDEF: Checking default SVD files' )
!
  IPT = INDEX ( NAMDEF, '*' )
  IF ( IPT .EQ. 0 ) STOP 'F-SVDDEF: Logical error'
!      
  DO ISPC = 1, NSPC
    DO IGAS = 1, NGAS                           ! Loop over gases
      IF ( IDXSFL(IGAS,ISPC) .GT. 0 ) CYCLE     ! SVD already found
!
! At this point, spectral range/gas found requiring LUT data but no data loaded
! Construct spec/gas part of SVD filename
!
      WRITE ( IDSTR, '(I2.2)' ) IDGOLD ( GAS(IGAS)%IDM ) 
! for historical reasons SF6 has ID=30 not 64 for MIPAS LUTs
      IF ( IDSTR .EQ. '64' ) IDSTR = '30'
      NAMSVD = NAMDEF(1:IPT-1) // TRIM ( SPC(ISPC)%LAB ) // '_' // & 
               IDSTR // NAMDEF(IPT+1:)
!
      CALL WRTLOG ( 'I-SVDDEF: looking for file:' // NAMSVD )
      IF ( LEXIST ( NAMSVD ) ) THEN
        CALL SVDFIL ( NAMSVD, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      ELSE
        CALL WRTLOG ( '... file not found' )
      END IF
    END DO
  END DO
!
END SUBROUTINE SVDDEF
END MODULE SVDDEF_SUB
