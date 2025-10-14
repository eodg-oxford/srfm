MODULE LUTDEF_SUB
CONTAINS
SUBROUTINE LUTDEF ( NAMDEF, FAIL, ERRMSG )
!
! VERSION
!   21JUN24 AD Checked.
!   20DEC17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Use default LUT filename to find any missing files
!   Called once by DRVLUT if filename template found in *LUT section of 
!   driver table.
!
! VARIABLE KINDS
    USE KIND_DAT
! 
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE SPCCOM_DAT ! Spectral range data
    USE GASCOM_DAT, ONLY: NGAS ! No. different molecules used
    USE LFLCOM_DAT, ONLY: IDXLFL ! Indices of assigned LFL files
!
! SUBROUTINES
    USE LEXIST_FNC ! Check if file exists
    USE LUTFIL_SUB ! Check if LUT file is required
    USE NAMGAS_FNC ! Return molecule name + (iso) associated with GASCOM index
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
    CHARACTER(LENREC) :: NAMLUT ! Name of LUT file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL WRTLOG ( 'I-LUTDEF: Checking default LUT files' )
!
  IPT = INDEX ( NAMDEF, '*' )
  IF ( IPT .EQ. 0 ) STOP 'F-LUTDEF: Logical error'
!      
  DO ISPC = 1, NSPC
    DO IGAS = 1, NGAS                           ! Loop over gases
      IF ( IDXLFL(IGAS,ISPC) .GT. 0 ) CYCLE     ! LUT already found
!
! At this point, spectral range/gas found requiring LUT data but no data loaded
! Construct spec/gas part of LUT filename
      NAMLUT = NAMDEF(1:IPT-1) // TRIM ( SPC(ISPC)%LAB ) // &
                 TRIM ( NAMGAS(IGAS) ) // NAMDEF(IPT+1:)
!
      CALL WRTLOG ( 'I-LUTDEF: looking for file:' // NAMLUT )
      IF ( LEXIST ( NAMLUT ) ) THEN
        CALL LUTFIL ( NAMLUT, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      ELSE
        CALL WRTLOG ( '... file not found' )
      END IF
    END DO
  END DO
!
END SUBROUTINE LUTDEF
END MODULE LUTDEF_SUB
