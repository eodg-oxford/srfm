MODULE XSCDEF_SUB
CONTAINS
SUBROUTINE XSCDEF ( NAMDEF, FAIL, ERRMSG )
!
! VERSION
!   08JUL24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Use default .xsc filename to find any missing files
!   Called once by DRVXSC if filename template found in *XSC section of 
!   driver table.
!
! VARIABLE KINDS
    USE KIND_DAT
! 
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE LENREC_DAT ! Max length of input text record
    USE SHPCON_DAT, ONLY: SHPXSC ! Line-shape for x/s molecule
!
! SUBROUTINES
    USE LEXIST_FNC ! Check if file exists
    USE XSCFIL_SUB ! Check if .xsc file required
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
    CHARACTER(LENREC) :: NAMXSC ! Name of XSC file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL WRTLOG ( 'I-XSCDEF: Checking default XSC files' )
!
  IPT = INDEX ( NAMDEF, '*' )
  IF ( IPT .EQ. 0 ) STOP 'F-XSCDEF: Logical error'
!      
  DO IGAS = 1, NGAS                           ! Loop over gases
    IF ( GAS(IGAS)%SHP .NE. SHPXSC ) CYCLE    ! Not a x/s molecule
    IF ( GAS(IGAS)%XSC ) CYCLE                ! x/s file already found
    NAMXSC = NAMDEF(1:IPT-1) // TRIM ( GAS(IGAS)%COD ) // NAMDEF(IPT+1:)
    CALL WRTLOG ( 'I-XSCDEF: looking for file: ' // NAMXSC )
    IF ( LEXIST ( NAMXSC ) ) THEN
      CALL XSCFIL ( NAMXSC, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    ELSE
      CALL WRTLOG ( '... file not found' )
    END IF
  END DO
!
END SUBROUTINE XSCDEF
END MODULE XSCDEF_SUB
