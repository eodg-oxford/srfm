MODULE NTEDEF_SUB
CONTAINS
SUBROUTINE NTEDEF ( NAMDEF, FAIL, ERRMSG )
!
! VERSION
!   07JUL24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Use default .nte filename to find any missing files
!   Called once by DRVNTE if filename template found in *NTE section of 
!   driver table.
!
! VARIABLE KINDS
    USE KIND_DAT
! 
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE LENREC_DAT ! Max length of input text record
!
! SUBROUTINES
    USE LEXIST_FNC ! Check if file exists
    USE NTEFIL_SUB ! Read NTE Vibrational Temperature data from file
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
    INTEGER(I4)       :: IDXMOL ! HITRAN/RFM Molecule ID
    INTEGER(I4)       :: IGAS   ! Absorber counter
    INTEGER(I4)       :: IPT    ! Location of '*' character in NAMDEF
    CHARACTER(LENREC) :: NAMNTE ! Name of NTE file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL WRTLOG ( 'I-NTEDEF: Checking default NTE files' )
!
  IPT = INDEX ( NAMDEF, '*' )
  IF ( IPT .EQ. 0 ) STOP 'F-NTEDEF: Logical error'
!      
  DO IGAS = 1, NGAS                   ! Loop over gases
    IDXMOL = GAS(IGAS)%IDM
    IF ( IDXMOL .GT. MAXHLN ) CYCLE   ! ignore any cross-section molecule
    IF ( GAS(IGAS)%NTE ) CYCLE        ! Vib Temp data already loaded
    NAMNTE = NAMDEF(1:IPT-1) // TRIM ( GAS(IGAS)%COD ) // NAMDEF(IPT+1:)
    CALL WRTLOG ( 'I-NTEDEF: looking for file:' // NAMNTE )
    IF ( LEXIST ( NAMNTE ) ) THEN
      CALL NTEFIL ( NAMNTE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    ELSE
      CALL WRTLOG ( '... file not found' )
    END IF
  END DO
!
END SUBROUTINE NTEDEF
END MODULE NTEDEF_SUB
