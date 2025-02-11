MODULE SHPGAS_SUB 
CONTAINS
SUBROUTINE SHPGAS ( ISHP, GASNAM, FAIL, ERRMSG )
!
! VERSION
!   27NOV23 AD Checked.
!   31MAY18 AD Bug#8: Allow CHI shape to replace MIX shape. Checked.
!   01MAY17 AD F90 conversion of part of inpshp.for. Tested.
!
! DESCRIPTION
!   Set lineshape for molecule
!   Called by DRVSHP for each molecule if SHP flag enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE SHPCON_DAT ! Line-shape indices
    USE IDXCON_DAT, ONLY: IDXCO2 ! HITRAN index of CO2 (=2)
!
! SUBROUTINE  
    USE CHKGAS_SUB ! Check for valid molecule name, isotope, Vib.Level
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: ISHP   ! Index of line shape to be applied
    CHARACTER(*),  INTENT(IN)  :: GASNAM ! Name of molecule
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IDXISO ! Isotope counter
    INTEGER(I4) :: IDXMOL ! HITRAN/RFM Molecule index
    INTEGER(I4) :: IDXVIB ! Global quantum index (should be zero)
    INTEGER(I4) :: IGAS   ! Absorber counter
    INTEGER(I4) :: JGAS   ! Index of isotopes in *GAS arrays
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL CHKGAS ( GASNAM, IDXMOL, IDXISO, IDXVIB, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN
  IGAS = IDXGAS ( IDXMOL ) 
  FAIL = .TRUE.           
  IF ( IDXISO .NE. 0 ) THEN
    ERRMSG = 'F-SHPGAS: Shapes cannot be applied just to isotopes: ' // GASNAM
  ELSE IF ( IDXVIB .NE. 0 ) THEN
    ERRMSG = 'F-SHPGAS: *SHP section appears to contain vib.temp. ' // &
             'instead of molecule name: ' // GASNAM
  ELSE IF ( IGAS .EQ. 0 ) THEN
    ERRMSG = 'F-SHPGAS: Line Shape specified for unlisted gas= ' // GASNAM
  ELSE IF ( GAS(IGAS)%SHP .EQ. SHPCTM ) THEN
    ERRMSG = 'F-SHPGAS: line shape specified for ' // GASNAM // &
             ' but already set to continuum-only'
  ELSE IF ( GAS(IGAS)%SHP .EQ. SHPXSC ) THEN
    ERRMSG = 'F-SHPGAS: line shape specified for ' // GASNAM // &
             ' but already set as cross-section'
  ELSE IF ( GAS(IGAS)%SHP .NE. 0 .AND. .NOT. & ! allow CHI shape to override MIX
            ( ISHP .EQ. SHPCHI .AND. GAS(IGAS)%SHP .EQ. SHPMIX ) ) THEN
    ERRMSG = 'F-SHPGAS: Two different line shapes specified for gas= ' // GASNAM
  ELSE IF ( ISHP .EQ. SHPCHI .AND. IDXMOL .NE. IDXCO2 ) THEN
    ERRMSG = 'F-SHPGAS: Chi-factor only appicable to CO2'
  ELSE
    FAIL = .FALSE.
  END IF
  IF ( FAIL ) RETURN
!
  GAS(IGAS)%SHP = ISHP
!
  DO IDXISO = 1, GAS(IGAS)%NIS        ! Also set any other isotopes
    JGAS = GAS(IGAS)%ISO(IDXISO)
    GAS(JGAS)%SHP = ISHP
  END DO
!
END SUBROUTINE SHPGAS
END MODULE SHPGAS_SUB
