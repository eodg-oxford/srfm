MODULE CIACHK_SUB
CONTAINS
SUBROUTINE CIACHK ( ID1, ID2, TEM, WNL, WNU, FAIL, ERRMSG )
!
! VERSION
!   21JUN24 AD Change combination with continua to warning rather than fatal    
!   24JUN19 AD Check CO2-CO2 not used with CO2 continuum, similarly H2O,N2,O2
!              Also check Air is unique broadener in spectral range. Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Check CIA data consistent with other RFM options.
!   Called by CIAFIL for each .cia file loaded
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CIACOM_DAT ! Collision-induced absorption data
    USE GASCOM_DAT ! Molecule and isotope data
    USE IDXCON_DAT ! RFM/HITRAN indices for specific molecules
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
    USE NAMGAS_FNC ! Return molecule name + (iso) associated with GASCOM index
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: ID1    ! Index of first molecule
    INTEGER(I4),   INTENT(IN)  :: ID2    ! Index of second molecule
    REAL(R4),      INTENT(IN)  :: TEM    ! Temperature [K] of dataset
    REAL(R8),      INTENT(IN)  :: WNL    ! Lower Wno [cm-1] of dataset
    REAL(R8),      INTENT(IN)  :: WNU    ! Upper Wno [cm-1] of dataset
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: TEMTOL = 0.1 ! Tolerance [K] to separate temperatures
!
! LOCAL VARIABLES
    LOGICAL     :: WRNCO2 = .TRUE.   ! Warn if CO2 CIA/CTM combination
    LOGICAL     :: WRNH2O = .TRUE.   ! Warn if H2O CIA/CTM combination
    LOGICAL     :: WRNN2  = .TRUE.   ! Warn if N2  CIA/CTM combination
    LOGICAL     :: WRNO2  = .TRUE.   ! Warn if O2  CIA/CTM combination
    INTEGER(I4) :: ICIA   ! Indices of CIA data sets
    INTEGER(I4) :: IGAS   ! Index of molecule in /GASCOM/ arrays
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  FAIL = .TRUE.
!
! Warn if CIA being combined with continua representing same feature
  IF ( ID1 .EQ. IDXH2O .AND. ID2 .EQ. IDXH2O .AND. WRNH2O ) THEN
    IGAS = IGSMOL(IDXH2O)
    IF ( GAS(IGAS)%CTM ) CALL WRTLOG ( &
      'W-CIACHK: H2O-H2O CIA may not be compatible with H2O continuum' )
    WRNH2O = .FALSE.
  END IF
  IF ( ID1 .EQ. IDXCO2 .AND. ID2 .EQ. IDXCO2 .AND. WRNCO2 ) THEN
    IGAS = IGSMOL(IDXCO2) 
    IF ( GAS(IGAS)%CTM ) CALL WRTLOG ( &
      'W-CIACHK: CO2-CO2 CIA may not be compatible with CO2 continuum' )
    WRNCO2 = .FALSE.
  END IF
  IF ( ID1 .EQ. IDXN2 .AND. ID2 .EQ. IDXN2 .AND. WRNN2 ) THEN
    IGAS = IGSMOL(IDXN2)
    IF ( GAS(IGAS)%CTM ) CALL WRTLOG ( &
      'W-CIACHK: N2-N2 CIA may not be compatible with N2 continuum' )
    WRNN2 = .FALSE.
  END IF 
  IF ( ID1 .EQ. IDXO2 .AND. ID2 .EQ. IDXO2 .AND. WRNO2 ) THEN
    IGAS = IGSMOL(IDXO2)
    IF ( GAS(IGAS)%CTM ) CALL WRTLOG ( &
      'W-CIACHK: O2-O2 CIA incompatible with O2 continuum' )
    WRNO2 = .FALSE.
  END IF
!
! Check no overlaps with previously loaded data
  DO ICIA = 1, NCIA
    IF ( ID1 .NE. CIA(ICIA)%ID1 .OR. ID2 .NE. CIA(ICIA)%ID2 ) CYCLE
    IF ( WNL .LE. CIA(ICIA)%WNU .AND. WNU .GT. CIA(ICIA)%WNL .AND. &
         ABS ( TEM - CIA(ICIA)%TEM ) .LT. TEMTOL ) THEN
      ERRMSG = 'F-CIACHK: Overlapping CIA datasets for ' // &
               TRIM ( C9REAL ( WNL ) ) // '-' // &
               TRIM ( C9REAL ( WNU ) ) //' cm-1, T=' // C9REAL ( TEM ) 
      RETURN
    END IF
  END DO     
!
! Check air and specific molecule(s) not both used as broadener in same 
! spectral range
  DO ICIA = 1, NCIA-1
    IF ( ID1 .EQ. CIA(ICIA)%ID1 .AND. &   ! look for matching molec#1 and 
         WNL .LE. CIA(ICIA)%WNU .AND. &   ! overlapping spectral range
         WNU .GE. CIA(ICIA)%WNL          ) THEN
      IF ( ID2 .EQ. IDXAIR .AND. CIA(ICIA)%ID2 .NE. IDXAIR ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-CIACHK: Cannot have both air and ' // &
                 TRIM ( NAMGAS ( CIA(ICIA)%IG2 ) )  // &
                 ' as CIA broadeners'
        RETURN
      ELSE IF ( ID2 .NE. IDXAIR .AND. CIA(ICIA)%ID2 .EQ. IDXAIR ) THEN
        FAIL = .TRUE. 
        ERRMSG = 'F-CIACHK: Cannot have both air and ' // &
                 TRIM ( NAMGAS ( IGSMOL(ID2) ) ) // &
                 ' as CIA broadeners'
        RETURN
      END IF
    END IF
  END DO
!
! Normal exit
  FAIL = .FALSE.
!
END SUBROUTINE CIACHK
END MODULE CIACHK_SUB

