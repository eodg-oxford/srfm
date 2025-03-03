MODULE DRVSHP_SUB
CONTAINS
SUBROUTINE DRVSHP ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   23DEC23 AD Checked.
!   01JUN17 AD F90 conversion of inpshp.for. Checked.
!
! DESCRIPTION
!   Read RFM driver table *SHP section
!   Called by RFMDRV if SHPFLG set TRUE.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE IDXCON_DAT ! RFM/HITRAN indices for specific molecules
    USE LENREC_DAT ! Max length of input text record
    USE SHPCON_DAT ! Line-shape indices
    USE FLGCOM_DAT, ONLY: VVWFLG ! apply Van Vleck-Weisskopf correction
!
! SUBROUTINES
    USE NXTREC_SUB ! Load next record from input file
    USE SHPGAS_SUB ! Set lineshape for molecule
    USE TXTFLD_SUB ! Identify start and end points of text field
    USE UPCASE_FNC ! Convert text string to upper case
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
    LOGICAL           :: ENDSEC ! T=End of driver table section reached.
    INTEGER(I4)       :: IEND   ! Pointer to end of filename in record
    INTEGER(I4)       :: IFLD   ! Text Field# within record
    INTEGER(I4)       :: IGAS   ! Absorber counter
    INTEGER(I4)       :: ISHP   ! Line shape index#
    INTEGER(I4)       :: ISTA   ! Pointer to start of filename in record
    INTEGER(I4)       :: SHPDEF = 0 ! Default Line shape to be used
    CHARACTER(LENGAS) :: FIELD  ! Absorber name (C*7 is longest name in use)
    CHARACTER(LENREC) :: RECORD ! Line of Driver file or External TAN file
    CHARACTER(3)      :: SHAPE  ! Line shape name
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Set all current Voigt lineshape gases to 0 to allow overwriting
  WHERE ( GAS%SHP .EQ. SHPVOI ) GAS%SHP = 0 
!
  DO                             ! Loop over records in driver table section
    CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    IF ( ENDSEC ) EXIT
    CALL TXTFLD ( RECORD, 1, ISTA, IEND )    ! Field#1 = Shape 
    IEND = MIN ( ISTA+2, IEND )                  ! Limit to C*3
    SHAPE = UPCASE ( RECORD(ISTA:IEND) ) 
    SELECT CASE ( SHAPE ) 
    CASE ( 'VOI' ) ; ISHP = SHPVOI
    CASE ( 'LOR' ) ; ISHP = SHPLOR
    CASE ( 'DOP' ) ; ISHP = SHPDOP
    CASE ( 'CHI' ) ; ISHP = SHPCHI
    CASE ( 'VAN' ) ; ISHP = SHPVVW
    CASE ( 'VVW' ) ; ISHP = SHPVVW
    CASE DEFAULT
      FAIL = .TRUE.                                  
      ERRMSG = 'F-DRVSHP: Unrecognised Line Shape: ' // SHAPE
      RETURN
    END SELECT
!
    IF ( ISHP .EQ. SHPVVW .AND. VVWFLG ) &
      CALL WRTLOG ( 'W-DRVSHP: VVW correction flg will be ignored ' // &
                    'for lines calculated with VVW shape' )
! Subsequent fields in same record should contain list of gases with spec.shape
    IFLD = 2
    DO
      CALL TXTFLD ( RECORD, IFLD, ISTA, IEND )
      IF ( ISTA .EQ. 0 ) EXIT
      FIELD = RECORD(ISTA:IEND)
      IF ( FIELD .EQ. '*' ) THEN
        IF ( SHPDEF .EQ. 0 ) THEN
          SHPDEF = ISHP
          CALL WRTLOG ( 'I-DRVSHP: Setting default ' // SHAPE // &
                        ' lineshape for any remaining gases' )
        ELSE
          FAIL = .TRUE.
          ERRMSG = 'F-DRVSHP: Two different default line shapes specified'
          RETURN
        END IF
      ELSE        ! assume molecule name
        CALL SHPGAS ( ISHP, FIELD, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        CALL WRTLOG ( 'I-DRVSHP: Using ' // SHAPE // ' lineshape for ' // &
                      'gas=' // FIELD )
      END IF
      IFLD = IFLD + 1
    END DO            ! Get next field in record
  END DO              ! Get next record in section
!
! Reached end of section - assign default shape to any gases not yet assigned
  IF ( SHPDEF .EQ. 0 ) SHPDEF = SHPVOI
  WHERE ( GAS%SHP .EQ. 0 ) GAS%SHP = SHPDEF
!
! Issue warning if inconsistent shape + continuum used for CO2 or H2O
  IF ( IGSMOL(IDXCO2) .NE. 0 ) THEN
    IGAS = IGSMOL(IDXCO2) 
    IF ( GAS(IGAS)%CTM .AND. &
         GAS(IGAS)%SHP .NE. SHPCTM .AND. &
         GAS(IGAS)%SHP .NE. SHPCHI          ) &
      CALL WRTLOG ( 'W-DRVSHP: CHI-FACTOR recommended ' // &
                    'for use with CO2 continuum' )
  ENDIF 
!
  IF ( IGSMOL(IDXH2O) .NE. 0 ) THEN
    IGAS = IGSMOL(IDXH2O) 
    IF ( GAS(IGAS)%CTM .AND. &
         GAS(IGAS)%SHP .NE. SHPCTM .AND. &
         GAS(IGAS)%SHP .NE. SHPLOR .AND. &
         GAS(IGAS)%SHP .NE. SHPVOI          ) &
      CALL WRTLOG ( 'W-DRVSHP: LORENTZ or VOIGT line-shape recommended ' // &
                    'for use with H2O continuum' ) 
  END IF
!
  FAIL = .FALSE.
!
END SUBROUTINE DRVSHP
END MODULE DRVSHP_SUB
