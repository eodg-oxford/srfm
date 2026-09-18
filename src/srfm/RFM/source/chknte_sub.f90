MODULE CHKNTE_SUB
CONTAINS
SUBROUTINE CHKNTE ( FAIL, ERRMSG )
!
! VERSION
!   24APR26 AD Also warn if no vib indices set for non-LTE.
!   21APR26 AD Original.
!
! DESCRIPTION
!   Check energy assigned to all required vibrational levels
!   Called by DRVCHK if NTE flag enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE HFLCOM_DAT ! HITRAN file data
    USE NTECOM_DAT ! Non-LTE data
    USE VIBCOM_DAT, ONLY: IDXMOL
!
! SUBROUTINES
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL     :: USENTE ! T=use non-LTE for this molecule
    INTEGER(I4) :: IDM    ! HITRAN ID for molecule
    INTEGER(I4) :: IGAS   ! Index in GASCOM for molecule 
    INTEGER(I4) :: IHFL   ! Index in HFLCOM of HITRAN data file for molecule
    INTEGER(I4) :: INTE   ! Index in NTE kinetic temperature profiles
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO INTE = 1, NNTE 
    IF ( NTE(INTE)%ENG .EQ. 0.0 ) THEN
      FAIL = .TRUE.
      ERRMSG = 'CHKNTE: No Energy assigned to ' // TRIM(NTE(INTE)%COD)
      RETURN
    END IF
  END DO
! 
  DO IGAS = 1, NGAS
    IF ( GAS(IGAS)%NTE ) THEN
      IDM = GAS(IGAS)%IDM
      IHFL = IFLIDM(IDM) 
      IF ( HFL(IHFL)%TYP .EQ. 'PAR' ) THEN
        IF ( ALLOCATED ( IDXMOL ) ) THEN
          USENTE = ANY ( IDXMOL .EQ. IDM ) 
        ELSE
          USENTE = .FALSE.
        END IF
      ELSE
        USENTE = HFL(IHFL)%TYP .EQ. 'BIN'
      END IF
      IF ( .NOT. USENTE ) THEN
        GAS(IGAS)%NTE = .FALSE.
        CALL WRTLOG ( 'W-CHKNTE: ' // TRIM(GAS(IGAS)%COD) // &
                      ': Switching to LTE - no vib.indices loaded' )
      END IF
    END IF
  END DO
!
  FAIL = .FALSE.
!
END SUBROUTINE CHKNTE
END MODULE CHKNTE_SUB
