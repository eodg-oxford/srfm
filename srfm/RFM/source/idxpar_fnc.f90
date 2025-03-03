MODULE IDXPAR_FNC
CONTAINS
INTEGER(I4) FUNCTION IDXPAR ( HITFLD )
!
! VERSION
!   16AUG24 AD Checked.
!   11AUG23 AD Use SETCOM instead of HDBCOM for IPAR_* values
!   30MAY23 AD Original.
!
! DESCRIPTION
!   Return IPAR value associated with HITRAN database line parameter
!   Called by HDBCHK
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE SETCOM_DAT ! HITRAN line parameter sets
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*), INTENT(IN) :: HITFLD ! Field from HITRAN file header record
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  SELECT CASE ( HITFLD ) 
  CASE ( 'nu' )           ; IDXPAR = IPAR_WNO
  CASE ( 'molec_id' )     ; IDXPAR = IPAR_IDM
  CASE ( 'local_iso_id' ) ; IDXPAR = IPAR_IDI
  CASE ( 'sw' )           ; IDXPAR = IPAR_STR
  CASE ( 'elower' )       ; IDXPAR = IPAR_ELS
  CASE DEFAULT
    IF ( INDEX ( HITFLD, '-err' ) .GT. 0 .OR. &
         INDEX ( HITFLD, '-ref' ) .GT. 0  ) THEN
      IDXPAR = 0
    ELSE IF ( HITFLD(1:6) .EQ. 'gamma_' ) THEN
      IF ( HITFLD(1:10) .EQ. 'gamma_self' ) THEN
        IDXPAR = IPAR_HWS
      ELSE
        IDXPAR = IPAR_HWA
      END IF
    ELSE IF ( HITFLD(1:2) .EQ. 'n_' ) THEN
      IF ( HITFLD(1:6) .EQ. 'n_self' ) THEN
        IDXPAR = IPAR_TCS
      ELSE
        IDXPAR = IPAR_TCA
      END IF
    ELSE IF ( HITFLD(1:6) .EQ. 'delta_' ) THEN
      IF ( HITFLD(1:10) .EQ. 'delta_self' ) THEN
        IDXPAR = IPAR_PSS
      ELSE
        IDXPAR = IPAR_PSA
      END IF
    ELSE IF ( HITFLD(1:2) .EQ. 'y_' ) THEN
      IF ( HITFLD(1:6) .EQ. 'y_self' ) THEN
        IDXPAR = IPAR_LMS
      ELSE
        IDXPAR = IPAR_LMA
      END IF
    ELSE
      IDXPAR = 0
    END IF
  END SELECT
!
END FUNCTION IDXPAR
END MODULE IDXPAR_FNC
