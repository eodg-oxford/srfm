MODULE TXTPAR_FNC
CONTAINS
FUNCTION TXTPAR ( IPAR )
!
! VERSION
!   27AUG24 AD Checked.
!   11AUG23 AD Use SETCOM instead of HDBCOM for IPAR values
!   31MAY23 AD Original.
!
! DESCRIPTION
!   Return text string describing HITRAN database line parameter
!   Called by HDBCHK.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE SETCOM_DAT ! HITRAN line parameter sets
!
  IMPLICIT NONE
!
! FUNCTION TYPE
   CHARACTER(:), ALLOCATABLE :: TXTPAR
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IPAR ! Index of HITRAN parameter
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  SELECT CASE ( IPAR ) 
    CASE ( IPAR_WNO ) ; TXTPAR = 'nu (wavenumber)'
    CASE ( IPAR_IDM ) ; TXTPAR = 'molec_id'
    CASE ( IPAR_IDI ) ; TXTPAR = 'local_iso_id'
    CASE ( IPAR_STR ) ; TXTPAR = 'sw (intensity)'
    CASE ( IPAR_ELS ) ; TXTPAR = 'elower (energy)'
    CASE ( IPAR_HWA ) ; TXTPAR = 'gamma_ (halfwidth)'
    CASE ( IPAR_TCA ) ; TXTPAR = 'n_ (T-coeff of h/w)'
    CASE ( IPAR_HWS ) ; TXTPAR = 'gamma_self (self h/w)'
    CASE ( IPAR_TCS ) ; TXTPAR = 'n_self (T-coeff of self h/w)'
    CASE ( IPAR_PSA ) ; TXTPAR = 'delta_ (press.shift)'
    CASE ( IPAR_PSS ) ; TXTPAR = 'delta_self (self press.shift)'
    CASE ( IPAR_LMA ) ; TXTPAR = 'y_ (line coupling)'
    CASE ( IPAR_LMS ) ; TXTPAR = 'y_self (self line coupling)'
    CASE DEFAULT ; STOP 'F-TXTPAR: Logical error'
  END SELECT
!
END FUNCTION TXTPAR
END MODULE TXTPAR_FNC
