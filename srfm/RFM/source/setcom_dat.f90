MODULE SETCOM_DAT
!
! VERSION
!   25AUG24 AD Checked.
!   11AUG23 AD Original.
!
! DESCRIPTION
!   HITRAN line parameter sets
!   Flags indicating which extra parameters are present in the HITRAN 
!   line data for different line parameter files.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL CONSTANTS
!   List of IPAR indices used for HITRAN line parameters 
!   Parameters 1-NRQPAR are considered mandatory for a line to be usable
!   Others can be set as self to air values, or to zero
!
! Mandatory parameters are those contained in the standard HITRAN .par format
    INTEGER(I4), PARAMETER :: IPAR_WNO = 1  ! 'nu' wavenumber
    INTEGER(I4), PARAMETER :: IPAR_IDM = 2  ! 'molec_id' molecule#
    INTEGER(I4), PARAMETER :: IPAR_IDI = 3  ! 'local_iso_id' isotope#
    INTEGER(I4), PARAMETER :: IPAR_STR = 4  ! 'sw' intenstity
    INTEGER(I4), PARAMETER :: IPAR_ELS = 5  ! 'elower' lower state energy
    INTEGER(I4), PARAMETER :: IPAR_HWA = 6  ! 'gamma_' air-broad half-width
    INTEGER(I4), PARAMETER :: IPAR_TCA = 7  ! 'n_' T-coeff of HWA
    INTEGER(I4), PARAMETER :: IPAR_HWS = 8  ! 'gamma_self' self.br.h-w
    INTEGER(I4), PARAMETER :: IPAR_PSA = 9  ! 'delta_' air pressure shift
    INTEGER(I4), PARAMETER :: NRQPAR   = 9  ! Number of mandatory parameters
!
    INTEGER(I4), PARAMETER :: IPAR_LMA = 10 ! 'y_' air line mixing
    INTEGER(I4), PARAMETER :: IPAR_LMS = 11 ! 'y_self' self line mixing
    INTEGER(I4), PARAMETER :: IPAR_PSS = 12 ! 'delta_self' self press.shift
    INTEGER(I4), PARAMETER :: IPAR_TCS = 13 ! 'n_self' T-coeff of HWS
    INTEGER(I4), PARAMETER :: MAXPAR   = 13 ! max number of parameters
!
! GLOBAL VARIABLES
    INTEGER              :: NSET     ! No. line parameter sets
    LOGICAL, TARGET, ALLOCATABLE :: SET(:,:) ! Line parameter sets
!
END MODULE SETCOM_DAT
