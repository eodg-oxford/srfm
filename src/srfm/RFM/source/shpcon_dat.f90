MODULE SHPCON_DAT
!
! VERSION
!   08FEB24 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Line-shape indices
!   Integer codes used for different lineshapes.
!   Number are arbitrary - just have to be different.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
!
! GLOBAL CONSTANTS
    INTEGER(I4), PARAMETER :: SHPVOI = 1  ! Voigt line shape 
    INTEGER(I4), PARAMETER :: SHPLOR = 2  ! Lorentz line shape 
    INTEGER(I4), PARAMETER :: SHPDOP = 3  ! Doppler line shape 
    INTEGER(I4), PARAMETER :: SHPCHI = 4  ! Chi-factor scaled Voigt lineshape 
    INTEGER(I4), PARAMETER :: SHPVVW = 5  ! Van Vleck-Weisskopf lineshape 
    INTEGER(I4), PARAMETER :: SHPMIX = 6  ! Voigt+line mixing
    INTEGER(I4), PARAMETER :: SHPCTM = 8  ! Continuum data only
    INTEGER(I4), PARAMETER :: SHPCIA = 9  ! CIA data only
    INTEGER(I4), PARAMETER :: SHPXSC = 10 ! Cross-section data 
    INTEGER(I4), PARAMETER :: SHPAIR = 12 ! Air (for Rayleigh extinction)
!
END MODULE SHPCON_DAT
