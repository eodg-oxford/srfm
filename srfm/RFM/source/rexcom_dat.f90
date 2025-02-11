MODULE REXCOM_DAT
!
! VERSION
!   23SEP24 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Rayleigh extinction path data
!   Filled by REXPTH
!   Used by SPCREX
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: REXTYP
    INTEGER(I4) :: ICL ! Index of corresponding calculated path
    REAL(R4)    :: AMT ! Absorber amount [kmol/cm^2]
    REAL(R4)    :: A0  ! Coefficient for King depolarisation factor
    REAL(R4)    :: A2  ! Coefficient for King depolarisation factor
    REAL(R4)    :: A4  ! Coefficient for King depolarisation factor
  END TYPE REXTYP
!
! GLOBAL VARIABLES
    TYPE(REXTYP), ALLOCATABLE :: REX(:)
!
    INTEGER(I4) :: NREX = 0 ! No. of REX calculated paths
!
END MODULE REXCOM_DAT
