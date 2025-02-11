MODULE ILSCOM_DAT
!
! VERSION
!   28MAY24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Instrument Lineshape data
!   Loaded by ILSFIL
!   Pointer FNC deallocated in RFMDAL.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: ILSTYP
    INTEGER(I4) :: NPT ! No. tabulation pts for each ILS Fn.
    REAL(R8)    :: PT1 ! Lower point for each ILS function
    REAL(R8)    :: PT2 ! Upper point for each ILS function
    REAL(R8)    :: PTD ! Point spacing for each ILS function
    REAL(R8)    :: WNL ! Lower wavenumber for ILS validity
    REAL(R8)    :: WNU ! Upper wavenumber for ILS validity
    REAL(R8), POINTER :: FNC(:) ! [NPT] Tabulated ILS function
  END TYPE ILSTYP
!
! GLOBAL VARIABLES
    TYPE(ILSTYP), ALLOCATABLE :: ILS(:)
!
    INTEGER(I4) :: NILS = 0   ! No. ILS functions stored
    INTEGER(I4) :: IDFILS = 0 ! Index of default ILS fn (0=none)
!
END MODULE ILSCOM_DAT

