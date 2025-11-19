MODULE ILSCOM_DAT
!
! VERSION
!   28MAY24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Instrument Lineshape data.
!   Loaded by ILSFIL. Pointers are deallocated by RFMDAL and explicitly by
!   ILSCOM_RESET to support consecutive runs without reloading the module.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
  PUBLIC :: ILSCOM_RESET
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
CONTAINS

  SUBROUTINE ILSCOM_RESET()
    INTEGER(I4) :: I

    IF ( ALLOCATED ( ILS ) ) THEN
      DO I = 1, SIZE ( ILS )
        IF ( ASSOCIATED ( ILS(I)%FNC ) ) DEALLOCATE ( ILS(I)%FNC )
      END DO
      DEALLOCATE ( ILS )
    END IF

    NILS   = 0
    IDFILS = 0
  END SUBROUTINE ILSCOM_RESET

END MODULE ILSCOM_DAT
