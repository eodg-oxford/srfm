MODULE VIBCOM_DAT
!
! VERSION
!   24APR26 AD Original
!
! DESCRIPTION
!   Vibrational level indices
!   List of C15 HITRAN labels for vib levels and associated indices.
!   Use for assigning vib indices to records from HITRAN .par files.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
  SAVE
  PUBLIC :: VIBCOM_RESET
!
! GLOBAL CONSTANTS
!
  TYPE :: VIBTYP
    INTEGER(I4) :: NLV  ! No. vib level labels for this molec. class
    INTEGER(I4),   POINTER :: IDX(:) ! List of Vib level indics
    CHARACTER(15), POINTER :: STR(:) ! List of HITRAN Vib level strings
  END TYPE VIBTYP
!
! GLOBAL VARIABLES
    TYPE(VIBTYP), ALLOCATABLE :: VIB(:)  ! Collection of molec. classes
!
    INTEGER :: NVIB  ! No.different molec. classes
    INTEGER :: NMOL  ! No.different molecules
    INTEGER(I4), ALLOCATABLE :: IDXMOL(:) ! HITRAN indices of molecules
    INTEGER(I4), ALLOCATABLE :: IVBMOL(:) ! Index of corresponding molec.class
!
CONTAINS

  SUBROUTINE VIBCOM_RESET()
    INTEGER(I4) :: I

    IF ( ALLOCATED ( VIB ) ) THEN
      DO I = 1, SIZE ( VIB )
        IF ( ASSOCIATED ( VIB(I)%IDX ) ) DEALLOCATE ( VIB(I)%IDX )
        IF ( ASSOCIATED ( VIB(I)%STR ) ) DEALLOCATE ( VIB(I)%STR )
      END DO
      DEALLOCATE ( VIB )
    END IF
    IF ( ALLOCATED ( IDXMOL ) ) DEALLOCATE ( IDXMOL )
    IF ( ALLOCATED ( IVBMOL ) ) DEALLOCATE ( IVBMOL )

    NVIB = 0
    NMOL = 0
  END SUBROUTINE VIBCOM_RESET

END MODULE VIBCOM_DAT
