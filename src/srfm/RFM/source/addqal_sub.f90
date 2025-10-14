MODULE ADDQAL_SUB
CONTAINS
SUBROUTINE ADDQAL ( IDXMOL, ISOLST, ILS, IUS )
!
! VERSION
!   26JUN24 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Add to list of molec iso/band qualifiers
!   Called by GASQAL for each qualifier
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE QALCOM_DAT ! Band/isotope selection qualifiers
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IDXMOL    ! HITRAN/RFM index of molecule
    INTEGER(I4), INTENT(IN) :: ISOLST(:) ! List of isotopes
    INTEGER(I4), INTENT(IN) :: ILS       ! Lower vib level index
    INTEGER(I4), INTENT(IN) :: IUS       ! Upper vib level index
!
! LOCAL VARIABLES
    INTEGER(I4) :: IISO ! Isotope counter
    INTEGER(I4) :: IQAL ! Qualifier counter
    INTEGER(I4) :: NISO ! No.of isotopes in ISOLST
    TYPE(QALTYP), ALLOCATABLE :: QALSAV(:) ! Saved QAL during reallocation
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  NISO = SIZE ( ISOLST )
  IQAL = NQAL
!
  IF ( ALLOCATED ( QAL ) ) CALL MOVE_ALLOC ( QAL, QALSAV )
  NQAL = NQAL + NISO
  ALLOCATE ( QAL(NQAL) )
  IF ( ALLOCATED ( QALSAV ) ) QAL(1:IQAL) = QALSAV
!
  DO IISO = 1, NISO
    IQAL = IQAL + 1
    QAL(IQAL)%IDM = IDXMOL
    QAL(IQAL)%ISO = ISOLST(IISO)
    QAL(IQAL)%ILS = ILS
    QAL(IQAL)%IUS = IUS
  END DO
!
END SUBROUTINE ADDQAL
END MODULE ADDQAL_SUB
