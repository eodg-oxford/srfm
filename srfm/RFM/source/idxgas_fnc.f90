MODULE IDXGAS_FNC
CONTAINS
INTEGER(I4) PURE FUNCTION IDXGAS ( IDXMOL, IDXISO )
!
! VERSION
!   22MAR23 AD Checked.
!   01MAY17 AD F90 conversion of igasid.for. Checked.
!
! DESCRIPTION
!   Index in GAS arrays of molecule,isotope
!   General purpose module.
!   Returns 0 if not assigned
!   If IDXISO not present
!           returns IGSMOL(IDXMOL) if IDXMOL used, else 0
!   If IDXISO is present and = -1
!           returns IGSMOL(IDXMOL only if gas not split into isotopes, else 0
!   If IDXISO is present and >= 0
!           returns GAS(IGAS)%ISO(IDXISO), ie index of isotope if split
!           or just IGSMOL(IDXMOL) if not split or IDXISO=0
!   (IDXISO=0 is equivalent to IDXISO not present)
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL VARIABLES
    USE GASCOM_DAT ! Molecule and isotope data
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),           INTENT(IN) :: IDXMOL ! HITRAN/RFM index of molecule
    INTEGER(I4), OPTIONAL, INTENT(IN) :: IDXISO ! HITRAN isotope#
!
! LOCAL VARIABLES
    INTEGER(I4) :: IGAS ! Temporary index
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  IGAS = IGSMOL(IDXMOL)
!
  IF ( PRESENT ( IDXISO ) ) THEN
    IF ( IDXISO .NE. 0 ) THEN
      IF ( IDXISO .EQ. -1 ) THEN  ! only return IGAS if not split 
        IF ( ISOMOL(IDXMOL) ) IGAS = 0
      ELSE
! can't have STOP statement in PURE function
!        IF ( IDXMOL .GT. MAXHLN ) STOP 'F-IDXGAS: logical error'
        IGAS = GAS(IGAS)%ISO(IDXISO)
      END IF
    END IF
  END IF
!
  IDXGAS = IGAS
!
END FUNCTION IDXGAS
END MODULE IDXGAS_FNC


