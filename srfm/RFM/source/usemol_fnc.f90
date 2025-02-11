MODULE USEMOL_FNC
CONTAINS
LOGICAL PURE FUNCTION USEMOL  ( IDXMOL )
!
! VERSION
!   06SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   T = molecule is required
!   General purpose function.
!   IDXMOL can be any integer.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT, ONLY: IGSMOL, MAXMOL ! Molecule and isotope data
!
  IMPLICIT NONE 
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IDXMOL ! HITRAN/RFM index of molecule
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  IF ( IDXMOL .LE. 0 .OR. IDXMOL .GT. MAXMOL ) THEN
    USEMOL = .FALSE.
  ELSE
    USEMOL = IGSMOL(IDXMOL) .GT. 0
  END IF
!
END FUNCTION USEMOL
END MODULE USEMOL_FNC

