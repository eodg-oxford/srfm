MODULE NAMMOL_FNC
CONTAINS
CHARACTER(LENGAS+2) PURE FUNCTION NAMMOL ( IDXMOL, IDXISO )
!
! VERSION
!   04SEP24 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Return molecule name + i[#iso]
!   General purpose module.
!   Assumes IDXMOL, IDXISO is currently loaded in GASCOM.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT, ONLY: LENGAS ! Max length of molecule name
!
! SUBROUTINES
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
    USE NAMGAS_FNC ! Return molecule name + (iso) associated with GASCOM index
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),           INTENT(IN) :: IDXMOL ! HITRAN/RFM index of molecule
    INTEGER(I4), OPTIONAL, INTENT(IN) :: IDXISO ! HITRAN isotope#
!
! LOCAL VARIABLES
    INTEGER :: IGAS ! Index of GAS in GASCOM 
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  IF ( PRESENT ( IDXISO ) ) THEN
    IGAS = IDXGAS ( IDXMOL, IDXISO ) 
  ELSE
    IGAS = IDXGAS ( IDXMOL ) 
  END IF
! can't have STOP statement in PURE function
!  IF ( IGAS .EQ. 0 ) STOP 'F-NAMMOL: Logical error'
!
  NAMMOL = NAMGAS ( IGAS ) 
!
END FUNCTION NAMMOL
END MODULE NAMMOL_FNC

