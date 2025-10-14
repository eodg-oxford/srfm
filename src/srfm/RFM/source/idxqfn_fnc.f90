MODULE IDXQFN_FNC
CONTAINS
INTEGER(I4) FUNCTION IDXQFN ( IDXMOL, IDXISO )
!
! VERSION
!   13SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Index in QFNCOM of molec,isotope
!   General purpose module.
!   Returns 0 if no vibrational partition function assigned.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE QFNCOM_DAT ! Non-LTE Vib.Partition Fn data
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IDXMOL ! HITRAN/RFM index of molecule
    INTEGER(I4), INTENT(IN) :: IDXISO ! HITRAN isotope#
!
! LOCAL VARIABLES
    INTEGER(I4) :: IQFN ! Counter
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  DO IQFN = 1, NQFN
    IF ( QFN(IQFN)%IDM .EQ. IDXMOL .AND. QFN(IQFN)%IDI .EQ. IDXISO ) THEN
      IDXQFN = IQFN
      RETURN
    END IF
  END DO
!
  IDXQFN = 0
!
END FUNCTION IDXQFN
END MODULE IDXQFN_FNC


