MODULE QALCOM_DAT
!
! VERSION
!   16APR23 AD Checked.
!   17NOV17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Band/isotope selection qualifiers
!   Set in ADDQAL
!   Accessed by USEQAL
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: QALTYP
    INTEGER(I4) :: IDM ! HITRAN/RFM molecule index
    INTEGER(I4) :: ISO ! Isotope# (or 0 for all isotopes)
    INTEGER(I4) :: ILS ! Lower vib level index (or 0 for all values)
    INTEGER(I4) :: IUS ! Upper vib level index (or 0 for all values)
  END TYPE QALTYP
!
! GLOBAL VARIABLES
    TYPE(QALTYP), ALLOCATABLE :: QAL(:)
!
    INTEGER(I4) :: NQAL = 0 ! No.of qualifiers set
!
END MODULE QALCOM_DAT
