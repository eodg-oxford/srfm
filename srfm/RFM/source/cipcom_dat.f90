MODULE CIPCOM_DAT
!
! VERSION
!   21JUL24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   CIA path data
!   Equivalent homogeneous path segments required for Collision Induced 
!   Absorption spectroscopic calculations.
!   Set in CIAPTH
!   Used in SPCCIA
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: CIPTYP
    INTEGER(I4) :: ICL ! Index of equivalent calc path in CLCCOM
    INTEGER(I4) :: ID1 ! RFM/HITRAN index of molec#1
    INTEGER(I4) :: ID2 ! RFM/HITRAN index of molec#2
    REAL(R4)    :: AM2 ! Absorber amount scale factor
    REAL(R4)    :: TEM ! Path Temperature [K] 
  END TYPE CIPTYP
!
! GLOBAL VARIABLES
    TYPE(CIPTYP), ALLOCATABLE :: CIP(:)
!
    INTEGER(I4) :: NCIP = 0   ! No. of CIP tables stored
!
END MODULE CIPCOM_DAT
