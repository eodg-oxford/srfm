MODULE CLCCOM_DAT
!
! VERSION
!   26JAN24 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Calculated path segments
!   Set of Curtis-Godson equivalent paths for which absorption coefficients
!   are calculated. These all have unit absorber amount, and the optical depths
!   of paths used for radiative transfer are scaled by absorber amount.
!   
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: CLCTYP
    LOGICAL     :: LBL ! T=line-by-line calculation, F=other
    INTEGER(I4) :: IGS ! Index in GAS of absorber
    REAL(R4)    :: AMT ! Absorber amount [kmol/cm^2]
    REAL(R4)    :: PPA ! Partial pressure [atm]
    REAL(R4)    :: PRE ! Pressure [atm] 
    REAL(R4)    :: TEM ! Temperature [K]
  END TYPE CLCTYP
!
! GLOBAL VARIABLES
    TYPE(CLCTYP), ALLOCATABLE :: CLC(:)
!
    INTEGER(I4) :: NCLC = 0 ! No. calculated path segments
!
END MODULE CLCCOM_DAT
