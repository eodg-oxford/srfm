MODULE CIACOM_DAT
!
! VERSION
!   22JUL24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Collision-induced absorption data
!   Loaded by REACIA
!   Pointers deallocated in RFMDAL
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: CIATYP
    INTEGER(I4) :: ID1 ! RFM/HITRAN index of molec#1
    INTEGER(I4) :: ID2 ! RFM/HITRAN index of molec#2
    INTEGER(I4) :: IG1 ! GASCOM index of molec#1
    INTEGER(I4) :: IG2 ! GASCOM index of molec#2
    INTEGER(I4) :: IGG ! Index in IGGCIA for molec#1,#2 combination
    INTEGER(I4) :: NPT ! No.Abs.Coeff values within each CIA dataset
    REAL(R4)    :: TEM ! Temperature [K] of each dataset
    REAL(R8)    :: WNL ! Lower Wno [cm-1] of each CIA table
    REAL(R8)    :: WNU ! Upper Wno [cm-1] of each CIA table
    REAL(R4), POINTER :: ABS(:) ! [NPT] Absorption data [cm^5/kmole^2]
    REAL(R8), POINTER :: WNO(:) ! [NPT] Wavenumber [cm-1] values for data
  END TYPE CIATYP
!
  TYPE(CIATYP), ALLOCATABLE :: CIA(:)
!
  INTEGER(I4) :: NCIA = 0   ! No. of CIA tables stored
  INTEGER(I4) :: NGGCIA = 0 ! No. different molec-molec combinations
  INTEGER(I4), ALLOCATABLE :: IGGCIA(:,:)  ! IG1,IG2 values of each combination
  INTEGER(I4), ALLOCATABLE :: IDDCIA(:,:)  ! ID1,ID2 values of each combination
!
END MODULE CIACOM_DAT
