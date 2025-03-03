MODULE QFNCOM_DAT
!
! VERSION
!   29JAN24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Non-LTE Vib.Partition Fn data
!   Loaded by ADDQFN
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: QFNTYP
    INTEGER(I4) :: IDM ! Nonlte molecule id 
    INTEGER(I4) :: IDI ! Nonlte isotope id
  END TYPE QFNTYP

! GLOBAL VARIABLES
    TYPE(QFNTYP), ALLOCATABLE :: QFN(:)
!
    INTEGER(I4) :: NQFN = 0   ! No. of nonLTE Vib Partition Fn profiles
!
END MODULE QFNCOM_DAT
