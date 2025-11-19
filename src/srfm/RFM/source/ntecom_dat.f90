MODULE NTECOM_DAT
!
! VERSION
!   15SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Non-LTE data
!   Loaded by ADDNTE
!   Non-LTE can be set by retrievals or climatology files
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL CONSTANTS
    INTEGER(I4), PARAMETER :: LENNTE = 14 ! Max length of NTE label
!
  TYPE :: NTETYP
    INTEGER(I4)       :: IDI ! Nonlte isotope id
    INTEGER(I4)       :: IDL ! Nonlte vib.level id
    INTEGER(I4)       :: IDM ! Nonlte molecule id 
    INTEGER(I4)       :: IDX ! Coded Gas/Isotope/IGQ info
    INTEGER(I4)       :: IQF ! Index of Vib.Partition Fn (else 0) 
    REAL(R4)          :: ENG ! Nonlte NTE level energy [cm-1]
    CHARACTER(LENNTE) :: COD ! Label for Vib Temp.
  END TYPE NTETYP
!
! GLOBAL VARIABLES
    TYPE(NTETYP), ALLOCATABLE :: NTE(:)
!
    INTEGER(I4) :: NNTE = 0   ! No. of nonLTE vib.temp profiles
!
END MODULE NTECOM_DAT
