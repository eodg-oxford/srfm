MODULE PTHCOM_DAT
!
! VERSION
!   24JAN24 AD Checked.   
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Path segment data
!   Initialised in RFMPTH
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: PTHTYP
    LOGICAL     :: NTE ! T=non-LTE molecule, F=LTE
    INTEGER(I4) :: IAT ! Atmospheric profile layer# for this path
    INTEGER(I4) :: ICL ! Index of corresponding CLC path
    INTEGER(I4) :: IDR ! 0=downward from obs, or irrelevant, 1=upward
    INTEGER(I4) :: IGS ! Absorbing gas# for this path
    INTEGER(I4) :: ITN ! Tangent Height# for this path
    REAL(R4)    :: AMT ! Path gas amount [ kmol / cm^2 ] 
    REAL(R4)    :: PPA ! Partial pressure [atm]
    REAL(R4)    :: PRE ! Pressure [atm]
    REAL(R4)    :: PSI ! Horiz.angle [deg] of low alt end of path (GRA flag)
    REAL(R4)    :: PSU ! Abs.weighted Horiz.angle [deg] of path
    REAL(R4)    :: RAY ! Ray length in path [km]
    REAL(R4)    :: TEM ! Temperature [K]	
  END TYPE PTHTYP
!
! GLOBAL VARIABLES
    TYPE(PTHTYP), ALLOCATABLE :: PTH(:)
!
    LOGICAL     :: USEDIR = .FALSE. ! T=use IDR to distinguish paths
    INTEGER(I4) :: NPTH             ! No. of paths used
!
END MODULE PTHCOM_DAT
