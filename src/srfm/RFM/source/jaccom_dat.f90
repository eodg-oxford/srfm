MODULE JACCOM_DAT
!
! VERSION
!   28JAN24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Jacobian data
!   Loaded by JACALT.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL CONSTANTS
    INTEGER(I4), PARAMETER :: JDXTEM = -1  ! Temperature Jacobian
    INTEGER(I4), PARAMETER :: JDXPRE = -2  ! Pressure Jacobian
    INTEGER(I4), PARAMETER :: JDXSFT = -3  ! Surface Temperature Jacobian
    INTEGER(I4), PARAMETER :: JDXSFE = -4  ! Surface Emissivity Jacobian
    INTEGER(I4), PARAMETER :: JDXVMR = 100 ! Assumed Max no.VMR profiles
    INTEGER(I4), PARAMETER :: LENJAC = 13  ! Max length of Jacobian code
!
! Note: VMR Jacobians are identified by JDX values 1:JDXVMR
!
  TYPE :: JACTYP
    LOGICAL           :: COL ! T=column retrieval, F=profile
    INTEGER(I4)       :: IAT ! Atmos level# of each Jacobian rtvl. element
    INTEGER(I4)       :: ILO ! Lower Atmos level# for perturbation=0
    INTEGER(I4)       :: ITN ! Tang.ray associated with Jacob (or 0 if none)
    INTEGER(I4)       :: IUP ! Upper Atmos level# for perturbation=0
    INTEGER(I4)       :: JDX ! Index of target species
    REAL(R4)          :: HGT ! Perturbation altitude [km]
    CHARACTER(LENJAC) :: COD ! Character code for Jacobian type
  END TYPE JACTYP
!
! GLOBAL VARIABLES
    TYPE(JACTYP), ALLOCATABLE :: JAC(:)
!
    INTEGER(I4) :: NJAC = 0    ! No. elements for Jacobian calc.
    INTEGER(I4) :: NTNJAC = 0  ! No. unptb tan.rays with defined Jacobians
    INTEGER(I4), ALLOCATABLE :: ITNJAC(:,:) ! [NTNJAC,NJAC] Indices of Jac rays
!
END MODULE JACCOM_DAT

