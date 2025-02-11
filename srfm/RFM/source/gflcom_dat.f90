MODULE GFLCOM_DAT
!
! VERSION
!   29OCT23 AD Checked.
!   01MAY17 AD F90 Conversion. Checked.
!
! DESCRIPTION
!   Irregular grid files
!   Loaded by ADDGFL.
!
! VARIABLE KINDS
    USE KIND_DAT   
!
! GLOBAL DATA
    USE NAMCOM_DAT, ONLY: LENNAM ! Max length of any filename
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: GFLTYP
    LOGICAL           :: BIN ! T=binary file, F=text file
    LOGICAL           :: GRD ! T=.grd file, F=spectral file
    REAL(R8)          :: WMN ! Lowest Wavenumber [cm-1] in file
    REAL(R8)          :: WMX ! Highest Wavenumber [cm-1] in file
    CHARACTER(LENNAM) :: NAM ! Names of GRD files
  END TYPE GFLTYP
!
! GLOBAL VARIABLES
    TYPE(GFLTYP), ALLOCATABLE :: GFL(:)
!
    INTEGER(I4) :: NGFL = 0 ! No. of GRD files to be used
!
END MODULE GFLCOM_DAT
