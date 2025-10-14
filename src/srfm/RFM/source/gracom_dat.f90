MODULE GRACOM_DAT
!
! VERSION
!   15JAN24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Atmospheric 2-D field
!   Mostly set by PRFGRA, missing values filled by INTGRA.
!   DNSGRA and RFRGRA set by ATMAUX.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: PSITOL = 1.0E-6 ! Tolerance [deg] for matching PSI
!                                            1e-6 deg = 900m on earth's surface
! GLOBAL VARIABLES
    INTEGER(I4) :: MPSI         ! Index of currently loaded profiles
    INTEGER(I4) :: NOMGRA       ! Index of ref (psi=0)  profile
    INTEGER(I4) :: NPSI = 0     ! No. of horizontal locations defined
    REAL(R4)    :: PSIATM = 0.0 ! Current angle for loading profiles
!
    LOGICAL,  ALLOCATABLE :: FLNGRA(:,:)   ! T=VT profiles set
    LOGICAL,  ALLOCATABLE :: FLPGRA(:)     ! T=pressure profile set
    LOGICAL,  ALLOCATABLE :: FLTGRA(:)     ! T=temperature profile set
    LOGICAL,  ALLOCATABLE :: FLVGRA(:,:)   ! T=vmr profiles set
    REAL(R4), ALLOCATABLE :: DNSGRA(:,:)   ! 2D Density fields
    REAL(R4), ALLOCATABLE :: DSHGRA(:,:)   ! 2D Density scale height fields
    REAL(R4), ALLOCATABLE :: EXTGRA(:,:)   ! 2D extinction field
    REAL(R4), ALLOCATABLE :: PREGRA(:,:)   ! 2D pressure field
    REAL(R4), ALLOCATABLE :: PSIGRA(:)     ! Psi values [deg] of each location
    REAL(R4), ALLOCATABLE :: RFRGRA(:,:)   ! 2D Refractivity fields
    REAL(R4), ALLOCATABLE :: TEMGRA(:,:)   ! 2D temperature field
    REAL(R4), ALLOCATABLE :: VIBGRA(:,:,:) ! 2D Vib.Temp fields
    REAL(R4), ALLOCATABLE :: VMRGRA(:,:,:) ! 2D VMR fields
!
END MODULE GRACOM_DAT
