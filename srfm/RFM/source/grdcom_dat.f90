MODULE GRDCOM_DAT
!
! VERSION
!   27FEB23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Irregular grid
!   Can be loaded either from .grd file by REAGRD or from irregular spectral
!   file by REASPC. 
!   Used in SPCFUL to set up fine grid, and ILSGRD for ILS convolution. 
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: GRDTYP
    INTEGER(I4) :: IFN ! Index of irreg grid point on reg fine grid
    REAL(R8)    :: WNO ! List of wavenumbers [cm-1] 
  END TYPE GRDTYP
!
! GLOBAL VARIABLES
    TYPE(GRDTYP), ALLOCATABLE :: GRD(:)
!
    LOGICAL      :: TWOGRD     ! T=grid uses two-point interpolation
    INTEGER(I4)  :: ISPGRD = 0 ! Spc.Rng index of currently loaded grid
    INTEGER(I4)  :: NGRD = 0   ! No. wavenumber points in grid
    CHARACTER(3) :: FNCGRD     ! Interpolation function
    REAL(R8)     :: WNRGRD     ! Grid spectral resolution [cm-1]
!
END MODULE GRDCOM_DAT

