MODULE WIDCOM_DAT
!
! VERSION
!   09APR23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Widemesh data
!   Initialised in SPCWID
!   Assumes NWID widemesh intervals, indexed 1:NWID,
!   NWID+1 widemesh grid points, indexed 0:NWID.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL VARIABLES
    LOGICAL     :: USECNT    ! T=NTE calc, so CNTWID defined
    INTEGER(I4) :: NLBL = 0  ! No. LBL calculated path segments
    INTEGER(I4) :: NWD2      ! No. widemesh half grid points (=NWID*2)
    INTEGER(I4) :: NWID      ! No. of wide mesh intervals 
    REAL(R8)    :: DELWID    ! Width of Widemesh intervals
    REAL(R8)    :: WN1WID    ! WNO at lower limit of widemesh grid
    REAL(R8)    :: WN2WID    ! WNO at upper limit of widemesh grid
    REAL(R8)    :: WNLWID    ! lowest WNO required for any lines 
    REAL(R8)    :: WNUWID    ! highest WNO required for any lines 
    INTEGER(I4), ALLOCATABLE :: IDXLBL(:)     ! [NLBL] Indices of CLC pths
    REAL(R4),    ALLOCATABLE :: ABSWID(:,:,:) ! Widemesh abs.
    REAL(R4),    ALLOCATABLE :: CNTWID(:,:,:) ! Widemesh nonlte factor
    REAL(R8),    ALLOCATABLE :: WNOWID(:)     ! [0:NWID] Wide mesh pts [/cm] 
    REAL(R8),    ALLOCATABLE :: WNOWD2(:)     ! [0:NWD2] wm parab. point
!
END MODULE WIDCOM_DAT
