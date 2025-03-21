MODULE TABCOM_DAT
!
! VERSION
!   10JUN24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Axes for tabulated absorption coefficients 
!   Also used when reading in these files (LUT option).
!   Loaded by DIMFIL, DIMGRD, TABPTH. Used by TABPTH, OPNTAB
!   V1TAB,V2TAB,NVTAB only used when TAB used as an input file with LUT option
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL VARIABLES
    LOGICAL     :: BINTAB ! T=output binary files, F=ASCII files
    LOGICAL     :: OFFTAB ! T=apply temperature offset profile
    INTEGER(I4) :: LUNTAB ! LUN for first .tab file 
    INTEGER(I4) :: NATAB  ! Total number of abs.coeff values
    INTEGER(I4) :: NPTAB  ! No. pressure points
    INTEGER(I4) :: NQTAB  ! No. VMR points
    INTEGER(I4) :: NTAB   ! No. different molec. tab files 
    INTEGER(I4) :: NTTAB  ! No. temperature points
    INTEGER(I4) :: NVTAB  ! No. spectral points
    REAL(R8)    :: DVTAB  ! Nominal wno. grid spacing for .tab file
    REAL(R8)    :: GHZTAB ! Scaling for output spectral axis (1=cm-1,c/1e7=GHz)
    REAL(R8)    :: V1TAB  ! Lower wavenumber limit
    REAL(R8)    :: V2TAB  ! Upper wavenumber limit
    INTEGER(I4), ALLOCATABLE :: IXORG(:)    ! Original pts if subsampling axes
    REAL(R4),    ALLOCATABLE :: PAXTAB(:)   ! [NPTAB] Irreg p-axis values [mb]
    REAL(R4),    ALLOCATABLE :: QAXTAB(:)   ! [NQTAB] Irreg q-axis values [%]
    REAL(R4),    ALLOCATABLE :: TAXTAB(:)   ! [NTTAB] Irreg T-axis values [K]
    REAL(R4),    ALLOCATABLE :: TEMTAB(:)   ! [NPTAB] Atm. T prof [K] on p-axis
    REAL(R4),    ALLOCATABLE :: VMRTAB(:,:) ! [NPTAB,NTAB] Atm. VMR profs[ppmv] 
!
END MODULE TABCOM_DAT
