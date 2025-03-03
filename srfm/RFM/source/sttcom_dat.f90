MODULE STTCOM_DAT
!
! VERSION
!   03MAR23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Widemesh statistics 
!   Associated with WID flag.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: STTTYP
    INTEGER(I4) :: ICT  ! 1=used ctm data for this interval, 0=no data
    INTEGER(I4) :: IXS  ! 1=used x/s data for this interval, 0=no data
    INTEGER(I4) :: NLF  ! No. lines contributing to finemesh interval
    INTEGER(I4) :: NLI  ! No. lines in widemesh interval
    INTEGER(I4) :: NLW  ! No. lines contributing to widemesh interval
    INTEGER(I4) :: NQD  ! No. paths using quadratic interpolation
  END TYPE STTTYP
!
! GLOBAL VARIABLES
    TYPE(STTTYP), ALLOCATABLE :: STT(:,:) ! [NSTW,NSTG] = [NWID,NGAS]
!
    INTEGER(I4) :: NSTW = 0 ! No. widemesh intervals
    INTEGER(I4) :: NSTG = 0 ! No. different absorbers
    REAL(R8)    :: WN1STT   ! Lower Wavenumber [cm-1] of widemesh grid
    REAL(R8)    :: WNDSTT   ! Wavenumber increment [cm-1] of widemesh grid
    REAL(R8)    :: WNLSTT   ! Lower output limit of current spectral range
    REAL(R8)    :: WNUSTT   ! Upper output limit of current spectral range
    INTEGER(I4), ALLOCATABLE :: NINSTT(:) ! [NSTG] No lines inside spc.range
    INTEGER(I4), ALLOCATABLE :: NEXSTT(:) ! [NSTG] No lines outside spc.range
!
END MODULE STTCOM_DAT
