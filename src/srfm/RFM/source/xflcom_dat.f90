MODULE XFLCOM_DAT
!
! VERSION
!   19FEB20 AD Original. Checked.
!
! DESCRIPTION
!   Contents of .xsc file
!   To save space, only x/s tabulations which overlap the RFM required range
!   are actually saved.
!   Read in by REAXSC, and allocated
!   Used by SAVXSC, and deallocated
!   
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: XFLTYP
    INTEGER(I4) :: IOF ! Offset for data in ABCXFL
    INTEGER(I4) :: ISP ! Spectral band#
    INTEGER(I4) :: NPT ! No. data points in tabulation
    REAL(R4)    :: PRE ! Pressure [Torr] of tabulation
    REAL(R4)    :: TEM ! Temperature [K] of tabulation
    REAL(R8)    :: WN1 ! Lower wavenumber of tabulation
    REAL(R8)    :: WN2 ! Upper wavenumber of tabulation
  END TYPE XFLTYP
!
! GLOBAL VARIABLES
    INTEGER(I4)               :: NXFL      ! No. stored tabulations
    REAL(R4),     ALLOCATABLE :: ABCXFL(:) ! Absorption coefficient data
    TYPE(XFLTYP), ALLOCATABLE :: XFL(:)    ! Tabulation info
!
END MODULE XFLCOM_DAT
