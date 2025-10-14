MODULE LUTCOM_DAT
!
! VERSION
!   23JUN24 AD Checked.
!   28MAR19 AD Bug#19 Add %WNL, %WNU. Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   TAB LUT data
!   Look-Up Tables of absorption coefficient are stored in external files
!   This module keeps track of data and interpolation weights to convert
!   file axes to calculated path CG (p,T,q) coordinates.
!   Pointers allocated in LUTPTH, deallocated in SPCDAL
! 
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: LUTTYP
    LOGICAL      :: BIN    ! T=binary file, F=ASCII file
    LOGICAL      :: GHZ    ! T=GHz spectral axis, F=Wno spectral axis
    INTEGER(I4)  :: IBU    ! Index of upper wavenumber buffer
    INTEGER(I4)  :: IGS    ! Index of absorber
    INTEGER(I4)  :: LUN    ! LUN for LUT data
    INTEGER(I4)  :: NLK    ! No. interpolation pts
    INTEGER(I4)  :: NLP    ! No. calc.paths for TAB LUT
    INTEGER(I4)  :: NX     ! No. p,T,q points in record
    REAL(R8)     :: WNL    ! Lower wavenumber limit [cm-1] of LUT data
    REAL(R8)     :: WNU    ! Upper wavenumber limit [cm-1] of LUT data
    REAL(R8)     :: WNO(2) ! Wavenos [cm-1] of each buffer
    INTEGER(I4), POINTER :: ICL(:)   ! [NLP] Indices of calc paths
    INTEGER(I4), POINTER :: IDX(:,:) ! [NLK,NLP] Interpolation p,T,q points
    REAL(R4),    POINTER :: WGT(:,:) ! [NLK,NLP] Interpolation weights
    REAL(R4),    POINTER :: TAB(:,:) ! [NX,2] Abs.coeff data records from file
  END TYPE LUTTYP
!
! GLOBAL VARIABLES
    TYPE(LUTTYP), ALLOCATABLE :: LUT(:)
!  
    INTEGER(I4) :: NLUT = 0 ! No.LUT LUTs used for current spectral range
!
END MODULE LUTCOM_DAT
