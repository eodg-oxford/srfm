MODULE SVDCOM_DAT
!
! VERSION
!   09JUN24 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   SVD-compressed LUT data
!   Pointers ICL, IDX, WGT allocated in SVDPTH; K, U in LUTSVD
!   deallocated in SPCDAL
! 
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: SVDTYP
    INTEGER(I4)  :: IGS ! Absorber index for SVD-LUT file
    INTEGER(I4)  :: IGO ! Offset for U Wavenumber axis rel. to full fine grid
    INTEGER(I4)  :: NL  ! No. basis vectors to use (may be less than array size)
    INTEGER(I4)  :: NLK ! No. p,T interpolation points (4)
    INTEGER(I4)  :: NSP ! No. calc paths using SVD-LUT
    INTEGER(I4)  :: NX  ! Combined no. (-lnp,T) coordinates
    CHARACTER(3) :: TAB ! Tabulation function      
    REAL(R4),    POINTER :: U(:,:)   ! [NG,NL] U-matrix, NG x NL
    REAL(R4),    POINTER :: K(:,:)   ! [NL,NX] K-matrix, NL x NX
    INTEGER(I4), POINTER :: ICL(:)   ! [NSP] Index of path in CLCCOM
    INTEGER(I4), POINTER :: IDX(:,:) ! [NLK,NSP] ! p,T interpolation points
    REAL(R4),    POINTER :: WGT(:,:) ! [NLK,NSP] interpolation weights
  END TYPE SVDTYP
!
! GLOBAL VARIABLES
    TYPE(SVDTYP), ALLOCATABLE :: SVD(:)
!  
    INTEGER(I4) :: NSVD = 0 ! No.SVD LUTs used for current spectral range
!
END MODULE SVDCOM_DAT
