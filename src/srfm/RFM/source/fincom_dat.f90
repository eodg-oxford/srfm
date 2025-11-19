MODULE FINCOM_DAT
!
! VERSION
!   02JUN24 AD Checked
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Finemesh data
!   Common variables associated with fine mesh calculation.
!   Assumes NFIN fine intervals, indexed 1:NFIN for particular IWID.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE RFMCON_DAT, ONLY: DEFFIN ! Default No.fine grid pts/cm-1
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL VARIABLES
    INTEGER(I4) :: IF1FIN ! Index of first pt on full fine grid
    INTEGER(I4) :: IF2FIN ! Index of last pt on full fine grid
    INTEGER(I4) :: NFIN   ! No. of fine mesh grid points
    INTEGER(I4) :: NOMFIN = DEFFIN ! Nominal No.pts/cm-1 for fine mesh
    REAL(R8)    :: WN1FIN ! Lower Wno at widemesh bounds
    REAL(R8)    :: WN2FIN ! Upper Wno at widemesh bounds
    REAL(R8)    :: WNLFIN ! Lower Wno for fine mesh calc.
    REAL(R8)    :: WNUFIN ! Upper Wno for fine mesh calc.
!
    REAL(R4), ALLOCATABLE :: ABSFIN(:,:) ! [NFIN,NCLC] Fine mesh absorption 
    REAL(R4), ALLOCATABLE :: CNTFIN(:,:) ! [NFIN,NCLC] Fine mesh non-LTE factor 
    REAL(R8), ALLOCATABLE :: WNOFIN(:)   ! [NFIN] Fine mesh Wno points [/cm] 
!
END MODULE FINCOM_DAT

