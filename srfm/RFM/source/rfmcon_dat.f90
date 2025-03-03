MODULE RFMCON_DAT
!
! VERSION
!   27JAN24 AD Checked.
!   01MAY17 AD Original. Checked
!
! DESCRIPTION
!   Constants associated with RFM calculations
! 
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
!
! GLOBAL CONSTANTS
    INTEGER(I4), PARAMETER :: DEFFIN = 2000   ! Default No.fine grid pts/cm-1 
    REAL(R4),    PARAMETER :: PCGMAX = 0.01   ! Max ln(p) dif for scaled path
    REAL(R4),    PARAMETER :: TCGMAX = 1.0    ! Max Tem dif [K] for scaled path
    REAL(R8),    PARAMETER :: FEXC   = 1.0D0  ! Window for fine mesh calcs
    REAL(R8),    PARAMETER :: FWIND  = 25.0D0 ! Window [cm-1] for Widemesh calc 
!
END MODULE RFMCON_DAT

