MODULE CTMO2_SUB
CONTAINS
SUBROUTINE CTMO2 ( ILBL ) 
!
! VERSION
!   13APR23 AD Checked.
!   16NOV17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   O2 continuum 
!   Called by SPCCTM for each path containing O2.
!   Calculates absorption across entire widemesh grid.
!
!   The continuum absorption/length s(/cm) is derived from
!
!      s(p,T) = (p/po)^2 (To/T)^2  v Bx10^-6  exp(beta(1/Tr - 1/T))
! 
!   where B and beta are tabulated functions of wavenumber, To=273K, po is
!   std.atmospheric pressure, V is the vmr of O2, Tr=296K
!   B is tabulated in (/cm /Am^2) where 1 Amagat is the ratio of the molar
!   density to the value at STP (hence the p/po To/T terms). 
!   Therefore path absorption k = s x, where x is the path length.
!   for consistency with other RFM modules, this is rearranged to be a product
!   of the path absorber amount u (kmol/cm^2), given by u = (vpx/(RT))
!   This gives a calculation 
!
!      k(p,T) = u p To^2 R / (po T) B exp(beta(1/Tr - 1/T))
!
!   where p is in atm, R is J/kmol/K, po is Pa. Note that the 10^6 scaling 
!   term for B is cancelled by the cm^3/m^3 conversion.
!
! REFERENCES
!   Thibault, F, V.Menoux, R.Le Doucen, L. Rosenmann, J-M. Hartmann and 
!   Ch. Boulet
!     Infrared collision-induced absorption by O2 near 6.4um for 
!     atmospheric applications: measurements and empirical modelling.
!     Applied Optics, 36, 563-567 (1997).
!
!   Note that the tabulation in the paper extends from 1365-1800, but here
!   values have been added at 1360 and 1805 (0 for B, end values for BETA)
!   to give a smooth transition to zero absorption.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE WIDCOM_DAT ! Widemesh data
    USE PHYCON_DAT, ONLY: RGAS ! Gas constant [J/K/kmole]
!
! SUBROUTINES
    USE LKPIDX_FNC ! Interpolate by real array index
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: ILBL ! Index of LBL calc path segment
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: NWNO   = 90     ! No. pts, upper Wno=1805.0
    REAL(R4),    PARAMETER :: WNOLOW = 1360.0 ! Lower Wno [cm-1] of ctm data 
    REAL(R4),    PARAMETER :: DELWNO = 5.0    ! Wno increment [cm-1] of ctm data
    REAL(R4),    PARAMETER :: TEMREF = 296.0  ! Reference temperature [K]
    REAL(R4),    PARAMETER :: B(NWNO) = &     ! [cm-1 Amagat-2] 
       (/ 0.000, 0.061, 0.074, 0.084, 0.096, 0.120, 0.162, 0.208, 0.246, &
          0.285, 0.314, 0.380, 0.444, 0.500, 0.571, 0.673, 0.768, 0.853, &
          0.966, 1.097, 1.214, 1.333, 1.466, 1.591, 1.693, 1.796, 1.922, &
          2.037, 2.154, 2.264, 2.375, 2.508, 2.671, 2.847, 3.066, 3.417, &
          3.828, 4.204, 4.453, 4.599, 4.528, 4.284, 3.955, 3.678, 3.477, &
          3.346, 3.290, 3.251, 3.231, 3.226, 3.212, 3.192, 3.108, 3.033, &
          2.911, 2.798, 2.646, 2.508, 2.322, 2.130, 1.928, 1.757, 1.588, &
          1.417, 1.253, 1.109, 0.990, 0.888, 0.791, 0.678, 0.587, 0.524, &
          0.464, 0.403, 0.357, 0.320, 0.290, 0.267, 0.242, 0.215, 0.182, &
          0.160, 0.146, 0.128, 0.103, 0.087, 0.081, 0.071, 0.064, 0.000  /)
    REAL(R4),    PARAMETER :: BETA(NWNO) = &  ! [K]
       (/ 467.0, 467.0, 400.0, 315.0, 379.0, 368.0, 475.0, 521.0, 531.0, &
          512.0, 442.0, 444.0, 430.0, 381.0, 335.0, 324.0, 296.0, 248.0, &
          215.0, 193.0, 158.0, 127.0, 101.0,  71.0,  31.0,  -6.0, -26.0, &
          -47.0, -63.0, -79.0, -88.0, -88.0, -87.0, -90.0, -98.0, -99.0, &
         -109.0,-134.0,-160.0,-167.0,-164.0,-158.0,-153.0,-151.0,-156.0, &
         -166.0,-168.0,-173.0,-170.0,-161.0,-145.0,-126.0,-108.0, -84.0, &
          -59.0, -29.0,   4.0,  41.0,  73.0,  97.0, 123.0, 159.0, 198.0, &
          220.0, 242.0, 256.0, 281.0, 311.0, 334.0, 319.0, 313.0, 321.0, &
          323.0, 310.0, 315.0, 320.0, 335.0, 361.0, 378.0, 373.0, 338.0, &
          319.0, 346.0, 322.0, 291.0, 290.0, 350.0, 371.0, 504.0, 504.0  /)
!
! LOCAL VARIABLES      
    INTEGER(I4) :: ICLC   ! Index of calc path segment
    INTEGER(I4) :: IQAD   ! Counter for quadratic interp.pts (1:3)
    INTEGER(I4) :: IWID   ! Counter for widemesh intervals (1:NWID)
    INTEGER(I4) :: IWD2   ! Counter for half-wide-mesh grid (0:NWID*2)
    REAL(R4)    :: AFACT  ! Absorber amount scaling for B
    REAL(R4)    :: BETAW  ! BETA term interpolated to wavenumber
    REAL(R4)    :: BW     ! B term interpolated to wavenumber
    REAL(R4)    :: TFACT  ! Temperature scaling for BETA
    REAL(R4)    :: XW     ! Position of half-WM point on O2 ctm tab. axis
    REAL(R4)    :: CTMPTH(0:NWD2) ! Continuum absorption for current path
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  ICLC = IDXLBL(ILBL)
!
  AFACT = CLC(ICLC)%AMT * CLC(ICLC)%PRE / CLC(ICLC)%TEM * &
          (273.0)**2 * RGAS / 101325.0 
  TFACT = 1.0 / TEMREF - 1.0 / CLC(ICLC)%TEM
!
  DO IWD2 = 0, 2*NWID
    XW = 1.0 + ( SNGL ( WNOWD2(IWD2) ) - WNOLOW ) / DELWNO 
!
! Interpolate in wavenumber for each tabulated temperature
    BW = LKPIDX ( XW, B )       ! LKPIDX returns 0 if XW outside range 1:N
    BETAW = LKPIDX ( XW, BETA ) 
    CTMPTH(IWD2) = AFACT * BW * EXP ( BETAW * TFACT ) 
  END DO
!
  DO IWID = 1, NWID
    IWD2 = 2 * IWID - 3       ! -1, 1, 3, ...
    DO IQAD = 1, 3
      IWD2 = IWD2 + 1         ! 0,1,2,  2,3,4,  4,5,6
      ABSWID(IQAD,IWID,ILBL) = ABSWID(IQAD,IWID,ILBL) + CTMPTH(IWD2)
      IF ( USECNT ) &
        CNTWID(IQAD,IWID,ILBL) = CNTWID(IQAD,IWID,ILBL) + CTMPTH(IWD2)
    END DO
  END DO
!
END SUBROUTINE CTMO2
END MODULE CTMO2_SUB

