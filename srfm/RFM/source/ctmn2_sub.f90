MODULE CTMN2_SUB
CONTAINS
SUBROUTINE CTMN2 ( ILBL ) 
!
! VERSION
!   15APR23 AD Checked.
!   16NOV17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   N2 continuum 
!   Called by SPCCTM
!   Calculates absorption across entire widemesh grid.
!
!   The continuum absorption/length s(/cm) is derived from
!
!      s(p,T) = v (p/po)^2 (To/T)^2 (v+(1-v)E(T)) Bx10^-6 exp(beta(1/Tr-1/T))
! 
!   where B and beta are tabulated functions of wavenumber, To=273K, po is
!   std.atmospheric pressure, V is the vmr of N2, Tr=296K, and E(T) is the
!   (wavenumber-independent) relative efficiency of O2 broadening to N2,
!
!      E(T) = 1.294 - 0.4545 (T/Tr)
!
!   B is tabulated in (/cm /Am^2) where 1 Amagat is the ratio of the molar
!   density to the value at STP (hence the p/po To/T terms). 
!   Therefore path absorption k = s x, where x is the path length.
!   For consistency with other RFM modules, this is rearranged to be a product
!   of the path absorber amount u (kmol/cm^2), given by u = (vpx/(RT))
!   This gives a calculation 
!
!      k(p,T) = u p To^2 R / (po T) (v-(1-v)E(T)) B 10^6 exp(beta(1/Tr - 1/T))
!
!   where p is in atm, R is J/kmol/K, po is Pa. Note that the extra 10^6 
!   scaling term for the cm^3/m^3 conversion.
!
! REFERENCES
!   LAFFERTY, W.J, A.M. SOLODOV, A. WEBER, W.B. OLSON and J-M HARTMANN
!     Infrared collision-induced absorption by N2 near 4.3um for 
!     atmospheric applications: measurements and empirical modelling.
!     Applied Optics, 35, 5911-5917 (1996).
!
!   Note that the tabulation in the paper extends from 2125-2600, but here
!   values have been added at 2120 and 2605 (0 for B, end values for BETA)
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
    INTEGER(I4), INTENT(IN) :: ILBL ! Index of LBL calc path
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: NWNO = 98       ! No.tab pts, upper Wno= 2605.0
    REAL(R4),    PARAMETER :: WNOLOW = 2120.0 ! Lower Wno [cm-1] of ctm data 
    REAL(R4),    PARAMETER :: DELWNO = 5.0    ! Wno increment [cm-1] of ctm data
    REAL(R4),    PARAMETER :: TEMREF = 296.0  ! Reference Temperature [K]
    REAL(R4),    PARAMETER :: B(NWNO) = &     ! [cm-1 Amagat-2] 
      (/ 0.000E-07, 0.445E-07, 0.522E-07, 0.646E-07, 0.775E-07, 0.903E-07, &   
         0.106E-06, 0.121E-06, 0.137E-06, 0.157E-06, 0.175E-06, 0.201E-06, &   
         0.230E-06, 0.259E-06, 0.295E-06, 0.326E-06, 0.366E-06, 0.405E-06, &  
         0.447E-06, 0.492E-06, 0.534E-06, 0.584E-06, 0.624E-06, 0.667E-06, &   
         0.714E-06, 0.726E-06, 0.754E-06, 0.784E-06, 0.809E-06, 0.842E-06, &  
         0.862E-06, 0.887E-06, 0.911E-06, 0.936E-06, 0.976E-06, 0.103E-05, &  
         0.111E-05, 0.123E-05, 0.139E-05, 0.161E-05, 0.176E-05, 0.194E-05, &  
         0.197E-05, 0.187E-05, 0.175E-05, 0.156E-05, 0.142E-05, 0.135E-05, &  
         0.132E-05, 0.129E-05, 0.129E-05, 0.129E-05, 0.130E-05, 0.132E-05, &  
         0.133E-05, 0.134E-05, 0.135E-05, 0.133E-05, 0.131E-05, 0.129E-05, &  
         0.124E-05, 0.120E-05, 0.116E-05, 0.110E-05, 0.104E-05, 0.996E-06, &    
         0.938E-06, 0.863E-06, 0.798E-06, 0.726E-06, 0.655E-06, 0.594E-06, &   
         0.535E-06, 0.474E-06, 0.424E-06, 0.377E-06, 0.333E-06, 0.296E-06, &   
         0.263E-06, 0.234E-06, 0.208E-06, 0.185E-06, 0.167E-06, 0.147E-06, &   
         0.132E-06, 0.120E-06, 0.109E-06, 0.985E-07, 0.908E-07, 0.818E-07, &   
         0.756E-07, 0.685E-07, 0.614E-07, 0.583E-07, 0.577E-07, 0.500E-07, &
         0.432E-07, 0.000E-07 /)
    REAL(R4),    PARAMETER :: BETA(NWNO) = &   !  [K]
      (/  802.0, 802.0, 761.0, 722.0, 679.0, 646.0, 609.0, 562.0, 511.0, &
          472.0, 436.0, 406.0, 377.0, 355.0, 338.0, 319.0, 299.0, 278.0, &
          255.0, 233.0, 208.0, 184.0, 149.0, 107.0,  66.0,  25.0, -13.0, &
          -49.0, -82.0,-104.0,-119.0,-130.0,-139.0,-144.0,-146.0,-146.0, &
         -147.0,-148.0,-150.0,-153.0,-160.0,-169.0,-181.0,-189.0,-195.0, &
         -200.0,-205.0,-209.0,-211.0,-210.0,-210.0,-209.0,-205.0,-199.0, &
         -190.0,-180.0,-168.0,-157.0,-143.0,-126.0,-108.0, -89.0, -63.0, &
          -32.0,   1.0,  35.0,  65.0,  95.0, 121.0, 141.0, 152.0, 161.0, &
          164.0, 164.0, 161.0, 155.0, 148.0, 143.0, 137.0, 133.0, 131.0, &
          133.0, 139.0, 150.0, 165.0, 187.0, 213.0, 248.0, 284.0, 321.0, &
          372.0, 449.0, 514.0, 569.0, 609.0, 642.0, 673.0, 673.0 /)
!
! LOCAL VARIABLES      
    INTEGER(I4) :: ICLC   ! Index of calc path segment
    INTEGER(I4) :: IQAD   ! Counter for quadratic interp.pts (1:3)
    INTEGER(I4) :: IWD2   ! Counter for half-wide-mesh grid (0:NWID*2)
    INTEGER(I4) :: IWID   ! Counter for widemesh intervals (1:NWID)
    REAL(R4)    :: AFACT  ! Absorber amount scaling for B
    REAL(R4)    :: BETAW  ! BETA term interpolated to wavenumber
    REAL(R4)    :: BW     ! B term interpolated to wavenumber
    REAL(R4)    :: EFACT  ! Factor to allow for Foreign (ie O2) broadening 
    REAL(R4)    :: TFACT  ! Temperature scaling for BETA
    REAL(R4)    :: VMR    ! N2 VMR for path
    REAL(R4)    :: XW     ! Position of half-WM point on CO2 continuum tab. axis
    REAL(R4)    :: CTMPTH(0:NWD2) ! Continuum absorption for current path
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  ICLC = IDXLBL(ILBL)
!
  AFACT = CLC(ICLC)%AMT * CLC(ICLC)%PRE / CLC(ICLC)%TEM * &
          (273.0)**2 * RGAS / 101325.0 * 1.0E6
  TFACT = 1.0/TEMREF - 1.0/CLC(ICLC)%TEM 
  VMR   = CLC(ICLC)%PPA / CLC(ICLC)%PRE
  EFACT = VMR + (1.0-VMR) * (1.294 - 0.4545 * CLC(ICLC)%TEM / TEMREF)
!
  DO IWD2 = 0, 2*NWID
    XW = 1.0 + ( SNGL ( WNOWD2(IWD2) ) - WNOLOW ) / DELWNO 
!
! Interpolate in wavenumber for each tabulated temperature
    BW = LKPIDX ( XW, B )     ! LKPIDX gives 0 if XW out of range 1:N
    BETAW = LKPIDX ( XW, BETA )
    CTMPTH(IWD2) = AFACT * BW * EFACT * EXP ( BETAW * TFACT ) 
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
END SUBROUTINE CTMN2
END MODULE CTMN2_SUB

