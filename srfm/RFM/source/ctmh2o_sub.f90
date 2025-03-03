MODULE CTMH2O_SUB
CONTAINS
SUBROUTINE CTMH2O ( ILBL ) 
!     
! VERSION
!   15MAR23 AD Rewritten/simplified for MT_CKD v4.1 continuum. Checked.
!   22JUL19 AD Rearrange to C11 = SNGL(...) to avoid compilation warnings. Checked.
!   31JAN19 AD Fix Bug#15 - indexing of ABSWID with ILBL instead of ICLC
!   02NOV18 AD Incorporated into distributed version of RFM
!              Correct XX evaluation for indexing XFCREV
!   14SEP18 KPS further bug fixes to get indexing of XFCREV and XFAC_RHU 
!               correct for RFM
!   25MAY18 KPS Bugfix
!   17MAY18 KPS new routine for MT_CKD3.2 - between ***  &&&
!   01MAY17 AD F90 conversion
  
! 
! DESCRIPTION    
!   H2O continuum MT_CKD v4.1
!   Called by SPCCTM for each path containing H2O if CTM flag enabled.
!   Calculate H2O continuum absorption across entire widemesh grid.
!
!   This version uses the MT_CKD v4.1 continuum, previous version of this 
!   subroutine using the MT_CKD v3.2 continuum is now renamed CTMC32.
!     
! REFERENCES
!   http://rtweb.aer.com/
!   Mlawer, M.J., D.C. Tobin, and S.A. Clough, 
!   A Revised Perspective on the Water Vapor Continuum:  The MT_CKD Model, 
!   in preparation for JQSRT, 2003.
!
!   The MT_CKD_1.1 differs from 1.0 in having an extra fudge factor applied
!   to the foreign broadening between approx 100-600 cm-1. 
!
!   The MT_CKD_2.5 differs from 2.5 in having different fudge factors applied
!   to both self and foreign broadening.
! 
!   The MT_CKD_3.2 has additional foreign broadening fudge factors plus new 
!   tabulated coefficients that are incorporated in H2OMTC_DAT, with older
!   v2.5 coefficients now in H2OC25_DAT 
! 
!   The MT_CKD_4.1 has temperature dependence for self-broadening via an 
!   array of coefficients. No fudge factors. Remove 1e-20 scaling.
!
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE H2OCTM_DAT ! MT_CKD H2O continuum data
    USE WIDCOM_DAT ! Widemesh data
    USE PHYCON_DAT, ONLY: ATMB, AVOG, C2
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
!
! LOCAL VARIABLES      
    INTEGER(I4) :: ICLC   ! Index of calc path segment
    INTEGER(I4) :: IQAD   ! Counter for quadratic interp.pts (1:3)
    INTEGER(I4) :: IWID   ! Counter for widemesh intervals (1:NWID)
    INTEGER(I4) :: IWD2   ! Counter for half-wide-mesh grid (0:NWID*2)
    REAL(R4)    :: C11    ! 0.5 C2/T [cm], where C2 is a radiation constant 
    REAL(R4)    :: CTWSLF ! Self-broadening interp. to Half-WM and T-adjusted
    REAL(R4)    :: CWFRN  ! Continuum Foreign Coeff. interp to Half-WM point.
    REAL(R4)    :: CWSLF  ! Continuum Self-Broad Coeff.interp to Half-WM 
    REAL(R4)    :: TCOEFF ! Self.broad. Temperature coeff interpolated to Half-WM.
    REAL(R4)    :: WNO    ! Wavenumber of current half-wide-mesh point
    REAL(R4)    :: XW     ! Position of half-WM pt on H2O ctm tabulation axis
    REAL(R4)    :: CTMPTH(0:NWD2) ! Continuum absorption for current path
!                                           
! EXECUTABLE CODE -------------------------------------------------------------
!
  ICLC = IDXLBL(ILBL)
!
  C11 = SNGL ( 0.5 * C2 / CLC(ICLC)%TEM )
!  
  CTMPTH = 0.0
!
  DO IWD2 = 0, NWD2
    WNO = SNGL ( WNOWD2(IWD2) )
    XW = 1.0 + ( WNO - WNOLOW ) / DELWNO 
!
! Interpolate in wavenumber 
    CWSLF = LKPIDX ( XW, H2OSLF )
    CWFRN = LKPIDX ( XW, H2OFRN )
    TCOEFF = LKPIDX ( XW, TCOSLF )
!
! Temperature adjustment to self-broadening
    CTWSLF = CWSLF * ( TREF/CLC(ICLC)%TEM )**TCOEFF
!
    CTMPTH(IWD2) = CLC(ICLC)%AMT * AVOG  & ! convert AMT from kmol to molecules
             * WNO * TANH ( C11 * WNO )  & ! radiation term in contnm.f
             * TREF / ( CLC(ICLC)%TEM * PREF / ATMB ) &   ! density scaling
     * ( CLC(ICLC)%PPA * CTWSLF + ( CLC(ICLC)%PRE - CLC(ICLC)%PPA ) * CWFRN )
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
END SUBROUTINE CTMH2O
END MODULE CTMH2O_SUB
