MODULE CTMC32_SUB
CONTAINS
SUBROUTINE CTMC32 ( ILBL ) 
!     
! VERSION
!   15MAR23 AD Renamed from CTMH2O to CTMC32
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
!   H2O continuum MT_CKD v3.2
!   Called by SPCCTM for each path containing H2O if CTM flag enabled.
!   Calculate H2O continuum absorption across entire widemesh grid.
!
!   This version uses the MT_CKD 3.2 continuum, previous version of this 
!   subroutine using the MT_CKD v2.5 continuum is now renamed CTMC25.
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
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE H2OMTC_DAT ! MT_CKD H2O continuum data
    USE WIDCOM_DAT ! Widemesh data
    USE PHYCON_DAT, ONLY: AVOG, C2, PREREF, TEMREF
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
! Parameters for self-broadening correction
    REAL(R4), PARAMETER :: BETA = 350.0 ! Central Wno [cm-1] for correction 
    REAL(R4), PARAMETER :: F1   = 0.25 
    REAL(R4), PARAMETER :: N_S  = 6.0  
! Parameters for foreign-broadening correction
    REAL(R4), PARAMETER :: BETA1 = 57.83 
    REAL(R4), PARAMETER :: BETA2 = 630.0 
    REAL(R4), PARAMETER :: C_1 = -0.42 
    REAL(R4), PARAMETER :: C_2 = 0.3 
    REAL(R4), PARAMETER :: F0 = 0.06 
    REAL(R4), PARAMETER :: HWSQ1 = 240.0**2 
    REAL(R4), PARAMETER :: N_1 = 8.0 
    REAL(R4), PARAMETER :: N_2 = 8.0 
    REAL(R4), PARAMETER :: V0F1 = 255.67    ! Central Wno [cm-1] for correction
!
! Self Correction from RHUBC-II    mt_ckd_3.0   Nov 2016
    REAL(R4), PARAMETER :: F1_RHU = 0.08
    REAL(R4), PARAMETER :: BETA_RHU  = 40.
!
! Self-broadening modification factor 820-960@10cm-1
! XFCREV is unchanged from v1.10 to v2.5
! "Error" in original RFM implementation - changed from 0:14 to 15 in XFCREV
    REAL(R4), PARAMETER :: XFCREV(15) = (/ &
      1.003, 1.009, 1.015, 1.023, 1.029, 1.033, 1.037, 1.039, 1.040, &
      1.046, 1.036, 1.027, 1.01, 1.002,  1.00 /)
!
! New in v3.2
! Foreign correction factors from joint RHUBC-II/RHUBC-I analysis (mt_ckd_3.0).
! Note XFAC_RHU changed from (-1:61) in MTCKD code to (63) - element 3 is wno=0!
    REAL(R4), PARAMETER :: XFAC_RHU(63) = (/ &
          0.7810,0.8330,                                             &
          0.8500,0.8330,0.7810,0.7540,0.8180,                        &
          0.9140,0.9980,0.9830,0.9330,0.8850,                        &
          0.8420,0.8070,0.8000,0.8010,0.8100,                        &
          0.8090,0.8320,0.8180,0.7970,0.8240,                        &
          0.8640,0.8830,0.8830,0.8470,0.8380,                        &
          0.8660,0.9410,1.0400,1.0680,1.1410,                        &
          1.0800,1.0340,1.1550,1.0990,1.0270,                        &
          0.9500,0.8950,0.8150,0.7830,0.7700,                        &
          0.7000,0.7650,0.7750,0.8500,0.9000,                        &
          0.9050,0.9540,1.0200,1.0200,1.0250,                        &
          1.0200,1.1000,1.1250,1.1200,1.1110,                        &
          1.1370,1.1600,1.1490,1.1070,1.0640,                        &
          1.0450 /)
!
! LOCAL VARIABLES      
    INTEGER(I4) :: ICLC   ! Index of calc path segment
    INTEGER(I4) :: IQAD   ! Counter for quadratic interp.pts (1:3)
    INTEGER(I4) :: IWID   ! Counter for widemesh intervals (1:NWID)
    INTEGER(I4) :: IWD2   ! Counter for half-wide-mesh grid (0:NWID*2)
    REAL(R4)    :: C11    ! 0.5 C2/T [cm], where C2 is a radiation constant 
    REAL(R4)    :: CTWFRN ! FASCOD-correction applied to Foreign-broad CWFRN 
    REAL(R4)    :: CTWSLF ! Self-broadening interp. to Half-WM and Path Temp.
    REAL(R4)    :: CW260  ! Continuum 260K Self-Broad Coeff.interp to Half-WM 
    REAL(R4)    :: CW296  ! Continuum 296K Self-Broad Coeff.interp to Half-WM 
    REAL(R4)    :: CWFRN  ! Continuum Foreign Coeff. interp to Half-WM point.
    REAL(R4)    :: FSCAL  ! Scaling factor for foreign broad 100-600cm-1 region
    REAL(R4)    :: SFAC   ! Interpolated value of XFAC for current wavenumber
    REAL(R4)    :: VDSQ1, VDMSQ1  ! fns of Wno association with f.b.correction
    REAL(R4)    :: VF1, VF2, VMF1 ! fns of Wno associated with f.b.correction
    REAL(R4)    :: WNO    ! Wavenumber of current half-wide-mesh point
    REAL(R4)    :: XW     ! Position of half-WM pt on H2O ctm tabulation axis
    REAL(R4)    :: XX     ! Position of half-WM pt on XFCREV tabulation axis
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
! Interpolate in wavenumber for each tabulated temperature
    CW296 = LKPIDX ( XW, H2O296 )
    CW260 = LKPIDX ( XW, H2O260 )
    CWFRN = LKPIDX ( XW, H2OFRN )
!
! Interpolate Self-CTM modification factor to current wavenumber
    IF ( WNO .GT. 820.0 .AND. WNO .LT. 960.0 ) THEN  ! same as v1.10
      XX = ( WNO - 810.0 ) / 10.0 
      SFAC = LKPIDX ( XX, XFCREV )
! Self-broad. correction factor for 2000-3190@10cm-1 removed in v3.2
    ELSE
      SFAC = 1.0
    END IF
!
! Additional modification for region near 350cm-1
    SFAC = SFAC * ( 1.0 + ( F1 / ( 1.0 + ( WNO / BETA )**N_S ) ) )  
!
!Additional Factor v3.2    
    SFAC = SFAC * ( 1.0 + ( F1_RHU / ( 1.0 + ( WNO /BETA_RHU )**N_S ) ) )
!
! Exponential extrapolation between temperatures
! 2.77777E-2 = 1/(296.0 - 260.0)
    CTWSLF = SFAC * CW296 * &
            (CW260/CW296)**( 2.77777E-2 * ( 296.0 - CLC(ICLC)%TEM ) )
	    
! Correction to foreign continuum 
!     Extensive changes to the foreign continuum were made for
!     mt_ckd_2.8.  Based on measurements by Baranov and Lafferty
!     (2012) from 850-1160 cm-1 and Mondelain et al.(2014) at
!     4255 cm-1, a revised MT_CKD foreign continuum formulation
!     was derived and has been implemented in window regions >
!     4000 cm-1 (blended with previous coefficients in transition
!     regions between windows and bands). For 1800-3000 cm-1,
!     this formulation guided the spectral shape of the foreign
!     coefficients, but the values were reduced to obtain
!     agreement with measurements in this window by Baranov and
!     Lafferty (2012) and IASI measurements from 1900-2150 cm-1.
!     Coefficients in this region were derived simultaneously
!     with N2-H2O CIA coefficients and water vapor self continuum
!     coefficents.
!
! Calculating scaling factor correction for foreign broadening ~100-600cm-1
!
! And adding RHU scaling for v3.2
! Note indexing differs from MTCKD - so that index 3 is wno=0.0
    IF (WNO .LE. 600.) THEN
          XX = ( WNO +30. )/ 10.     ! XX=3 is XFAC_RHU @ 0.0 cm-1
          FSCAL = LKPIDX ( XX, XFAC_RHU )
    ELSE
          VDSQ1  = ( WNO - V0F1 )**2
          VDMSQ1 = ( WNO + V0F1 )**2
          VF1    = ( ( WNO - V0F1 ) / BETA1 )**N_1
          VMF1   = ( ( WNO + V0F1 ) / BETA1 )**N_1
          VF2    = ( WNO / BETA2 )**N_2
          FSCAL = 1.0 + ( F0 + C_1 * &
                ( ( HWSQ1 / ( VDSQ1 + HWSQ1 + VF1 ) ) +  &
                  ( HWSQ1 / ( VDMSQ1 + HWSQ1 + VMF1 ) )  ) ) &
               / ( 1.0 + C_2 * VF2 )
    ENDIF
!
    CTWFRN = CWFRN * FSCAL
!
    CTMPTH(IWD2) = CLC(ICLC)%AMT * AVOG * 1.0E-20 &
             * WNO * TANH ( C11 * WNO )  & ! radiation term in contnm.f
             * TEMREF / ( CLC(ICLC)%TEM * PREREF ) & 
     * ( CLC(ICLC)%PPA * CTWSLF + ( CLC(ICLC)%PRE - CLC(ICLC)%PPA ) * CTWFRN )
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
END SUBROUTINE CTMC32
END MODULE CTMC32_SUB
