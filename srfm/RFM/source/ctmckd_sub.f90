MODULE CTMCKD_SUB
CONTAINS
SUBROUTINE CTMCKD ( ILBL ) 
!
! VERSION
!   15SEP24 AD Checked.
!   22JUL19 AD Rearrange to C11 = SNGL(...) to avoid compilation warnings. Checked.
!   16NOV17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   CKD H2O continuum
!   Called by SPCCTM if USECKD flag set TRUE (default: FALSE)
!   Calculate CKD H2O continuum absorption across entire widemesh grid.
!   NB: The RFM default is to use the new MT_CKD continuum in CTMH2O.
!     
!   This uses CKD v2.41. To use CKD v2.1 set the flag USE241 = FALSE
!     
! REFERENCES
!   CLOUGH S.A.,KNEIZYS F.X.,DAVIES R.,GAMACHE R., TIPPING R. (1980)
!     Theoretical line shape for H2O vapour: Application to the continuum
!     Atmospheric Water Vapour, Eds. A.Deepak,T.D.Wilkeson, L.H.Ruhnke, 
!     Academic Press, New York
!
!   CLOUGH S.A.,KNEIZYS F.X.,ROTHMAN L.S.,GALLERY W.O. (1981), 
!     Atmospheric spectral transmittance and radiance: FASCOD1B
!     SPIE 277 Atmospheric Transmission  152-166
!
!   CLOUGH S.A.,KNEIZYS F.X.,ROTHMAN L.S.,ANDERSON G.P.,SHETTLE E.P. (1987)
!     Current issues in infrared atmospheric transparency
!     International meeting on Atmospheric Transparency for Satellite 
!     Applications, 15-19 Sept. 1986 Capri, Italy. Ed. G.V. Silvestrini.CUEN.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE H2OCKD_DAT ! CKD H2O continuum tabulation data 
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
    LOGICAL,  PARAMETER :: USE241 = .TRUE. ! T= Use v2.41, F= use v2.1
    REAL(R4), PARAMETER :: XFAC(0:50) = & ! Self-broad modific.700-1200cm-1
       (/ 1.00000,1.01792,1.03767,1.05749,1.07730,1.09708, &
          1.10489,1.11268,1.12047,1.12822,1.13597,1.14367, &
          1.15135,1.15904,1.16669,1.17431,1.18786,1.20134, &
          1.21479,1.22821,1.24158,1.26580,1.28991,1.28295, &
          1.27600,1.26896,1.25550,1.24213,1.22879,1.21560, &
          1.20230,1.18162,1.16112,1.14063,1.12016,1.10195, &
          1.09207,1.08622,1.08105,1.07765,1.07398,1.06620, &
          1.05791,1.04905,1.03976,1.02981,1.00985,1.00000, &
          1.00000,1.00000,1.00000 /)
!
! LOCAL VARIABLES      
    INTEGER(I4) :: ICLC   ! Index of calc path
    INTEGER(I4) :: IQAD   ! Counter for quadratic interp.pts (1:3)
    INTEGER(I4) :: IWD2   ! Counter for half-wide-mesh grid (0:NWID*2)
    INTEGER(I4) :: IWID   ! Counter for Widemesh intervals
    REAL(R4)    :: A1,A2,A3 ! Factors in computing absorption coefficient
    REAL(R4)    :: C11    ! 0.5 C2/T [cm], where C2 is a radiation constant 
    REAL(R4)    :: CTWFRN ! FASCOD-correction applied to Foreign-broad CWFRN 
    REAL(R4)    :: CTWSLF ! Self-broad interp. to Half-WM and Path Temperature.
    REAL(R4)    :: CW260  ! Continuum 260K Self-Broad Coeff.interp to Half-WM 
    REAL(R4)    :: CW296  ! Continuum 296K Self-Broad Coeff.interp to Half-WM 
    REAL(R4)    :: CWFRN  ! Continuum Foreign Coeff. interpolated to Half-WM pt.
    REAL(R4)    :: SFAC   ! Interpolated value of XFAC for current wavenumber
    REAL(R4)    :: V2     ! (Wavenumber-offset)**2
    REAL(R4)    :: WNO    ! Wavenumber of current half-wide-mesh point
    REAL(R4)    :: XW     ! Position of half-WM point on H2O ctm tabulation axis
    REAL(R4)    :: XX     ! Position of half-WM point on XFAC tabulation axis
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
    IF ( WNO .LE. 700.0 .OR. WNO .GE. 1200.0 ) THEN
      SFAC = 1.0
    ELSE
      XX = ( WNO - 700.0 ) / 10.0 
      SFAC = LKPIDX ( XX, XFAC ) 
    END IF
!
! Clough et al. use exponential extrapolation between temperatures
! 2.77777E-2 = 1/(296.0 - 260.0)
    CTWSLF = SFAC * CW296 * &
             (CW260/CW296)**( 2.77777E-2 * ( 296.0 - CLC(ICLC)%TEM ) )
!
! FASCOD correction to self continuum (1.9.85); factor of 0.7667 at 1050 cm-1
! VOS2=1050, FACTRS2=-0.2333, HWSQ2 = 200**2
    V2 = ( WNO - 1050.0 )**2
    CTWSLF = CTWSLF * ( 1.0 - 0.2333 * 4.0E4 / ( V2 + 4.0E4 ) )
!
! CKD 2.1 Additional correction terms for self-continuum
! V0S3=1310, FACTRS3=0.15, HWSQ3=120**2, BETAS3=5E-6
    V2 = ( WNO - 1310.0 )**2
    CTWSLF = CTWSLF * ( 1.0 - 0.15 * 14400.0 / &
            ( V2 * ( 1.0 + 5.0E-6 * V2 ) + 14400.0 ) )
!
! CKD 2.41 further correction (mostly affecting < 100cm-1)
! V0S1=0, FACTRS1=0.688 (CKD 2.2 used FACTRS1=0.3), HWSQ1=100**2, BETAS1=1E-4
    IF ( USE241 ) THEN
      V2 = WNO**2
      CTWSLF = CTWSLF * ( 1.0 + 0.688 * 1.0E4 / &
              ( V2 * ( 1.0 + 1.0E-4 * V2 ) + 1.0E4 ) )
    END IF
!
! Correction to foreign continuum 
!
    CTWFRN = CWFRN
! V0F=1130, FACTRF=-0.97, HWSQF=330**2=1.089E5, BETAF=8E-11, 
    V2 = ( WNO - 1130.0 )**2
    CTWFRN = CTWFRN * ( 1.0 - 0.97 * 1.089E5 / &
             ( V2 * ( 1.0 + 8.0E-11 * V2 * V2 ) + 1.089E5 ) )
!
    IF ( USE241 ) THEN
! CKD 2.41 correction VOF1=350.0, FACTRF1=-0.7, HWSQF1=200**2, BETAF1=5E-9 
      V2 = ( WNO - 350.0 )**2
      CTWFRN = CTWFRN * ( 1.0 - 0.7 * 4.0E4 / &
              ( V2 * ( 1.0 + 5.0E-9 * V2 * V2 ) + 4.0E4 ) )
! CKD 2.41 correction V0F1A=630.0, FACTRF1a=0.75, HWSQF1a=65**2, BETAF1a=2E-8
      V2 = ( WNO - 630.0 )**2
      CTWFRN = CTWFRN * ( 1.0 + 0.75 * 4225.0 / &
               ( V2 * ( 1.0 + 2.0E-8 * V2 * V2 ) + 4225.0 ) )
! CKD 2.41 correction V0F3=1975.0, FACTRF3=-0.65, HWSQF3=250**2, BETAF3=5E-6
      V2 = ( WNO - 1975.0 )**2
      CTWFRN = CTWFRN * ( 1.0 - 0.65 * 62500.0 / &
              ( V2 * ( 1.0 + 5.0E-6 * V2 ) + 62500.0 ) )     ! NB V^4, not V^6
    ELSE
! CKD 2.1 correction VOF2=1900, FACTRF2=-0.6, HWSQF2=150**2, BETAF2=3E-6
      V2 = ( WNO - 1900.0 )**2
      CTWFRN = CTWFRN * ( 1.0 - 0.6 * 22500.0 / &
              ( V2 * ( 1.0 + 3.0E-6 * V2 ) + 22500.0 ) )     ! NB V^4, not V^6
    END IF
!
    A1 = WNO * CLC(ICLC)%AMT * AVOG * TANH ( C11 * WNO )
    A2 = TEMREF / ( CLC(ICLC)%TEM * PREREF )
    A3 = 1.0E-20 * ( CLC(ICLC)%PPA * CTWSLF + &
                     ( CLC(ICLC)%PRE - CLC(ICLC)%PPA ) * CTWFRN )
    CTMPTH(IWD2) = A1 * A2 * A3          
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
END SUBROUTINE CTMCKD
END MODULE CTMCKD_SUB
