MODULE CTMC25_SUB
CONTAINS
SUBROUTINE CTMC25 ( ILBL ) 
!     
! VERSION
!   14SEP24 AD Checked.
!   22JUL19 AD Rearrange to C11 = SNGL(...) to avoid compilation warnings. Checked.
!   02NOV18 AD Renamed from CTMH2O to CTMC25.
!   16NOV17 AD F90 conversion
! 
! DESCRIPTION    
!   H2O continuum MT_CKD v2.5
!   Called by SPCCTM for each path containing H2O if CTM flag enabled.
!   Calculate H2O continuum absorption across entire widemesh grid.
!
!   This version uses the MT_CKD continuum, previous version of this 
!   subroutine using the CKD continuum is now renamed CTMCKD.
!   v2.5 algorithms extracted from AER program contnm.f 
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
!   The MT_CKD_2.5 differs from 1.1 in having different fudge factors applied
!   to both self and foreign broadening.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE H2OC25_DAT ! MT_CKD v2.5 H2O continuum data
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
! Self-broadening modification factor 820-960@10cm-1
! XFCREV is unchanged from v1.10 to v2.5
    REAL(R4), PARAMETER :: XFCREV(15) = (/ &
      1.003, 1.009, 1.015, 1.023, 1.029, 1.033, 1.037, 1.039, 1.040, &
      1.046, 1.036, 1.027, 1.01, 1.002,  1.00 /)
!
! Self-broad. factor 2000-3190@10cm-1
! XFCRV1 is new with v2.5
    REAL(R4), PARAMETER :: XFCRV1(120) = (/ &
      1.000, 1.040, 1.080, 1.120, 1.160, 1.200, 1.240, 1.280, 1.318, &
      1.357, 1.404, 1.453, 1.499, 1.553, 1.608, 1.674, 1.746, 1.818, &
      1.899, 1.984, 2.078, 2.174, 2.276, 2.385, 2.502, 2.624, 2.747, &
      2.883, 3.018, 3.170, 3.321, 3.473, 3.635, 3.803, 3.974, 4.144, &
      4.327, 4.500, 4.703, 4.887, 5.102, 5.286, 5.498, 5.701, 5.935, &
      6.155, 6.405, 6.633, 6.892, 7.115, 7.397, 7.650, 7.917, 8.177, &
      8.437, 8.704, 8.953, 9.192, 9.428, 9.644, 9.821, 9.954, 10.11, &
      10.17, 10.21, 10.26, 10.29, 10.28, 10.26, 10.20, 10.15, 10.16, &
      10.25, 10.02, 9.965, 10.01, 9.934, 9.847, 9.744, 9.566, 9.436, &
      9.181, 8.872, 8.547, 8.155, 7.730, 7.261, 6.777, 6.271, 5.807, &
      5.313, 4.845, 4.444, 4.074, 3.677, 3.362, 3.087, 2.826, 2.615, &
      2.385, 2.238, 2.148, 1.979, 1.939, 1.773, 1.696, 1.642, 1.569, &
      1.510, 1.474, 1.425, 1.375, 1.322, 1.272, 1.230, 1.180, 1.130, &
      1.080, 1.040, 1.000 /)
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
      XX = ( WNO - 810.0 ) / 10.0   ! XX=1 is XFCREV(1) @ 820cm-1
      SFAC = LKPIDX ( XX, XFCREV )
    ELSE IF ( WNO .GT. 2000.0 .AND. WNO .LT. 3190.0 ) THEN
      XX = ( WNO - 1990.0 ) / 10.0
      SFAC = LKPIDX ( XX, XFCRV1 )
    ELSE
      SFAC = 1.0
    END IF
!
! Additional modification for region near 350cm-1
    SFAC = SFAC * ( 1.0 + ( F1 / ( 1.0 + ( WNO / BETA )**N_S ) ) )    
!
! Exponential extrapolation between temperatures
! 2.77777E-2 = 1/(296.0 - 260.0)
    CTWSLF = SFAC * CW296 * &
             (CW260/CW296)**( 2.77777E-2 * ( 296.0 - CLC(ICLC)%TEM ) )
!
! Correction to foreign continuum 
! Calculating scaling factor correction for foreign broadening ~100-600cm-1
    VDSQ1  = ( WNO - V0F1 )**2
    VDMSQ1 = ( WNO + V0F1 )**2
    VF1    = ( ( WNO - V0F1 ) / BETA1 )**N_1
    VMF1   = ( ( WNO + V0F1 ) / BETA1 )**N_1
    VF2    = ( WNO / BETA2 )**N_2
    FSCAL = 1.0 + ( F0 + C_1 * &
                ( ( HWSQ1 / ( VDSQ1 + HWSQ1 + VF1 ) ) +  &
                  ( HWSQ1 / ( VDMSQ1 + HWSQ1 + VMF1 ) )  ) ) &
               / ( 1.0 + C_2 * VF2 )
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
END SUBROUTINE CTMC25
END MODULE CTMC25_SUB
