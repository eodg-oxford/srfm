MODULE VRTSUM_SUB
CONTAINS
SUBROUTINE VRTSUM ( IVMR, IATM, TEM, PRE, PPA, AMT, LEN ) 
!
! VERSION
!   05JAN23 AD Checked.
!   01MAY17 AD F90 conversion of rfmnad.for. Checked.
!
! DESCRIPTION
!   Curtis-Godson integrals for vertical path segment
!   Called by VRTPTH, FLXPTH
!   Calculate CG quantities for each path = absorber & layer
!   NB: CG-pressure calculation also used in subroutine DIMPRE
!
!   Use notation Dx = x_l - x_u where x_l,x_u are x evaluated at the lower
!   and upper boundaries respectively.
!   v = vmr, p = pressure, T = temperature, e=partial pressure, u=abs.amount
!   Assumes T varies linearly with altitude in all cases.
!
!   LIN flag disabled: 
!     Assumes v varies within a layer as v = v_l(p/p_l)^b
!     where b=ln(v_u/v_l)/ln(p_u/p_l)
!        u =  1/(gM(1+b)) Dpv         
!       pu =  1/(gM(2+b)) Dp^2v      
!       eu = 1/(2gM(1+b)) Dp^2v^2    
!       Tu = 1/(gM(1+b)) ( DpTv - DTDpv/((1+b)Dlnp) ) 
!     special cases 
!       b=-1,   u = (1/gM)(v_l p_l Dlnp)     
!              eu = (1/gM)(v_l^2 p_l^2 Dlnp )
!              Tu = (1/gM)(v_l p_l (T_u+T_l)/2 Dlnp )
!       b=-2,  pu = (1/gM)(v_l p_l^2 Dlnp )
!   LIN flag enabled
!     Assumes v varies within a layer as v = v_l + a * (z - z_l)
!        u = (1/gM)(Dpv - DpDv/Dlnp)
!       pu = (1/gM)(Dp^2v - Dp^2Dv/(2Dlnp))
!        eu = (1/gM)(Dp^2v^2 - Dp^2vDv/Dlnp + Dp^2(Dv)^2/(2Dlnp)^2) )
!       Tu = (1/gM)(DpTv - (DvDpT+DTDpv)/(Dlnp) + 2DpDTDv/(Dlnp)^2 )
!
!   To avoid making assumptions about g or M, for hydrostatic equilbm assume
!                   M.g = R <T> ln(P2/P1) / ( Z2-Z1 )
!   However, to make the calculation independent of altitude scale the PRE
!   flag can be used, with the value of M.g explicitly inserted.
!
!   It is assumed that the CG pressure, part.press and temperature do not
!   vary with airmass (ie zenith angle) so calculated for first angle listed
!   in *TAN section only.
!
! VARIABLE KINDS
    USE KIND_DAT  
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE PHYADJ_DAT, ONLY: GRAVTY, WGTAIR
    USE PHYCON_DAT, ONLY: ATMB, RGAS 
    USE FLGCOM_DAT, ONLY: HYDFLG
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN)  :: IVMR ! Index of absorber
    INTEGER(I4), INTENT(IN)  :: IATM ! Atmospheric layer
    REAL(R4),    INTENT(OUT) :: TEM  ! CG Temperature [K]
    REAL(R4),    INTENT(OUT) :: PRE  ! CG Pressure [atm]
    REAL(R4),    INTENT(OUT) :: PPA  ! CG Part.Press [atm]
    REAL(R4),    INTENT(OUT) :: AMT  ! Absorb.amount [kmol/cm2]
    REAL(R4),    INTENT(OUT) :: LEN  ! Path length [km]
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: ARGMIN = TINY ( 1.0 ) ! smallest real no.
! Note: thin layers can cause problems in evaluating DPV etc, so avoid
    REAL(R4), PARAMETER :: DLPMIN = 1.0E-3 ! Min value for dlnp for expl.diff.
                                       ! 0.1%, about 7m min layer thickness!
! LOCAL VARIABLES
    REAL(R4) :: B          ! VMR profile parameter: VMR = VMR1*(P/P1)^X
    REAL(R4) :: DLNP       ! ln(p_l/p_u)
    REAL(R4) :: DP         ! p_l - p_u
    REAL(R4) :: DP2        ! p_l^2 - p_u^2
    REAL(R4) :: DP2V       ! p_l^2 v_l - p_u^2 v_u
    REAL(R4) :: DP2V2      ! p_l^2 v_l^2 - p_u^2 v_u^2
    REAL(R4) :: DPT        ! p_l T_l - p_u T_u 
    REAL(R4) :: DPTV       ! p_l T_l v_l - p_u T_u v_u
    REAL(R4) :: DPV        ! p_l v_l - p_u v_u
    REAL(R4) :: DT         ! T_l - T_u
    REAL(R4) :: DV         ! v_l - v_u
    REAL(R4) :: ECG        ! Curtis-Godson partial pressure [mb]
    REAL(R4) :: HGT1, HGT2 ! Hgts [km] at lower,upper segment bounds
    REAL(R4) :: PAVG       ! Average pressure [mb] sqrt ( p_l * p_u )
    REAL(R4) :: PCG        ! Curtis-Godson pressure [mb]
    REAL(R4) :: PRE1, PRE2 ! Pressure [mb] at low,upp segment bounds
    REAL(R4) :: PVMAX      ! Max value of ( p_l v_l or p_u v_u )
    REAL(R4) :: PVMIN      ! Min value of ( p_l v_l or p_u v_u )
    REAL(R4) :: TAVG       ! Mid-point temperature [K] of segment 
    REAL(R4) :: TCG        ! Curtis-Godson temperature [K]
    REAL(R4) :: TEM1, TEM2 ! Temperatures [K] at low,upp seg.bounds
    REAL(R4) :: USUM       ! Absorber Amount in layer [10^5*kmole/m^2]
    REAL(R4) :: VMR1, VMR2 ! Vmr [ppv] at low,upp seg.bounds
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Layer lower boundary
  HGT1 = HGTATM(IATM)
  TEM1 = TEMATM(IATM)
  PRE1 = PREATM(IATM)
  VMR1 = VMRATM(IATM,IVMR) * 1.0E-6     ! Convert [ppmv] to [ppv]
! Layer upper boundary
  HGT2 = HGTATM(IATM+1)
  TEM2 = TEMATM(IATM+1)
  PRE2 = PREATM(IATM+1)
  VMR2 = VMRATM(IATM+1,IVMR) * 1.0E-6   ! Convert [ppmv] to [ppv]
!
! Calculate summation terms for CG integrals (allow for VMR1=0 and/or VMR2=0)
  DLNP  = LOG(PRE1/PRE2)
  DP2V  = PRE1*PRE1*VMR1 - PRE2*PRE2*VMR2
  DP2V2 = PRE1*PRE1*VMR1*VMR1 - PRE2*PRE2*VMR2*VMR2
  DPTV  = PRE1*TEM1*VMR1 - PRE2*TEM2*VMR2
  DPV   = PRE1*VMR1 - PRE2*VMR2
  DT    = TEM1 - TEM2
  TAVG = 0.5 * ( TEM1 + TEM2 )
  PAVG = SQRT ( PRE1 * PRE2 )
  IF ( LINVMR(IVMR) .OR. VMR1 .LT. ARGMIN .OR. VMR2 .LT. ARGMIN &
      .OR. ABS(DPV) .LT. ARGMIN .OR. ABS(DP2V) .LT. ARGMIN ) THEN  
    DP    = PRE1 - PRE2
    DP2   = PRE1*PRE1 - PRE2*PRE2
    DPT   = PRE1*TEM1 - PRE2*TEM2
    DV    = VMR1 - VMR2
    USUM  = DPV - DP * DV / DLNP
    IF ( USUM .GT. ARGMIN .AND. DLNP .GE. DLPMIN ) THEN  ! Bug#105
      PCG = 0.5 * ( DP2V - 0.5*DP2*DV/DLNP ) / USUM
      ECG = 0.5 * ( DP2V2 - ( DP2V*DV - 0.5*DP2*DV*DV/DLNP ) / DLNP ) / USUM
      TCG = ( DPTV - ( DV*DPT + DT*DPV - 2.0*DP*DT*DV/DLNP ) / DLNP ) / USUM
    ELSE 
      PCG = PAVG
      TCG = TAVG
      ECG = 0.0
    END IF
  ELSE
    B = LOG(VMR1/VMR2) / DLNP
    IF ( ABS ( 1.0 + B ) .LE. 1.0E-2 .OR. DLNP .LT. DLPMIN ) THEN
      USUM = PRE1*VMR1*DLNP
      PCG  = DP2V / ( 2.0 + B ) / USUM
      ECG  = PRE1*VMR1
      TCG  = TAVG 
    ELSE 
      USUM = DPV / ( 1.0 + B )
      IF ( ABS ( 2.0 + B ) .LE. 1.0E-3 ) THEN
        PCG = PRE1*PRE1*VMR1*DLNP / USUM          ! fix Bug#102
      ELSE
        PCG = DP2V / ( 2.0 + B ) / USUM
      END IF
      ECG = 0.5 * DP2V2 / DPV
      TCG =  DPTV / DPV - DT / ( 1.0 + B ) / DLNP 
    END IF
  END IF
!
! Add some extra tests in case any numerical problems have been missed
  IF ( TCG .GT. MAX ( TEM1, TEM2 ) .OR. &
       TCG .LT. MIN ( TEM1, TEM2 )      ) THEN  ! eg TEM1=TEM2
    TCG = TAVG
  END IF
  IF ( PCG .GT. PRE1 .OR. PCG .LT. PRE2 ) PCG = PAVG  ! PRE1 always > PRE2
  PVMIN = MIN ( PRE1 * VMR1, PRE2 * VMR2 )
  PVMAX = MAX ( PRE1 * VMR1, PRE2 * VMR2 )
  IF ( ECG .GT. PVMAX .OR. ECG .LT. PVMIN ) ECG = SQRT(PVMIN*PVMAX) !eg ECG=0
  IF ( USUM .GT. PVMAX * DLNP .OR. USUM .LT. 0.0 ) &
    USUM = MAX ( SQRT ( PVMIN * PVMAX ) * DLNP, 0.0 ) 
!
! So far USUM is (absorber amount)*(gM) and is in [mb]. 
! Required units are kmole.cm^-2. 
! So multiply by 100 (convert mb to Pa) for SI units, then divide by 10000
! to convert m^-2 to cm^-2, net factor 0.01
!
! Remove gM factor from USUM - see notes at top of subroutine.
! Factor 10 represents net factor 0.01 * 1000 m/km (convert HGT to SI)
  IF ( HYDFLG ) THEN
    USUM = USUM * ABS ( HGT2 - HGT1 ) / RGAS / TAVG / DLNP*10.0
  ELSE
! GRAVTY and WGTAIR in phycom.inc
    USUM = USUM / ( GRAVTY * WGTAIR ) * 0.01
  END IF
!
! Convert p from mb to atm.
  TEM = TCG
  PRE = PCG / ATMB
  PPA = ECG / ATMB
  AMT = USUM
  LEN = HGT2 - HGT1
!
END SUBROUTINE VRTSUM
END MODULE VRTSUM_SUB
