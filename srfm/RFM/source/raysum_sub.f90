MODULE RAYSUM_SUB 
CONTAINS
SUBROUTINE RAYSUM ( NVMR, ZBOT, DSBOT, ZTOP, DSTOP, TCG, PCG, ECG, AMT, LEN )
!
! VERSION
!   23MAY24 AD Checked.
!   29APR21 AD Bug#29 - correct WGT2 for v.thin sub-layer
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   CG integrals in ray path segment
!   Called by LIMPTH
!   Profile layers are sub-divided into subsegments (according to NP) purely 
!   for the ray-tracing within this module.
!
!   Each subsegment has a lower, mid.pt (by altitude) and upper boundary denoted
!   by 1,2,3 respectively. Further quarter-points, denoted by M, may be 
!   required near the tangent point for calculating COS term.
!
!   The algorithm is derived from GENLN2, which is itself taken from FASCODE.
!   The CG integration is new, and derived from integrating a quadratic
!   fitted to three points (y1,y2,y3) at (f,0,g).
!   The integral is:
!      I = (w1*y1 + w2*y2 + w3*y3)*h    (h is the total width, h=g-f).
!   where  w1 = (g+2*f)/(6*f), 
!          w3 = (f+2*g)/(6*g)
!          w2 = -h^2/(6*f*g)
!   NB: if f=-g, ie interval h is split equally, this reduces to 
!   (w1,w2,w3) = (1,4,1)/6, which is Simpson's rule. However, it is the 
!   altitude which is split in half, not the path length, so integrating  
!   along the path the split is not at the mid-point.
!
! REFERENCE
!   GALLERY, W.O, F.X. KNEIZYS and S.A. CLOUGH
!   Air Mass Computer Program for Atmospheric Transmittance/Radiance 
!   Calculation: FSCATM.
!   Report: AFGL-TR-83-0065, Hanscom AFB, Mass. 01731.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE PHYADJ_DAT, ONLY: RADCRV ! Local radius of curvature [km]
    USE PHYCON_DAT, ONLY: AVOG,ATMB ! [molec/kmole], [mb/atm] conversions
!
! SUBROUTINES
    USE DSHVAL_FNC ! Interpolate Density Scale Height
    USE RFRVAL_FNC ! Interpolate Refractivity profile
    USE SNELL_SUB  ! Apply Snell's Law
    USE VALATM_SUB ! Interpolate values from atmospheric profiles
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN)  :: NVMR    ! No. different absorbers
    REAL(R4),    INTENT(IN)  :: ZBOT    ! Initial lowest altitude
    REAL(R8),    INTENT(IN)  :: DSBOT   ! Initial sine(zenith angle)
    REAL(R4),    INTENT(IN)  :: ZTOP    ! Final upper altitude
    REAL(R8),    INTENT(OUT) :: DSTOP   ! Final sine(zenith angle)
    REAL(R4),    INTENT(OUT) :: TCG(:)  ! CG Temperatures [K]
    REAL(R4),    INTENT(OUT) :: PCG(:)  ! CG Pressures [atm]
    REAL(R4),    INTENT(OUT) :: ECG(:)  ! CG Part.Press [atm]
    REAL(R4),    INTENT(OUT) :: AMT(:)  ! Absorb.amounts [kmol/cm2]
    REAL(R4),    INTENT(OUT) :: LEN     ! Path length [km]
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: NP = 2  ! Min.no. path segs for integration [NP]
!
! LOCAL VARIABLES
    LOGICAL  :: FINAL            ! TRUE = last segment in current path 
    REAL(R4) :: DNS2, DNS3       ! Air molecular density [/cm^3]
    REAL(R4) :: PRE2, PRE3       ! Pressure [mb]
    REAL(R4) :: PREBOT           ! Pressure [mb] at lower boundary
    REAL(R4) :: TEM2, TEM3       ! Temperatures [K]
    REAL(R4) :: TEMBOT           ! Temperature [K] at lower boundary
    REAL(R4) :: VMR2(NVMR), VMR3(NVMR) ! Absorber vmrs [ppmv]
    REAL(R4) :: Z1, Z2, Z3, ZM   ! Altitudes [km] segment bounds & midpts
!  
    REAL(R8) :: DCC                        ! n.r.sin(theta)
    REAL(R8) :: DCOS1, DCOS2, DCOS3        ! Cos of zenith angle 
    REAL(R8) :: DELS12, DELS23, DELS13     ! Path len. 1-2,2-3,1-3 [km]
    REAL(R8) :: DELX12, DELX23             ! DX increments 1-2, and 2-3
    REAL(R8) :: DELZ                       ! Vert thickness of seg [km]
    REAL(R8) :: DRAD1, DRAD2, DRAD3, DRADM ! Rad.of DZn from cent. earth
    REAL(R8) :: DRAT1, DRAT2, DRAT3, DRATM ! R/(n/n') (R=DRAD)
    REAL(R8) :: DRFR2, DRFR3, DRFRM        ! Refractivity
    REAL(R8) :: DSDX1, DSDX2, DSDX3        ! dS/dX, where S=path length
    REAL(R8) :: DSH                        ! Scale ht [km] dens or Refr
    REAL(R8) :: DSIN1, DSIN2, DSIN3, DSINM ! Sine of zenith angle 
    REAL(R8) :: DSMAX                      ! Max length of segment[km]
    REAL(R8) :: DSRAY                      ! Ray path length [km]
    REAL(R8) :: DSUMPP(NVMR)        ! Pre[mb]*VMR[ppmv]*abs amount
    REAL(R8) :: DSUMPU(NVMR)        ! Pre[mb]*abs amount
    REAL(R8) :: DSUMTU(NVMR)        ! Tem*abs amount
    REAL(R8) :: DSUMU(NVMR)         ! Absorber amount [cm-2]
    REAL(R8) :: DZ1, DZ2, DZ3, DZM  ! Altitudes [km] 
    REAL(R8) :: DZMIN               ! Min segment thickness [km]
    REAL(R8) :: DX1, DX2, DX3       ! Intermediate variables
    REAL(R8) :: DY1, DY2, DY3       ! Integrals cos near tp
    REAL(R8) :: FG1(NVMR),  FG2(NVMR),  FG3(NVMR)  ! Gas molec dens [/cm^3]
    REAL(R8) :: FP1(NVMR),  FP2(NVMR),  FP3(NVMR)  ! Gas mol.dens*P[mb/cm^3]
    REAL(R8) :: FPP1(NVMR), FPP2(NVMR), FPP3(NVMR) ! mol.dns *PP [mb/cm^3]
    REAL(R8) :: FT1(NVMR),  FT2(NVMR),  FT3(NVMR)  ! mol.dens * Temp [K/cm^3]
    REAL(R8) :: WGT1, WGT2, WGT3    ! Weights for CG integration
!
! EXECUTABLE CODE -------------------------------------------------------------
! 
  DSH = DSHVAL ( 0.5 * ( ZBOT + ZTOP ) )    ! Scale Ht for refract.
!
! Nominal path segment length based on layer thickness/NP 
! (Limit to .GE. 10 metre in case t.p. ZBOT close to ZTOP )
  DSMAX = DBLE ( MAX ( ( ZTOP - ZBOT ) / NP, 0.01 ) ) 
!
! Initialise variables at lowest altitude
  Z3 = ZBOT                   
  DZ3 = DBLE ( Z3 )           
  DRAD3 = DZ3 + RADCRV     ! Distance from centre of earth
  DSIN3 = DSBOT
!
! Set a minimum vertical thickness of a segment based on a horizontal
! distance DSMAX through a layer (from ds^2 = 2.R.dz approximately)
  DZMIN = 0.5D0 * DSMAX**2 / DRAD3
!
! Interpolate required quantities at Z3
  DRFR3 = RFRVAL ( Z3 )
!
! Calculate snells law constant R = - r/(n/n'), =0 if refractivity=0
  CALL SNELL ( .FALSE., DRAD3, DRFR3, DSH, DCC, DSIN3, DRAT3 )
!
  CALL VALATM ( Z3, TEM3, PRE3, DNS3, VMR3 ) 
  PREBOT = PRE3
  TEMBOT = TEM3
  FG3  = DBLE ( DNS3 * VMR3 ) 
  FT3  = FG3 * DBLE ( TEM3 )
  FP3  = FG3 * DBLE ( PRE3 )
  FPP3 = FP3 * DBLE ( VMR3 )
!
! Cos of initial zenith angle
  DCOS3 = DSQRT ( MAX ( 1.D0 - DSIN3**2, 0.D0 ) ) 
  DY3 = MAX ( 1.0D0 - DSIN3, 0.D0 )        ! Matches DCOS3 definition
  DX3 = DRAD3 * DCOS3                      ! Change to intermediate variable X
  DSDX3 = 1.D0/(1.D0 - DRAT3 * DSIN3**2 )  ! ds/dx
!
! Initialise layer lower boundary weighted quantities for integration
  DSRAY = 0.D0                             ! Ray length is zero at this stage
  DSUMU = 0.D0                             ! Path gas amount
  DSUMTU = 0.D0                            ! Path average temperature
  DSUMPU = 0.D0                            ! Path average pressure
  DSUMPP = 0.D0                            ! Path average partial pressure
!
  FINAL = .FALSE.
  DO WHILE ( .NOT. FINAL )                   ! Loop over sub-layers
!
! Update quantites at new sublevel 1 from old sublevel 3
    DZ1 = DZ3
    Z1 = Z3
    DRAD1 = DRAD3
    DSIN1 = DSIN3
    DCOS1 = DCOS3
    DX1 = DX3
    DRAT1 = DRAT3
    DSDX1 = DSDX3
    FG1 = FG3
    FT1 = FT3
    FP1 = FP3
    FPP1 = FPP3
!
! Calculate altitude increment for next segment
    DELZ = MAX ( DSMAX * DCOS1, DZMIN )
    DZ3 = DZ1 + DELZ
!    write (*,*) sngl(dz3), sngl(delz), ztop, sngl(dzmin)
    IF ( DZ3 .GE. DBLE(ZTOP) - DZMIN ) THEN ! top of layer reached
      DZ3 = DBLE ( ZTOP )
      DELZ = DZ3 - DZ1
      FINAL = .TRUE.
    END IF
    DZ2 = 0.5D0 * ( DZ1 + DZ3 )
    Z2 = SNGL ( DZ2 )
    Z3 = SNGL ( DZ3 )
!
! Calculate refracted path parameters at Z2
    DRAD2 = RADCRV + DZ2
    DRFR2 = RFRVAL ( Z2 )
    CALL SNELL ( .TRUE., DRAD2, DRFR2, DSH, DCC, DSIN2, DRAT2 )
!
! Interpolate required quantities at Z2 and evaluate required functions
    CALL VALATM ( Z2, TEM2, PRE2, DNS2, VMR2 ) 
    FG2  = DBLE ( DNS2 * VMR2 )
    FT2  = FG2 * DBLE ( TEM2 )
    FP2  = FG2 * DBLE ( PRE2 )
    FPP2 = FP2 * DBLE ( VMR2 ) 
!
! Calculate refracted path parameters at Z3
    DRAD3 = RADCRV + DZ3
    DRFR3 = RFRVAL ( Z3 )
    CALL SNELL ( .TRUE., DRAD3, DRFR3, DSH, DCC, DSIN3, DRAT3 )
!
! Interpolate required quantities at Z3 and evaluate required functions
    CALL VALATM ( Z3, TEM3, PRE3, DNS3, VMR3 ) 
    FG3  = DBLE ( DNS3 * VMR3 )
    FT3  = FG3 * DBLE ( TEM3 )
    FP3  = FG3 * DBLE ( PRE3 )
    FPP3 = FP3 * DBLE ( VMR3 )
!
! Need the cosine of the angle, leads to precision problems. near the tangent 
! height
    IF ( (1.D0 - DSIN2) .GT. 1.0D-5 ) THEN    ! Not too close to t.p.
      DCOS2 = DSQRT ( 1.D0 - DSIN2**2 )
      DCOS3 = DSQRT ( 1.D0 - DSIN3**2 )
    ELSE                                      ! Close to t.p.
      DY1 = DY3
      DZM = 0.5D0 * ( DZ1 + DZ2 )
      ZM = SNGL ( DZM )
      DRADM = DZM + RADCRV
      DRFRM = RFRVAL ( ZM ) 
      CALL SNELL ( .TRUE., DRADM, DRFRM, DSH, DCC, DSINM, DRATM )
      DY2 = DY1 + ( DZM - DZ1 ) / 3.D0 * (                    &
                        ( 1.D0 - DRAT1 ) * DSIN1 / DRAD1 +   &
                 4.D0 * ( 1.D0 - DRATM ) * DSINM / DRADM +   &
                        ( 1.D0 - DRAT2 ) * DSIN2 / DRAD2   )
      DCOS2 = DSQRT ( 2.D0*DY2 - DY2**2 )
!
      DZM = 0.5D0 * ( DZ2 + DZ3 )
      ZM = SNGL ( DZM ) 
      DRADM = DZM + RADCRV
      DRFRM = RFRVAL ( ZM )
      CALL SNELL ( .TRUE., DRADM, DRFRM, DSH, DCC, DSINM, DRATM )
      DY3 = DY2 + ( DZM - DZ2 ) / 3.D0 * (                    &
                        ( 1.D0 - DRAT2 ) * DSIN2 / DRAD2 +   &
                 4.D0 * ( 1.D0 - DRATM ) * DSINM / DRADM +   &
                        ( 1.D0 - DRAT3 ) * DSIN3 / DRAD3   )
      DCOS3 = DSQRT ( 2.D0 * DY3 - DY3**2 )
    ENDIF            ! End special case for calculating Cos
!
    DX2 = DRAD2 * DCOS2
    DX3 = DRAD3 * DCOS3
    DELX12 = DX2 - DX1
    DELX23 = DX3 - DX2
!
! Ray path length
    DSDX2 = 1.D0 / ( 1.D0 - DRAT2 * DSIN2**2 )
    DSDX3 = 1.D0 / ( 1.D0 - DRAT3 * DSIN3**2 )
    DELS12 = 0.5D0 * DELX12 * ( DSDX1 + DSDX2 )
    DELS23 = 0.5D0 * DELX23 * ( DSDX2 + DSDX3 )
    DELS13 = DELS12 + DELS23
    DSRAY = DSRAY + DELS13
!
! Path gas amount, Av.Temp, Av. Pres, and Av.P.Pres 
! The following integral is for a quadratic fitted to three points
! (f=-DELS12, g=DELS23 h=DELS13 from notes at top).
! Factor 1.0E5 to convert distance [km] to [cm]
    IF ( ABS ( DELS12 ) .GT. 0.001 .AND. &
         ABS ( DELS23 ) .GT. 0.001       ) THEN
      WGT1 = 1.0E5 * ( DELS13 * (2.0D0*DELS12-DELS23) / (6.0D0*DELS12) )
      WGT3 = 1.0E5 * ( DELS13 * (2.0D0*DELS23-DELS12) / (6.0D0*DELS23) )
      WGT2 = 1.0E5 * ( DELS13**3 / (6.0D0*DELS12*DELS23) )
    ELSE                            ! very thin layer
      WGT1 = 0.0
      WGT2 = 1.0E5 * DELS13                          ! Bug#29
      WGT3 = 0.0
      IF ( .NOT. FINAL ) THEN
        write (*,*) dz3, dzmin, dsmax, dels12,dels23,dels13, ztop, zbot
        STOP 'F-RAYSUM: logical error'   ! start near top of layer
      END IF 
    END IF
!
    DSUMU = DSUMU  + WGT1*FG1 + WGT2*FG2 + WGT3*FG3 
    DSUMTU= DSUMTU + WGT1*FT1 + WGT2*FT2 + WGT3*FT3 
    DSUMPU= DSUMPU + WGT1*FP1 + WGT2*FP2 + WGT3*FP3 
    DSUMPP= DSUMPP + WGT1*FPP1+ WGT2*FPP2+ WGT3*FPP3
!
  END DO
!
  WHERE ( DSUMU .GT. TINY(1.0) ) 
    TCG = SNGL ( DSUMTU / DSUMU )
    PCG = SNGL ( DSUMPU / DSUMU ) / ATMB   ! [atm]
    ECG = SNGL ( DSUMPP / DSUMU ) * 1.0E-6 / ATMB   ! [ppv] & [atm] 
  ELSEWHERE
    TCG = 0.5 * ( TEMBOT + TEM3 ) 
    PCG = SQRT ( PREBOT * PRE3 ) / ATMB
    ECG = 0.0
  ENDWHERE
  DSTOP = DSIN3
  AMT = SNGL ( DSUMU * 1.0E-6 ) / AVOG   ! [kmol/cm2] (1e-6: ppmv to ppv)
  LEN = SNGL ( DSRAY ) 
!
END SUBROUTINE RAYSUM
END MODULE RAYSUM_SUB
