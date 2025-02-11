MODULE GRASUM_SUB
CONTAINS
SUBROUTINE GRASUM ( NVMR, ZTAN, ZBOT, PBOT, TBOT, ZTOP, &
                    PTOP, TTOP, TCG, PCG, ECG, AMT, LEN, PSI )
!
! VERSION
!   01JUN24 AD Checked.
!   01MAY17 AD F90 conversion of rfmgra.for. Checked.
!
! DESCRIPTION
!   Construct ray paths through 2D atmosphere.
!   Called once by GRAPTH if GRA flag enabled.
!   Calculate CG quantities for each path, defined as a combination of a 
!   single absorber and single layer.
!   Profile layers are sub-divided into segments (according to NP) purely for 
!   the ray-tracing within this module.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE PHYADJ_DAT, ONLY: RADCRV ! Local radius of curvature [km]
    USE PHYCON_DAT, ONLY: ATMB, AVOG ! [mb]/[atm], [molec]/[kmol]
!
! SUBROUTINES
    USE RAYGRA_SUB ! Ray-tracing in a 2-D atmosphere (z,psi)
    USE VALGRA_SUB ! Interpolate values from 2D atmospheric field
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN)  :: NVMR ! No. different absorbers
    REAL(R4),    INTENT(IN)  :: ZTAN   ! Tangent altitude [km]
    REAL(R4),    INTENT(IN)  :: ZBOT   ! Initial lowest altitude [km]
    REAL(R4),    INTENT(IN)  :: PBOT   ! Initial horiz.angle [deg]
    REAL(R4),    INTENT(IN)  :: TBOT   ! Initial zenith angle [deg]
    REAL(R4),    INTENT(IN)  :: ZTOP   ! Final upper altitude [km]
    REAL(R4),    INTENT(OUT) :: PTOP   ! Final horiz.angle [deg]
    REAL(R4),    INTENT(OUT) :: TTOP   ! Final zenith angle [deg]
    REAL(R4),    INTENT(OUT) :: TCG(:) ! CG Temperatures [K]
    REAL(R4),    INTENT(OUT) :: PCG(:) ! CG Pressures [atm]
    REAL(R4),    INTENT(OUT) :: ECG(:) ! CG Part.Press [atm]
    REAL(R4),    INTENT(OUT) :: AMT(:) ! Absorb.amounts [kmol/cm2]
    REAL(R4),    INTENT(OUT) :: LEN    ! Path length [km]
    REAL(R4),    INTENT(OUT) :: PSI(:) ! Absorb.weighted mean horiz angle [deg]
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: NP = 10 ! No.of path segments/layer for CG int.
!
! LOCAL VARIABLES
    LOGICAL  :: FINAL          ! T = final segment for integration
    REAL(R4) :: AMT1(NVMR)     ! Absorber amounts [/cm^2] at segment lower bound
    REAL(R4) :: AMT2(NVMR)     ! Absorber amounts [/cm^2] at segment upper bound
    REAL(R4) :: AMTM(NVMR)     ! Absorber amounts [/cm^2] at segment mid-point
    REAL(R4) :: DNS1,DNS2,DNSM ! Air molec density [/cm^3] across segment
    REAL(R4) :: DS             ! Path length [km] of segment
    REAL(R4) :: DSHALF         ! Half of DS
    REAL(R4) :: DZNOM          ! Nominal altitude increment [km] of segment
    REAL(R4) :: P1, P2, PM     ! Horiz.angles [deg] across segment
    REAL(R4) :: PRE1,PRE2,PREM ! Pressures [mb] across segment
    REAL(R4) :: PREBOT         ! Pressure [mb] at layer lower bound
    REAL(R4) :: T1, T2, TM     ! Zenith angles [deg] across segment
    REAL(R4) :: TEM1,TEM2,TEMM ! Temperatures [K] across segment
    REAL(R4) :: TEMBOT         ! Temperature [K] at layer lower bound
    REAL(R4) :: VMR1(NVMR)     ! VMRs [ppmv] at segment lower bound
    REAL(R4) :: VMR2(NVMR)     ! VMRs [ppmv] at segment upper bound
    REAL(R4) :: VMRM(NVMR)     ! VMRs [ppmv] at segment mid-point
    REAL(R4) :: Z1, Z2, ZM     ! Altitudes [km] across segment
!
    REAL(R8) :: DSUMDS         ! Cumulative layer path length [km]
    REAL(R8) :: DSUMU(NVMR)    ! Cumulative abs. amount U in layer [molec/cm^2]
    REAL(R8) :: DSUMTU(NVMR)   ! Sum of Temp*U in layer [K*U]
    REAL(R8) :: DSUMPA(NVMR)   ! Sum of Psi*U in layer [deg*U]
    REAL(R8) :: DSUMPU(NVMR)   ! Sum of Pres*U in layer [mb*U]
    REAL(R8) :: DSUMPP(NVMR)   ! Sum of P.Pres*U in layer [mb*U/10^6]
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  T2 = TBOT
  Z2 = ZBOT
  P2 = PBOT
  CALL VALGRA ( Z2, P2, TEM2, PRE2, DNS2, VMR2 )
  TEMBOT = TEM2
  PREBOT = PRE2
! Nominal path length DS based on approx NP intervals through each layer
! (in the case of tangent layer: NP intervals if t.p. were at base of layer)
! If ZBOT=ZTAN reduces to DS = SQRT ( 2.0 * SNGL(RADCRV) * DZMAX ) 
  DZNOM = ( ZTOP - ZBOT ) / NP
  DS = SQRT ( 2.0 * SNGL(RADCRV) ) &
         * ( SQRT (ZBOT+DZNOM-ZTAN) - SQRT (ZBOT-ZTAN) )
!
! Initialise CG integrations for this layer
  DSUMDS = 0.0D0               ! Path length
  DSUMU  = 0.0D0               ! Path gas amount
  DSUMTU = 0.0D0               ! Path average temperature
  DSUMPU = 0.0D0               ! Path average pressure
  DSUMPP = 0.0D0               ! Path average partial pressure
  DSUMPA = 0.0D0               ! Path average horiz.angle
!
  FINAL = .FALSE.
!
  DO WHILE ( .NOT. FINAL )     ! Repeat for each path segment
    Z1 = Z2
    P1 = P2
    T1 = T2
    DNS1 = DNS2
    PRE1 = PRE2
    TEM1 = TEM2
    VMR1 = VMR2
    CALL RAYGRA ( Z1, P1, T1, 0, DS, Z2, P2, T2 )
    IF ( Z2 .GT. ZTOP ) THEN      ! Stop at top of layer
      Z2 = ZTOP
      CALL RAYGRA ( Z1, P1, T1, 1, DS, Z2, P2, T2 ) ! DS will be reduced
      FINAL = .TRUE.
    END IF
    CALL VALGRA ( Z2, P2, TEM2, PRE2, DNS2, VMR2 )
    DSHALF = DS/2.0
    CALL RAYGRA ( Z1, P1, T1, 0, DSHALF, ZM, PM, TM )
    CALL VALGRA ( ZM, PM, TEMM, PREM, DNSM, VMRM )
!
! Add to CG integrals using Simpson's rule (1,4,1)/6
! 0.1 = 10^5 cm/km * 10^-6 ppv/ppmv
    DSUMDS = DSUMDS + DBLE ( DS ) 
    AMT1 = DS * DNS1 * VMR1 / 6.0 * 0.1
    AMT2 = DS * DNS2 * VMR2 / 6.0 * 0.1
    AMTM = DS * DNSM * VMRM / 6.0 * 0.1 * 4.0
! For trapezoidal rule, uncomment the next 3 lines
!      amt1 = 0.5 * ds * dns1 * vmr1
!      amt2 = 0.5 * ds * dns2 * vmr2
!      amtm = 0.0 
    DSUMU  = DSUMU  + DBLE ( AMT1 + AMTM + AMT2 )
    DSUMTU = DSUMTU + DBLE ( AMT1*TEM1 + AMTM*TEMM + AMT2*TEM2 )
    DSUMPU = DSUMPU + DBLE ( AMT1*PRE1 + AMTM*PREM + AMT2*PRE2 )
    DSUMPP = DSUMPP + DBLE ( AMT1*PRE1*VMR1 + AMTM*PREM*VMRM + AMT2*PRE2*VMR2 )
    DSUMPA = DSUMPA + DBLE ( AMT1*P1 + AMTM*PM + AMT2*P2 ) 
!
  END DO
!
! Form final integrated quantites
  PTOP = P2
  TTOP = T2
  WHERE ( DSUMU .GT. TINY(1.0) ) 
    TCG = SNGL ( DSUMTU / DSUMU ) 
    PCG = SNGL ( DSUMPU / DSUMU ) / ATMB   ! [atm]
    ECG = SNGL ( DSUMPP / DSUMU ) * 1.0E-6 / ATMB  ! [ppv] & [atm]
    PSI = SNGL ( DSUMPA / DSUMU )
! If gas amount is zero return averages for pressure and temperature
  ELSEWHERE
    TCG = 0.5 * ( TEM2 + TEMBOT ) 
    PCG = SQRT( PRE2 * PREBOT ) / ATMB
    ECG = 0.0
    PSI = 0.5 * ( PTOP + PBOT ) 
  ENDWHERE
  AMT = SNGL ( DSUMU ) / AVOG !  [kmol/cm2] 
  LEN = SNGL ( DSUMDS )
!
END SUBROUTINE GRASUM
END MODULE GRASUM_SUB
