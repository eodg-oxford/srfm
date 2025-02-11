MODULE NTECLC_SUB
CONTAINS
SUBROUTINE NTECLC ( PRE, TEM, GAMMA, ANLTE, CNLTE, SQ )
!
! VERSION
!   22OCT23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Calculate various non-LTE parameters for line
!   Called by ADJUST.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE HITCOM_DAT ! HITRAN line data
    USE NTECOM_DAT ! Non-LTE data
    USE PHYCON_DAT, ONLY: ATMB, C2 ! Mathematical and physical constants
!
! SUBROUTINES
    USE IDXQFN_FNC ! Index in QFNCOM of molec,isotope
    USE VAL1DI_GEN ! Interpolate value from 1D array
    USE QTFCT_FNC  ! Calculate total internal partition sum
    USE QTNTE_FNC  ! Calculate total internal partition sums for non-LTE
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN)  :: PRE   ! Path pressure [atm]
    REAL(R4), INTENT(IN)  :: TEM   ! Path kinetic temperature [K]
    REAL(R8), INTENT(IN)  :: GAMMA ! exp - ( hcv/kT )
    REAL(R8), INTENT(OUT) :: ANLTE ! Correction factor for k absorption
    REAL(R8), INTENT(OUT) :: CNLTE ! Correction factor for c absorption
    REAL(R4), INTENT(OUT) :: SQ    ! Ratio TIPS @296K/@path.tem
!
! LOCAL CONSTANTS 
    REAL(R4), PARAMETER :: PRELIM = 1.0E-5 ! Pres [atm] defining high alt limit
    REAL(R4), PARAMETER :: TEMLIM = 275.0  ! Default Max Vib.T for high alt.emis
    REAL(R4), PARAMETER :: FACLIM = 0.01   ! Min VT/KT ratio equiv to VT=0
!
! LOCAL VARIABLES
    LOGICAL     :: BOTDAT ! T = Vib.Temp available for lower level
    LOGICAL     :: TOPDAT ! T = Vib.Temp available for upper level
    INTEGER(I4) :: IQFN   ! Index of Rot.Part.Fn.
    REAL(R4)    :: ENBOT  ! Energy [/cm] of lower level
    REAL(R4)    :: ENTOP  ! Energy [/cm] of upper level
    REAL(R4)    :: LNP    ! Log(p/mb) for interpolating Vib.Temp. profiles
    REAL(R4)    :: QVTEM  ! Non-LTE Vib. partition fn.
    REAL(R4)    :: TEMDEF ! Default vibrational temperature
    REAL(R4)    :: VTBOT  ! Vib.Temp of lower level      
    REAL(R4)    :: VTTOP  ! Vib.Temp of upper level      
    REAL(R8)    :: RLOW   ! Ratio of lower state population to LTE value
    REAL(R8)    :: RUPP   ! Ratio of upper state population to LTE value
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( TEM .GT. TEMLIM .AND. PRE .LT. PRELIM ) THEN
    TEMDEF = TEMLIM
  ELSE
    TEMDEF = TEM
  END IF
  LNP = LOG ( PRE * ATMB )
!
  BOTDAT = HIT%ILV .NE. 0
  IF ( BOTDAT ) THEN
    VTBOT = TEM + VAL1DI ( LNPATM, LNP, VIBATM(:,HIT%ILV), .FALSE., .FALSE. )
    ENBOT  = NTE(HIT%ILV)%ENG
  ELSE
    VTBOT = TEMDEF
    ENBOT = HIT%ELS
  END IF               
  IF ( VTBOT .GT. FACLIM * TEM ) THEN
    RLOW = DEXP ( -C2 * DBLE(ENBOT) *  DBLE(TEM-VTBOT) / DBLE(TEM*VTBOT) )
  ELSE
    RLOW = 0.0D0
  END IF
!
! Note that this routine is only called if HIT%IUV or HIT%ILV .NE. 0, so either 
! TOPDAT or BOTDAT or both will be TRUE.
! If there is no data for the top level, there must be data for the bottom
! level so assume the top level is in LTE with the bottom (RUPP=RLOW)
!
  TOPDAT = HIT%IUV .NE. 0
  IF ( TOPDAT ) THEN
    VTTOP = TEM + VAL1DI ( LNPATM, LNP, VIBATM(:,HIT%IUV), .FALSE., .FALSE. )
    ENTOP  = NTE(HIT%IUV)%ENG
    IF ( VTTOP .GT. FACLIM * TEM ) THEN
      RUPP = DEXP ( -C2 * DBLE(ENTOP) * DBLE(TEM-VTTOP) / DBLE(TEM*VTTOP) )
    ELSE 
      RUPP = 0.0D0
    END IF
  ELSE
    RUPP = RLOW
  END IF
!
  ANLTE = ( RLOW - RUPP * GAMMA ) / ( 1.D0 - GAMMA )
  CNLTE = RUPP
!
! Partition functions (just use LTE functions for the moment)
  IQFN = IDXQFN ( HIT%IDM, HIT%IDI )
  IF ( IQFN .EQ. 0 ) THEN               ! No rot.part.fn defined
    SQ = QTFCT ( HIT%IDM, HIT%IDI, TEM )
  ELSE 
    QVTEM = VAL1DI ( LNPATM, LNP, QFNATM(:,IQFN), .FALSE., .FALSE. )
    SQ = QTNTE ( HIT%IDM, HIT%IDI, TEM, QVTEM )
  END IF
!
END SUBROUTINE NTECLC
END MODULE NTECLC_SUB

