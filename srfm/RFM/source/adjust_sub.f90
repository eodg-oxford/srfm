MODULE ADJUST_SUB
CONTAINS
SUBROUTINE ADJUST ( TEM, PRE, PPA, AMT, ANTE, CNTE )
!
! VERSION
!   04AUG24 AD Checked.
!   11AUG23 AD Allow for different line parameter sets   
!   31MAY23 AD Allow for HITRAN database entries
!   29JAN20 AD New HITCOM variable names. Checked.
!              Allow for line mixing data in RFM basic hitran format
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Adjust line parameters for path conditions
!   Called by SPCWID, SPCFIN
!   The line strength and width at 296K (TEMREF) and 1 atm (PREF) read from 
!   the HITRAN line data file are adjusted for atmospheric path conditions. 
!   The doppler halfwidth is also calculated. These parameters are then used 
!   in the lineshape formulation. 
!
!   Adjusts data in HITCOM and loads into ADJCOM
! 
!   Line mixing data and temperature dependence interpolation due to 
!   L.L.Strow (private communication with Dave Edwards)
!
!   Line strengths for different isotopes of a gas on the HITRAN data base 
!   weighted by atmospheric abundance of the particular isotope. Absolute 
!   strengths may be obtained by dividing by this abundance i.e. 
!   STREN/ABUN(IDGAS,ISO) where IDGAS is the gas ID and ISO the isotope ID. 
!   This will be important when performing calculations for other planetary 
!   atmospheres.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE ADJCOM_DAT ! Path-adjusted line data
    USE HITCOM_DAT ! HITRAN line data
    USE SETCOM_DAT ! HITRAN line parameter sets
    USE IDXCON_DAT, ONLY: IDXCO2  ! RFM/HITRAN index for CO2
    USE PHYCON_DAT, ONLY: C2, PREREF, RGAS, TEMREF, VLIGHT ! Physical constants
    USE FLGCOM_DAT, ONLY: MIXFLG ! T = use line-mixing
!
! SUBROUTINES
    USE NTECLC_SUB ! Calculate various non-LTE parameters for line
    USE QTFCT_FNC  ! Calculate total internal partition sum
    USE YMIX_FNC   ! Calculate line-mixing y-coefficient
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN)  :: TEM  ! Path temperature [K]
    REAL(R4), INTENT(IN)  :: PRE  ! Path pressure [atm]
    REAL(R4), INTENT(IN)  :: PPA  ! Path partial pressure [atm]
    REAL(R4), INTENT(IN)  :: AMT  ! Path amount [kmol/cm2]
    REAL(R4), INTENT(OUT) :: ANTE ! Non-lte factor for k abs
    REAL(R4), INTENT(OUT) :: CNTE ! Non-lte factor for c abs
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: R2 = 2.0 * LOG(2.0) * RGAS ! 2ln2 k N = 11526.3
!
! LOCAL VARIABLES
    LOGICAL, POINTER :: SET_LMA ! T=value for air line-mixing 
    LOGICAL, POINTER :: SET_LMS ! T=value for self line-mixing
    LOGICAL, POINTER :: SET_PSS ! T=value for self pressure-shift
    LOGICAL, POINTER :: SET_TCS ! T=value for self-broad. Temp coeff.
    INTEGER(I4) :: ISET   ! Index for line parameter sets  
    REAL(R4)    :: SQ     ! Ratio of tps@296K/tps@path_temp
    REAL(R4)    :: TFACT  ! TEMREF/TEM - temperature scale factor
    REAL(R8)    :: ANLTE  ! Non-LTE Correction factor for k absorption
    REAL(R8)    :: CNLTE  ! Non-LTE Correction factor for c absorption
    REAL(R8)    :: GAMMA  ! exp ( -hcv/kT )
    REAL(R8)    :: GAMREF ! exp ( -hcv/kT_ref )
    REAL(R8)    :: SB     ! exp( -hcE_l/kT_path ) / exp( -hcE_l/kT_ref )
    REAL(R8)    :: SE     ! Ratio of stimulated emission @path/@ref
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Save path parameters - could be required by CHISHP
  TEMADJ = TEM
  PREADJ = PRE
  PPAADJ = PPA
  TFACT = TEMREF/TEM
!
  ISET = HIT%IST
  SET_LMA => SET(IPAR_LMA,ISET)
  SET_LMS => SET(IPAR_LMS,ISET)
  SET_PSS => SET(IPAR_PSS,ISET)
  SET_TCS => SET(IPAR_TCS,ISET)
!
! Pressure shift (often 0)  
  IF ( SET_PSS ) THEN
    WNOADJ = HIT%WNO + DBLE ( (PRE-PPA) * HIT%PSA + PPA * HIT%PSS )
  ELSE
    WNOADJ = HIT%WNO + DBLE ( PRE * HIT%PSA )
  END IF
!
! Convert for line width in cm-1 at 296K and 1atm.
! Assume air, self half-width (HWA,HWS) & Tcoeff (TCA) are always present
  IF ( SET_TCS ) THEN  ! Only use TCS if HWS is present
    WIDADJ = ( HIT%HWA * ( PRE - PPA ) * TFACT**HIT%TCA + &
               HIT%HWS *         PPA   * TFACT**HIT%TCS     ) / PREREF
  ELSE
    WIDADJ = ( HIT%HWA * ( PRE - PPA ) + &
               HIT%HWS *         PPA       ) * TFACT**HIT%TCA / PREREF
  END IF
!
! Calculate Doppler half-width at half-max HWHM in cm-1. 
  DOPADJ = SNGL ( HIT%WNO / VLIGHT ) * SQRT ( R2 * TEM / HIT%WGT )
!
! Calculate the line mixing y coefficient
  IF ( MIXFLG ) THEN
    IF ( SET_LMS ) THEN
      YMXADJ = ( HIT%LMA * ( PRE - PPA ) * TFACT**HIT%TCA + &
                 HIT%LMS *         PPA   * TFACT**HIT%TCS   ) / PREREF
    ELSE IF ( SET_LMA ) THEN
      YMXADJ = HIT%LMA * PRE * TFACT**HIT%TCA / PREREF
    ELSE IF ( HIT%IDM .EQ. IDXCO2 ) THEN ! Internally stored CO2 line-mixing 
      YMXADJ = YMIX ( TEM, PRE, PPA )
    ELSE
      YMXADJ = 0.0
    END IF
  ENDIF
!
! Convert for line strength in cm-1.(mol.cm-2)-1 at 296K.
!
! Boltzman factor for lower state energy
  SB = DEXP ( DBLE(HIT%ELS) * C2 * DBLE(TEM-TEMREF)/DBLE(TEM*TEMREF) )
!
! Stimulated emission 
  GAMMA  = DEXP ( -C2 * HIT%WNO / DBLE ( TEM ) )
  GAMREF = DEXP ( -C2 * HIT%WNO / DBLE ( TEMREF ) )
  SE = ( 1.D0 - GAMMA ) / ( 1.D0 - GAMREF )
!
! Nonlte calculation of absorption coefficient modifiers
  SQ = 1.0                       
  IF ( HIT%IUV .NE. 0 .OR. HIT%ILV .NE. 0 ) THEN
!    IF ( PTH(IPTH)%IVJ .GT. 0 ) CALL PTBVIB ( PTH(IPTH)%IVJ, .TRUE. ) 
    CALL NTECLC ( PRE, TEM, GAMMA, ANLTE, CNLTE, SQ )
!    IF ( PTH(IPTH)%IVJ .GT. 0 ) CALL PTBVIB ( PTH(IPTH)%IVJ, .FALSE. ) 
  ELSE
    ANLTE = 1.0D0
    CNLTE = 1.0D0       
    SQ = QTFCT ( HIT%IDM, HIT%IDI, TEM )
  ENDIF
!
! SB can be larger than allowed for SNGL (eg for mid-IR NO lines in HITRAN2012)
! so combine all factors together before converting to SNGL in the hope that it
! will be small enough to fit into STRADJ
  STRADJ = SNGL ( AMT * HIT%STR * SB * SE * SQ )   ! Bug#112
  ANTE = SNGL ( ANLTE )
  CNTE = SNGL ( CNLTE )
!
END SUBROUTINE ADJUST
END MODULE ADJUST_SUB
