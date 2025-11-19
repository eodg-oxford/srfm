MODULE SPCRAD_SUB
CONTAINS
SUBROUTINE SPCRAD 
!
! VERSION
!   03MAR20 AD Reset OPT=0 at surface for nadir-viewing. Checked.
!   05MAR19 AD Add RADSKY for Diffuse surface reflection.
!              Change arguments for RADLAY
!   01JUL17 AD F90 version. Checked.
!
! DESCRIPTION
!   Radiative transfer calculation
!   Called by RFMSPC
!   Calculation starts at remote side of atmosphere and moves towards observer
!   (NB reverse direction to RFM v4)
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FINCOM_DAT ! Finemesh data
    USE FLGCOM_DAT ! Option flags
    USE FULCOM_DAT ! Full grid data
    USE TANCOM_DAT ! Tangent path data
    USE ATMCOM_DAT, ONLY: NATM   ! No. atmos.profile levels (1+No.layers)
    USE OBSCOM_DAT, ONLY: IATOBS ! Observer level (or 0 if unspecified)
    USE PHYADJ_DAT, ONLY: TEMSPA ! Cosmic background temperature
!
! SUBROUTINES
    USE PLANCK_FNC ! Planck Function
    USE RADLAY_SUB ! Radiative transfer calculation through a layer
    USE RADLEV_SUB ! Save radiances at intermediate output levels
    USE RADSFC_SUB ! Surface contribution to radiative transfer calculation
    USE RADSKY_FNC ! Diffuse sky radiance at base of atmosphere
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IATM         ! Index of atmospheric level
    INTEGER(I4) :: ILOW         ! Lower atm.level for calculation
    INTEGER(I4) :: ISKY         ! Index of TAN path used for diffuse sky rad
    INTEGER(I4) :: ITAN         ! Counter for tangent paths
    INTEGER(I4) :: IUPP         ! Upper atm.level for calculation
    REAL(R8)    :: BBFSPA(NFIN) ! Blackbody space emitted radiance
    REAL(R8)    :: OPTLAY(NFIN) ! Optical thickness of atmos layer
    REAL(R8), POINTER :: OPT(:) ! Optical depth
    REAL(R8), POINTER :: RAD(:) ! Radiance
    REAL(R8), POINTER :: TRA(:) ! Path transmittance
!       
! EXECUTABLE CODE -------------------------------------------------------------
!
  BBFSPA = PLANCK ( TEMSPA, WNOFIN )
!
! Ensure any diffuse sky radiance paths calculated first
  DO ITAN = 1, MTAN
    IF ( TAN(ITAN)%SKY ) RADFUL(IFUL1:IFUL2,ITAN) = RADSKY ( ITAN, BBFSPA )
  END DO
!
  DO ITAN = 1, MTAN
    IF ( .NOT. TAN(ITAN)%CLC ) CYCLE   ! No radiance calculation required
    IF ( TAN(ITAN)%SKY ) CYCLE         ! Already calculated
!
    OPT => OPTFUL(IFUL1:IFUL2,ITAN) ; OPT = 0.0D0
    TRA => TRAFUL(IFUL1:IFUL2,ITAN) ; TRA = 1.0D0
    RAD => RADFUL(IFUL1:IFUL2,ITAN) ; RAD = BBFSPA
!
! If diffuse surface reflection then set radiance to precalc. sky radiance
    IF ( TAN(ITAN)%ISK .GT. 0 ) THEN 
      ISKY = TAN(ITAN)%ISK
      RAD = RADFUL(IFUL1:IFUL2,ISKY) 
!
! Downward path calculation from space into atmosphere lowest point
    ELSE
      IUPP = MAX ( NATM-1, 1 )   ! Max to allow for NATM=1 for HOM flag
      ILOW = TAN(ITAN)%IAT      
!
      IF ( LEVFLG ) CALL RADLEV ( ITAN, NATM, -1 ) 
!
      DO IATM = IUPP, ILOW, -1   ! Loop over layers
        CALL RADLAY ( ITAN, IATM, -1, RAD, OPTLAY )  
        OPT = OPT + OPTLAY
        TRA = EXP ( -OPT )
        IF ( LEVFLG ) CALL RADLEV ( ITAN, IATM, -1 ) 
      END DO
    END IF
!
! Skip upward calculation if not required
    IF ( ZENFLG ) CYCLE                      ! zenith view
    IF ( HOMFLG .AND. .NOT. SFCFLG ) CYCLE   ! single pass for HOM flag
!
! Surface contribution
    IF ( TAN(ITAN)%SFC ) THEN
      CALL RADSFC ( OPT, RAD, TAN(ITAN)%JDX ) 
      IF ( NADFLG ) OPT = 0.0D0              ! Ensure OPT,TRA for atm.only
      TRA = EXP ( -OPT ) 
    END IF
!
! Upward path from lowest point to observer
    ILOW = TAN(ITAN)%IAT
    IF ( OBSFLG ) THEN
      IUPP = IATOBS - 1
    ELSE 
      IUPP = MAX ( NATM-1, 1 )     ! max to allow for NATM=1 for HOM flag
    END IF
!
    DO IATM = ILOW, IUPP
      IF ( LEVFLG ) CALL RADLEV ( ITAN, IATM, 1 ) 
      CALL RADLAY ( ITAN, IATM, 1, RAD, OPTLAY )  
      OPT = OPT + OPTLAY
      TRA = EXP ( -OPT )
    END DO      
    IF ( LEVFLG ) CALL RADLEV ( ITAN, IATM, 1 ) 
!
  END DO
!
  NULLIFY ( OPT, RAD, TRA ) 
!
END SUBROUTINE SPCRAD
END MODULE SPCRAD_SUB
