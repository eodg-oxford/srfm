MODULE ATMAUX_SUB
CONTAINS
SUBROUTINE ATMAUX 
!
! VERSION
!   11MAY24 AD Checked.
!   12MAR19 AD Allow for WNORFR. Checked.
!   30DEC17 AD F90 version. Checked.
!
! DESCRIPTION
!   Set up auxiliary profiles of atmospheric parameters
!   Called by DRVCHK and PTBATM.
!   Assumes TEM, PRE, VMR, EXT profiles are specified.
!   Derive profiles of Refractivity, Density, Log(Density), Log(Pressure),
!   Log(VMR) and Density Scale Height.
!
! VARIABLE KINDS
    USE KIND_DAT 
! 
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data    
    USE GRACOM_DAT ! Atmospheric 2-D field
    USE FLGCOM_DAT, ONLY: GEOFLG          ! T = geom. (ie non-refrac) ray paths
    USE PHYADJ_DAT, ONLY: GRAVTY, WGTAIR, WNORFR  ! g, WgtAir, Wno.refrac.
    USE PHYCON_DAT, ONLY: AVOG, RGAS      ! Avog.No & Molar gas constant
!
! SUBROUTINES
    USE REFRAC_FNC ! Calculate refractivity
!
  IMPLICIT NONE
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: ARGMIN = TINY ( 1.0 )  ! Smallest Real 
    REAL(R4), PARAMETER :: LOGMIN = LOG ( ARGMIN )  
!
! LOCAL VARIABLES
    INTEGER(I4) :: IATM   ! Counter for profile levels
    INTEGER(I4) :: IPSI   ! Counter for horizontal locations
    REAL(R4), ALLOCATABLE :: LNDPSI(:) ! Log(dens) profiles at horiz locations
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( ALLOCATED ( DNSATM ) ) THEN
    IF ( SIZE(DNSATM) .NE. NATM ) DEALLOCATE ( DNSATM, LNDATM, RFRATM, DSHATM ) 
  END IF
!
  IF ( .NOT. ALLOCATED ( DNSATM ) ) &
    ALLOCATE ( DNSATM(NATM), LNDATM(NATM), RFRATM(NATM), DSHATM(NATM) ) 
! 
! Lnp profile may be separately allocated from other profiles
  IF ( ALLOCATED ( LNPATM ) ) THEN
    IF ( SIZE(LNPATM) .NE. NATM ) DEALLOCATE ( LNPATM ) 
  END IF
  IF ( .NOT. ALLOCATED ( LNPATM ) ) ALLOCATE ( LNPATM(NATM) ) 
!
  LNPATM = LOG ( MAX ( PREATM, ARGMIN ) )
!
! Convert p [mb] to [Pa]: *100 
! Convert dens [m-3] to [cm-3]: *1E-6
  DNSATM = PREATM * AVOG * 1.0E-4 / RGAS / TEMATM 
  LNDATM = LOG ( MAX ( DNSATM, ARGMIN ) )
!
  IF ( GEOFLG ) THEN
    RFRATM = 0.0
  ELSE 
    RFRATM = REFRAC ( PREATM, TEMATM, WNORFR )
  END IF
!
! Convert ext [km-1] to [cm-1]: *1E-5
! Convert [ppv] tp [ppmv]: *1E6
  IF ( IAXVMR .GT. 0 ) VMRATM(:,IAXVMR) = EXTATM * 10.0 / DNSATM
!
! If layers are too narrow cannot calculate scale height accurately
! from difference, so use approximation RT/gM (R/gM=8.3/(9.8*29)=0.029km/K)
! Also prevent any negative density gradients
  DO IATM = 2, NATM
    IF ( LNDATM(IATM-1) - LNDATM(IATM) .LT. 0.01 ) THEN
      DSHATM(IATM-1) = RGAS / GRAVTY / WGTAIR * 1.0E-3 * TEMATM(IATM-1) 
!
! Use refractivity difference to calculate scale height if possible 
! NB similar calculation in VALGRA for 'DSH' option
    ELSE IF (  RFRATM(IATM) .GT. ARGMIN*1.0E6 ) THEN
      DSHATM(IATM-1) = ( HGTATM(IATM) - HGTATM(IATM-1) ) / &
                   LOG ( RFRATM(IATM-1)/RFRATM(IATM) )
! If refractivity too small to calculate accurately, use density scale height
    ELSE IF ( LNDATM(IATM) .GT. LOGMIN ) THEN
      DSHATM(IATM-1) = ( HGTATM(IATM) - HGTATM(IATM-1) ) / &
                       ( LNDATM(IATM-1) - LNDATM(IATM) )
    ELSE
      DSHATM(IATM-1) = 1.0E20
      WRITE (*,*) 'F-ATMAUX: Setting DSHATM large'
      stop
    END IF
  END DO
!
! Convert 2D profile quantities
  IF ( NPSI .GT. 0 ) THEN
    DNSGRA = PREGRA * AVOG * 1.0E-4 / RGAS / TEMGRA
    IF ( GEOFLG ) THEN
      RFRGRA = 0.0
    ELSE
      DO IPSI = 1, NPSI    ! NPSI = 0 unless GRA flag enabled
        RFRGRA(:,IPSI) = REFRAC ( PREGRA(:,IPSI), TEMGRA(:,IPSI), WNORFR ) 
      END DO
    END IF
    ALLOCATE ( LNDPSI(NATM) ) 
    DO IPSI = 1, NPSI    ! NPSI = 0 unless GRA flag enabled
      RFRGRA(:,IPSI) = REFRAC ( PREGRA(:,IPSI), TEMGRA(:,IPSI), WNORFR ) 
      LNDPSI = LOG ( MAX ( DNSGRA(:,IPSI), ARGMIN ) )
      DO IATM = 2, NATM
        IF ( LNDPSI(IATM-1) - LNDPSI(IATM) .LT. 0.01 ) THEN
          DSHGRA(IATM-1,IPSI) = 0.029 * TEMGRA(IATM-1,IPSI) 
        ELSE IF ( RFRGRA(IATM,IPSI) .GT. ARGMIN*1.0E6 ) THEN
          DSHGRA(IATM-1,IPSI) = ( HGTATM(IATM) - HGTATM(IATM-1) ) / &
                   LOG ( RFRGRA(IATM-1,IPSI)/RFRGRA(IATM,IPSI) )
        ELSE IF ( LNDPSI(IATM) .GT. LOGMIN ) THEN
          DSHGRA(IATM-1,IPSI) = ( HGTATM(IATM) - HGTATM(IATM-1) ) / &
                       ( LNDPSI(IATM-1) - LNDPSI(IATM) )
        ELSE 
          DSHGRA(IATM-1,IPSI) = 1.0E20
          STOP 'F-ATMAUX: Setting DSHGRA large'
        END IF
      END DO
    END DO
    IF ( IAXVMR .GT. 0 ) VMRGRA(:,:,IAXVMR) = EXTGRA * 10.0 / DNSGRA
  END IF
!
END SUBROUTINE ATMAUX
END MODULE ATMAUX_SUB
