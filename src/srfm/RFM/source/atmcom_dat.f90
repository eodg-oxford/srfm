MODULE ATMCOM_DAT
!
! VERSION
!   01OCT24 AD Checked.
!   18APR22 AD Add IAIVMR. Checked.
!   28JUN18 AD Bug#9 Initialise NATM=0 to indicate 'not yet defined'
!   21JUN17 AD Add IATSFC, FIXPRE, SETHGT, LEVATM
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Atmospheric profile data.
!   Values are defined at profile levels 1:NATM. ATMCOM_RESET deallocates all
!   associated storage so the model can be reinitialised between successive
!   Python-driven runs.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
  PUBLIC :: ATMCOM_RESET
!
! GLOBAL CONSTANTS
     INTEGER(I4), PARAMETER :: LENATM = 20 ! Max length of .atm profile label
!
! GLOBAL VARIABLES
    LOGICAL     :: FIXPRE = .FALSE. ! T=fixed lnp grid, F=fixed hgt grid
    LOGICAL     :: SETHGT = .FALSE. ! T = Altitude profile set
    LOGICAL     :: SETPRE = .FALSE. ! T = Pressure profile set
    LOGICAL     :: SETTEM = .FALSE. ! T = Temperature profile set
    INTEGER(I4) :: IATSFC = 1       ! Index of surface level
    INTEGER(I4) :: IAIVMR           ! Index of air profile
    INTEGER(I4) :: IAXVMR           ! Index of aerosol extinction profile
    INTEGER(I4) :: NATM = 0         ! No. atmospheric layers for profiles
    INTEGER(I4) :: NVIB             ! No. diff. Vibrational temp profiles
    INTEGER(I4) :: NVMR             ! No. different VMR profiles
    REAL(R4)    :: HGTSFC           ! Height [km] of bottom of atmosphere
    REAL(R4)    :: HGTTOA           ! Height [km] of top of atmosphere
    REAL(R4)    :: PRESFC           ! Pressure [mb] of bottom of atmosphere
!
    LOGICAL,     ALLOCATABLE :: LINVMR(:)   ! T = linear interpolation with alt
    LOGICAL,     ALLOCATABLE :: NTEVMR(:)   ! T = non-LTE species
    LOGICAL,     ALLOCATABLE :: SETVMR(:)   ! T = VMR profiles set
    INTEGER(I4), ALLOCATABLE :: ICOATM(:,:) ! Indices for cooling rates
    INTEGER(I4), ALLOCATABLE :: ITNATM(:)   ! Indices of Flux ouput levels, or 0
    INTEGER(I4), ALLOCATABLE :: NCOATM(:)   ! No.cooling rate indices
    REAL(R4),    ALLOCATABLE :: DNSATM(:)   ! Number Density [/cm^3]
    REAL(R4),    ALLOCATABLE :: DSHATM(:)   ! Density Scale Height [km]
    REAL(R4),    ALLOCATABLE :: EXTATM(:)   ! Extinction profile [km-1]
    REAL(R4), TARGET, &
                 ALLOCATABLE :: HGTATM(:)   ! Altitude profile [km]
    REAL(R4),    ALLOCATABLE :: LNDATM(:)   ! Ln ( DNSATM )
    REAL(R4), TARGET, &
                 ALLOCATABLE :: LNPATM(:)   ! Ln ( PREATM )
    REAL(R4),    ALLOCATABLE :: PREATM(:)   ! Pressure profile [mb]
    REAL(R4),    ALLOCATABLE :: QFNATM(:,:) ! [NATM,NVIB] Vib.Part.Fn profiles
    REAL(R4),    ALLOCATABLE :: RFRATM(:)   ! Refractivity profile
    REAL(R4),    ALLOCATABLE :: TEMATM(:)   ! Temperature profile [K]
    REAL(R4),    ALLOCATABLE :: VIBATM(:,:) ! [NATM,NVIB] Vib Temp profiles [K]
    REAL(R4),    ALLOCATABLE :: VMRATM(:,:) ! [NATM,NVMR] VMR profiles [ppmv]
    REAL(R4),    POINTER     :: LEVATM(:)   ! Profile levels (hgt or lnp)   
    REAL(R8),    ALLOCATABLE :: WCOATM(:,:) ! Weights for cooling rates
    CHARACTER(LENATM), ALLOCATABLE :: NAMVMR(:) ! Names of VMR profiles
!
CONTAINS

  SUBROUTINE ATMCOM_RESET()
    IF ( ALLOCATED ( LINVMR ) )    DEALLOCATE ( LINVMR )
    IF ( ALLOCATED ( NTEVMR ) )    DEALLOCATE ( NTEVMR )
    IF ( ALLOCATED ( SETVMR ) )    DEALLOCATE ( SETVMR )
    IF ( ALLOCATED ( ICOATM ) )    DEALLOCATE ( ICOATM )
    IF ( ALLOCATED ( ITNATM ) )    DEALLOCATE ( ITNATM )
    IF ( ALLOCATED ( NCOATM ) )    DEALLOCATE ( NCOATM )
    IF ( ALLOCATED ( DNSATM ) )    DEALLOCATE ( DNSATM )
    IF ( ALLOCATED ( DSHATM ) )    DEALLOCATE ( DSHATM )
    IF ( ALLOCATED ( EXTATM ) )    DEALLOCATE ( EXTATM )
    IF ( ALLOCATED ( HGTATM ) )    DEALLOCATE ( HGTATM )
    IF ( ALLOCATED ( LNDATM ) )    DEALLOCATE ( LNDATM )
    IF ( ALLOCATED ( LNPATM ) )    DEALLOCATE ( LNPATM )
    IF ( ALLOCATED ( PREATM ) )    DEALLOCATE ( PREATM )
    IF ( ALLOCATED ( QFNATM ) )    DEALLOCATE ( QFNATM )
    IF ( ALLOCATED ( RFRATM ) )    DEALLOCATE ( RFRATM )
    IF ( ALLOCATED ( TEMATM ) )    DEALLOCATE ( TEMATM )
    IF ( ALLOCATED ( VIBATM ) )    DEALLOCATE ( VIBATM )
    IF ( ALLOCATED ( VMRATM ) )    DEALLOCATE ( VMRATM )
    IF ( ASSOCIATED ( LEVATM ) )   NULLIFY ( LEVATM )
    IF ( ALLOCATED ( WCOATM ) )    DEALLOCATE ( WCOATM )
    IF ( ALLOCATED ( NAMVMR ) )    DEALLOCATE ( NAMVMR )

    FIXPRE = .FALSE.
    SETHGT = .FALSE.
    SETPRE = .FALSE.
    SETTEM = .FALSE.
    IATSFC = 1
    IAIVMR = 0
    IAXVMR = 0
    NATM   = 0
    NVIB   = 0
    NVMR   = 0
    HGTSFC = 0.0_R4
    HGTTOA = 0.0_R4
    PRESFC = 0.0_R4
  END SUBROUTINE ATMCOM_RESET

END MODULE ATMCOM_DAT
