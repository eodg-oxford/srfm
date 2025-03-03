MODULE RADSKY_FNC
CONTAINS
FUNCTION RADSKY ( ISKTAN, RADTOA ) 
!
! VERSION
!   02MAY24 AD Checked.
!   05MAR19 AD Original. Checked.
!
! DESCRIPTION
!   Diffuse sky radiance at base of atmosphere
!   Called by SPCRAD
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FINCOM_DAT ! Finemesh data
    USE PTHCOM_DAT ! Path segment data
    USE QADCOM_DAT ! Gaussian quadrature data
    USE ATMCOM_DAT, ONLY: IATSFC, NATM, NVMR ! Atmospheric profile data
    USE FLGCOM_DAT, ONLY: BFXFLG ! T = Source Fn varies with layer optical dept!h
    USE PHYCON_DAT, ONLY: PI
!
! SUBROUTINES
    USE IDXPTH_FNC ! Index in PTHCOM of tan/atm/gas/dir
    USE PLANCK_FNC ! Planck Function
    USE SRCBFX_FNC ! Modified source function for BFX flag
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: ISKTAN       ! Index of vertical tangent path
    REAL(R8),    INTENT(IN) :: RADTOA(NFIN) ! Radiance at top of atmosphere
!
! FUNCTION TYPE
    REAL(R8) :: RADSKY ( NFIN ) 
!
! LOCAL VARIABLES
    INTEGER(I4) :: IATM            ! Index of atmospheric level
    INTEGER(I4) :: ICLC            ! Index of calculated path
    INTEGER(I4) :: IPTH            ! Path counter
    INTEGER(I4) :: IQAD            ! Counter over sec theta quadrature angles
    INTEGER(I4) :: IVMR            ! Absorber counter
    REAL(R8)    :: OPTGAS(NFIN)    ! Optical thickness of one absorber
    REAL(R8)    :: OPTLAY(NFIN)    ! Optical thickness of atmos layer
    REAL(R8)    :: RQAD(NFIN,NQAD) ! Radiances for quadrature paths
    REAL(R8)    :: SRCLAY(NFIN)    ! Layer source function
    REAL(R8)    :: TRALAY(NFIN)    ! Layer transmittance
!       
! EXECUTABLE CODE -------------------------------------------------------------
!
! Initialise to space radiance
  DO IQAD = 1, NQAD
    RQAD(:,IQAD) = RADTOA
  END DO
!
! This effectively reproduces FLXLAY but scales AMT to get vertical path amount
  DO IATM = NATM-1, IATSFC, -1            ! Loop over layers downwards
    OPTLAY = 0.0D0
    SRCLAY = 0.0D0
    DO IVMR = 1, NVMR
      IPTH = IDXPTH ( ISKTAN, IATM, IVMR )
      ICLC = PTH(IPTH)%ICL
      OPTGAS = DBLE ( ABSFIN(:,ICLC) * PTH(IPTH)%AMT ) 
      SRCLAY = SRCLAY + OPTGAS * PLANCK ( PTH(IPTH)%TEM, WNOFIN )
      OPTLAY = OPTLAY + OPTGAS
    END DO
!
    SRCLAY = SRCLAY / SIGN ( MAX ( ABS(OPTLAY), TINY(1.0D0) ), OPTLAY )
!
    DO IQAD = 1, NQAD
      TRALAY = EXP ( -OPTLAY / XQAD(IQAD) )
      IF ( BFXFLG ) &
        SRCLAY = SRCBFX ( ISKTAN, IATM, -1, WNOFIN, OPTLAY/XQAD(IQAD), SRCLAY )
      RQAD(:,IQAD) = RQAD(:,IQAD) * TRALAY + SRCLAY * ( 1.0D0 - TRALAY )
    END DO
  END DO
!
  RADSKY = MATMUL ( RQAD, WQAD ) / PI  ! /PI to convert irrad to rad
!
END FUNCTION RADSKY
END MODULE RADSKY_FNC


