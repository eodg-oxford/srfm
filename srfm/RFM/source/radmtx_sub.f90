MODULE RADMTX_SUB
CONTAINS
SUBROUTINE RADMTX 
!
! VERSION
!   01MAY24 AD Checked.
!   05MAR19 AD Use FLXATM. Checked.
!   01JUN17 AD F90 conversion of rfmflx.for. Checked.
!
! DESCRIPTION
!   Calculate radiance matrix
!   Called by SPCFLX if MTX and (RAD or COO) flags enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE FULCOM_DAT ! Full grid data
    USE LEVCOM_DAT ! Intermediate output levels
    USE QADCOM_DAT ! Gaussian quadrature data
    USE FINCOM_DAT, ONLY: NFIN, WNOFIN ! Finemesh grid
    USE TANCOM_DAT, ONLY: NTAN         ! No. nominal output levels
!
! SUBROUTINES
    USE FLXATM_SUB ! Atmospheric flux calculation
    USE FLXSFC_SUB ! Surface radiance flux
    USE FLXSPA_SUB ! Space radiance flux
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    LOGICAL     :: PTBSRC ! T=perturbed source function
    INTEGER(I4) :: ILEV   ! Counter for output levels
    INTEGER(I4) :: ITAN   ! Counter for output levels
    INTEGER(I4) :: JTAN   ! Index of flux matrix elements
    REAL(R8)    :: RQAD(NFIN,NQAD) ! Radiances for quadrature paths
!       
! EXECUTABLE CODE -------------------------------------------------------------
!
  PTBSRC = .FALSE.
  DO ILEV = 1, NLEV         ! Perturbed source fn levels
    CALL FLXSPA ( WNOFIN, RQAD )    ! Initialise with space radiance
! Downward path
    CALL FLXATM ( .TRUE., RQAD, ILEV ) 
    IF ( ZENFLG ) CYCLE   ! only consider downwelling radiances
!
! Incorporate surface contribution to bottom-of-atmosphere reflected radiances
    CALL FLXSFC ( WNOFIN, WQAD, RQAD ) 
!
! Upward path
    CALL FLXATM ( .FALSE., RQAD, ILEV )    
!
  END DO
!
! Subtract unperturbed path values (stored 1:NTAN)
  DO ITAN = 1, NTAN
    DO ILEV = 1, NLEV
      JTAN = ITNLEV(ITAN,ILEV)
      RADFUL(IFUL1:IFUL2,JTAN) = RADFUL(IFUL1:IFUL2,JTAN) - &
                                 RADFUL(IFUL1:IFUL2,ITAN) 
      IF ( COOFLG ) COOFUL(IFUL1:IFUL2,JTAN) = COOFUL(IFUL1:IFUL2,JTAN) - &
                                               COOFUL(IFUL1:IFUL2,ITAN) 
    END DO
  END DO
!
END SUBROUTINE RADMTX
END MODULE RADMTX_SUB
