MODULE RADMTX_SUB
CONTAINS
SUBROUTINE RADMTX 
!
! VERSION
!   22JUL26 AD Checked.
!   01AUG25 AD Use ITNOFF to index TAN paths for outputs 
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
    USE QADCOM_DAT ! Gaussian quadrature data
    USE TANCOM_DAT ! Tangent path data
    USE FINCOM_DAT, ONLY: NFIN, WNOFIN ! Finemesh grid
!
! SUBROUTINES
    USE FLXATM_SUB ! Atmospheric flux calculation
    USE FLXSFC_SUB ! Surface radiance flux
    USE FLXSPA_SUB ! Space radiance flux
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IPTB   ! Index of perturbed atm layer
    INTEGER(I4) :: ITAN   ! Counter for output levels
    INTEGER(I4) :: ITNOFF ! Offset in TANCOM for outputs
    REAL(R8)    :: RQAD(NFIN,NQAD) ! Radiances for quadrature paths
!       
! EXECUTABLE CODE -------------------------------------------------------------
!
  ITNOFF = 0
  DO ITAN = 1, NTAN                 ! Perturbed source fn at NTAN levels
    ITNOFF = ITNOFF + NTAN
    IPTB = TAN(ITAN)%IAT
!
! Initialise with subtracted unperturbed paths
    RADFUL(IFUL1:IFUL2,ITNOFF+1:ITNOFF+NTAN) = -RADFUL(IFUL1:IFUL2,1:NTAN) 
    IF ( COOFLG ) &
      COOFUL(IFUL1:IFUL2,ITNOFF+1:ITNOFF+NTAN) = -COOFUL(IFUL1:IFUL2,1:NTAN)
!
    CALL FLXSPA ( WNOFIN, RQAD )    ! Initialise with space radiance
! Downward path
    CALL FLXATM ( .TRUE., RQAD, IPTB, ITNOFF ) 
    IF ( ZENFLG ) CYCLE   ! only consider downwelling radiances
!
! Incorporate surface contribution to bottom-of-atmosphere reflected radiances
    CALL FLXSFC ( WNOFIN, WQAD, RQAD ) 
!
! Upward path
    CALL FLXATM ( .FALSE., RQAD, IPTB, ITNOFF )
!
  END DO
!
END SUBROUTINE RADMTX
END MODULE RADMTX_SUB
