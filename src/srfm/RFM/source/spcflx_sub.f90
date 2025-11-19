MODULE SPCFLX_SUB
CONTAINS
SUBROUTINE SPCFLX 
!
! VERSION
!   04MAY24 AD Checked.
!   05MAR19 AD Simplify by using FLXATM. Checked.
!   01JUN17 AD F90 conversion of rfmflx.for. Checked.
!
! DESCRIPTION
!   Spectral flux calculation
!   Called by RFMSPC if FLX flag enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE QADCOM_DAT ! Gaussian quadrature data
    USE FINCOM_DAT, ONLY: NFIN, WNOFIN ! Finemesh data
    USE SFCCOM_DAT, ONLY: RFLSFC ! T=reflective surface
    USE PHYADJ_DAT, ONLY: TEMSPA ! Cosmic background temperature
!
! SUBROUTINES
    USE FLXATM_SUB ! Atmospheric flux calculation
    USE FLXSFC_SUB ! Surface radiance flux
    USE PLANCK_FNC ! Planck function    
    USE RADMTX_SUB ! Calculate radiance matrix
    USE TRAMTX_SUB ! Calculate transmittance matrix
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IQAD            ! Counter for quadrature points
    REAL(R8)    :: BBFSPA(NFIN)    ! Space radiance
    REAL(R8)    :: RQAD(NFIN,NQAD) ! Radiances for quadrature paths
!       
! EXECUTABLE CODE -------------------------------------------------------------
!
! Initialise space radiance
  BBFSPA = PLANCK ( TEMSPA, WNOFIN )
  DO IQAD = 1, NQAD
    RQAD(:,IQAD) = BBFSPA
  END DO
!
! Downward path - not required if NADFLG and non-reflecting surface
  IF ( RFLSFC .OR. .NOT. NADFLG ) CALL FLXATM ( .TRUE., RQAD )
!
  IF ( ZENFLG ) RETURN   ! only consider downwelling radiances
!
! Surface contribution to reflected radiance
  CALL FLXSFC ( WNOFIN, WQAD, RQAD ) 
!
! Upward path
  CALL FLXATM ( .FALSE., RQAD )
!
  IF ( MTXFLG ) THEN
    IF ( ABSFLG .OR. TRAFLG ) CALL TRAMTX
    IF ( RADFLG .OR. COOFLG ) CALL RADMTX
  END IF
!
END SUBROUTINE SPCFLX
END MODULE SPCFLX_SUB
