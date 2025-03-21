MODULE FLXSFC_SUB
CONTAINS
SUBROUTINE FLXSFC ( WNO, WQAD, RQAD )
!
! VERSION
!   21SEP24 AD Checked.
!   01MAY17 AD F90 conversion of rfmflx.for. Checked.
!
! DESCRIPTION
!   Surface radiance flux
!   Called by RADMTX, SPCFLX
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE SFCCOM_DAT ! Surface parameters
    USE PHYCON_DAT, ONLY: PI
!
! SUBROUTINES
    USE INTERP_GEN ! Interpolate array
    USE PLANCK_FNC ! Planck function
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(IN)    :: WNO(:)    ! Wavenumber axis [cm-1]
    REAL(R8), INTENT(IN)    :: WQAD(:)   ! Quadrature weights
    REAL(R8), INTENT(INOUT) :: RQAD(:,:) ! Radiance at various angles
!
! LOCAL VARIABLES
    INTEGER(I4) :: IQAD              ! Counter for quadrature points
    INTEGER(I4) :: NQAD              ! No. quadrature points
    REAL(R8)    :: BBFSFC(SIZE(WNO)) ! Surface radiance
    REAL(R8)    :: EMS(SIZE(WNO))    ! Surface emissivity
    REAL(R8)    :: RADBOA(SIZE(WNO)) ! Bottom of atmosphere radiance
!       
! EXECUTABLE CODE -------------------------------------------------------------
!
  NQAD = SIZE ( RQAD, 2 ) 
  BBFSFC = PLANCK ( TEMSFC, WNO )
!
  IF ( RFLSFC ) THEN
    EMS = INTERP ( WNOSFC, WNO, EMSSFC )
! Assume Lambertian surface, so integrate radiances prior to reflection
    RADBOA = MATMUL ( RQAD, WQAD ) / PI  ! /PI to convert irrad to rad
    DO IQAD = 1, NQAD
      RQAD(:,IQAD) = RADBOA * (1.0D0-EMS) + EMS * BBFSFC
    END DO
  ELSE
    DO IQAD = 1, NQAD
      RQAD(:,IQAD) = BBFSFC
    END DO
  END IF
!
END SUBROUTINE FLXSFC
END MODULE FLXSFC_SUB
