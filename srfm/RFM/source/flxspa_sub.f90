MODULE FLXSPA_SUB
CONTAINS
SUBROUTINE FLXSPA ( WNO, RQAD ) 
!
! VERSION
!   13SEP24 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Space radiance flux
!   Called by RADMTX.
!   Initialises the array RQAD containing pencil beams at various angles used
!   for hemispherical integration with Planck radiance appropriate to the
!   temperature of cold space.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE PHYADJ_DAT, ONLY: TEMSPA ! Cosmic background temperature
!
! SUBROUTINES
    USE PLANCK_FNC ! Planck function
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(IN)  :: WNO(:)    ! Wavenumber axis [cm-1]
    REAL(R8), INTENT(OUT) :: RQAD(:,:) ! Radiance at various angles
!
! LOCAL VARIABLES
    INTEGER(I4) :: IQAD              ! Counter for quadrature points
    REAL(R8)    :: BBFSPA(SIZE(WNO)) ! Space radiance
!       
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( TEMSPA .GT. TINY(0.0) ) THEN
    BBFSPA = PLANCK ( TEMSPA, WNO )
    DO IQAD = 1, SIZE ( RQAD, 2 ) 
      RQAD(:,IQAD) = BBFSPA
    END DO
  ELSE
    RQAD = 0.0D0
  END IF
!
END SUBROUTINE FLXSPA
END MODULE FLXSPA_SUB
