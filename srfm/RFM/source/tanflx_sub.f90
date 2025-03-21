MODULE TANFLX_SUB
CONTAINS
SUBROUTINE TANFLX 
!
! VERSION
!   30AUG24 AD Checked.
!   24JUN19 AD Remove COOWGT, called by DRVCHK instead. Checked.
!   05MAR19 AD Use INIQAD
!   01JUN17 AD F90 original. Checked.
!
! DESCRIPTION
!   Initialise spectral flux calculation
!   Called by DRVTAN if FLX flag enabled.
!   Assumes that atmospheric profile now fixed in size.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE TANCOM_DAT ! Tangent path data
    USE ATMCOM_DAT, ONLY: NATM, ITNATM ! Indices of Flux ouput levels (or 0)
!
! SUBROUTINES
    USE INIQAD_SUB ! Initialise Gaussian quadrature for flux calculations
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IATM ! Index of atmospheric level
    INTEGER(I4) :: ITAN ! Counter for output levels
!       
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL INIQAD
!
  ALLOCATE ( ITNATM(NATM) )
  ITNATM = 0
  DO ITAN = 1, NTAN 
    IATM = TAN(ITAN)%IAT
    ITNATM(IATM) = ITAN
  END DO
!
END SUBROUTINE TANFLX
END MODULE TANFLX_SUB
