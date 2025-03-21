MODULE TEMVAL_FNC
CONTAINS
REAL(R4) FUNCTION TEMVAL ( HGT, PSI )
!
! VERSION
!   05JAN24 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Interpolate Temperature from internal profile
!   Called by SRCBFX if BFX flag enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GRACOM_DAT, ONLY: PSIGRA, TEMGRA ! Atmospheric 2-D field
    USE ATMCOM_DAT, ONLY: HGTATM, TEMATM ! Atmospheric profile data
!
! SUBROUTINES
    USE VAL1DI_GEN ! Interpolate value from 1D array
    USE VAL2DI_FNC ! Interpolate value from 2D array
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4),           INTENT(IN) :: HGT ! Altitude [km]
    REAL(R4), OPTIONAL, INTENT(IN) :: PSI ! LOS angle [deg]
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( PRESENT ( PSI ) ) THEN
    TEMVAL = VAL2DI ( HGTATM, HGT, PSIGRA, PSI, TEMGRA )
  ELSE
    TEMVAL = VAL1DI ( HGTATM, HGT, TEMATM )
  END IF
!
END FUNCTION TEMVAL
END MODULE TEMVAL_FNC

