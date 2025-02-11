MODULE RFRVAL_FNC
CONTAINS
REAL(R8) FUNCTION RFRVAL ( HGT, PSI )
!
! VERSION
!   01JAN24 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Interpolate Refractivity profile
!   General purpose function
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE GRACOM_DAT ! Atmospheric 2-D field
!
! SUBROUTINES
    USE VAL1DI_GEN ! Interpolate value from 1D array
    USE VAL2DI_FNC ! Interpolate value from 2D array
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4),           INTENT(IN) :: HGT ! Altitude [km]
    REAL(R4), OPTIONAL, INTENT(IN) :: PSI ! Horiz.angle [deg]
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( HGT .GE. HGTTOA ) THEN   ! above top of atmosphere
    RFRVAL = 0.0D0
  ELSE IF ( PRESENT ( PSI ) ) THEN 
    RFRVAL = DBLE ( VAL2DI ( HGTATM, HGT, PSIGRA, PSI, RFRGRA, .TRUE. ) )
  ELSE
    RFRVAL = DBLE ( VAL1DI ( HGTATM, HGT, RFRATM, .TRUE. ) )
  END IF
!
END FUNCTION RFRVAL
END MODULE RFRVAL_FNC

