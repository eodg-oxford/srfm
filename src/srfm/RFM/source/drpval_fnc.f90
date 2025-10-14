MODULE DRPVAL_FNC
CONTAINS
REAL(R4) PURE FUNCTION DRPVAL ( HGT, PSI )
!
! VERSION
!   14JAN24 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Evaluate d(Refractivity)/d(Psi) [/rad]
!   Called by GRADVS
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GRACOM_DAT ! Atmospheric 2-D field
    USE ATMCOM_DAT, ONLY: HGTATM ! Profile altitude levels
    USE PHYCON_DAT, ONLY: DG2RAD ! pi/180 [rad]/[deg] conversion factor
!
! SUBROUTINES
    USE VAL2DI_FNC ! Interpolate value from 2D array
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: HGT ! Altitude [km]
    REAL(R4), INTENT(IN) :: PSI ! LOS angle [deg]
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( NPSI .EQ. 1 ) THEN ! No horizontal gradient
    DRPVAL = 0.0
  ELSE                    ! Since horiz.grad is +/-, don't do log interp.
    DRPVAL = VAL2DI ( HGTATM, HGT, PSIGRA, PSI, RFRGRA, .FALSE. )
  END IF
!
END FUNCTION DRPVAL
END MODULE DRPVAL_FNC
