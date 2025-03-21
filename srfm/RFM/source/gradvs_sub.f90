MODULE GRADVS_SUB
CONTAINS
SUBROUTINE GRADVS ( P, R, T, DPDS, DRDS, DTDS )
!
! VERSION
!   13JAN24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Derivatives for 2-D ray-tracing
!   Called by RAYGRA.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE PHYADJ_DAT, ONLY: RADCRV ! Local radius of curvature [km]
    USE FLGCOM_DAT, ONLY: GEOFLG ! T = geometric (ie non-refrac) ray paths
    USE PHYCON_DAT, ONLY: DG2RAD ! [rad]/[deg] conversion factor
!
! SUBROUTINES
    USE DRPVAL_FNC ! Evaluate d(Refractivity)/d(Psi) [/rad]
    USE DSHVAL_FNC ! Interpolate Density Scale Height
    USE RFRVAL_FNC ! Interpolate Refractivity profile
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(IN)  :: P    ! D.P. horizontal angle [rad]
    REAL(R8), INTENT(IN)  :: R    ! Radius coordinate of path [km]
    REAL(R8), INTENT(IN)  :: T    ! D.P. zenith angle [rad]
    REAL(R8), INTENT(OUT) :: DPDS ! Rate of change of P with S
    REAL(R8), INTENT(OUT) :: DRDS ! Rate of change of R with S
    REAL(R8), INTENT(OUT) :: DTDS ! Rate of change of T with S
!
! LOCAL VARIABLES
    REAL(R4) :: HGT    ! Height [km] of path segment
    REAL(R4) :: PSI    ! Horizontal angle [deg] of path
    REAL(R8) :: DNDP   ! Rate of change of refractive index with P
    REAL(R8) :: DNDR   ! Rate of change of refractive index with R
    REAL(R8) :: REFRAC ! Refractivity (=Ref.index - 1)
    REAL(R8) :: RFRIDX ! Refractive index
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  HGT = SNGL ( R - RADCRV )
  PSI = SNGL ( P ) / DG2RAD
  IF ( GEOFLG ) THEN
    DNDR = 0.0D0
    DNDP = 0.0D0
    RFRIDX = 1.0D0
  ELSE
    REFRAC = RFRVAL ( HGT, PSI ) 
    RFRIDX = 1.0D0 + REFRAC
    DNDR = - REFRAC / DBLE ( DSHVAL ( HGT, PSI ) )
    DNDP = DBLE ( DRPVAL ( HGT, PSI ) )
  END IF
!
  DRDS = COS ( T )
  DPDS = SIN ( T ) / R
  DTDS = COS ( T ) * DNDP / ( RFRIDX * R ) &
       - SIN ( T ) * ( 1.0D0/R + DNDR/RFRIDX ) 
!
END SUBROUTINE GRADVS
END MODULE GRADVS_SUB

