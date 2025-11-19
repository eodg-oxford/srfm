MODULE SRCBFX_FNC
CONTAINS
FUNCTION SRCBFX ( ITAN, IATM, IDIR, WNOFIN, OPTLAY, SRCLAY, PSI )
!
! VERSION
!   05NOV23 AD Checked.
!   01MAY17 AD F90 original, adapted from parts of radtra.for. Checked.
!
! DESCRIPTION
!   Modified source function for BFX flag
!   General purpose module.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE TANCOM_DAT ! Tangent path data
    USE ATMCOM_DAT, ONLY: HGTATM ! Atmospheric profile altitudes [km]
    USE PHYADJ_DAT, ONLY: RADCRV ! Local radius of curvature [km]
!
! SUBROUTINES
    USE BBFWGT_FNC ! Interpolation weights for BFX flag
    USE PLANCK_FNC ! Planck Function
    USE TEMVAL_FNC ! Interpolate Temperature from internal profile
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: ITAN ! Index of ray path
    INTEGER(I4), INTENT(IN) :: IATM ! Index of atmospheric layer
    INTEGER(I4), INTENT(IN) :: IDIR ! +1 = upward path, -1 = downward path
    REAL(R8),    INTENT(IN) :: WNOFIN(:) ! Wavenumber [cm-1] array
    REAL(R8),    INTENT(IN) :: OPTLAY(:) ! Layer optical depth array
    REAL(R8),    INTENT(IN) :: SRCLAY(:) ! Absorber-weighted source function
    REAL(R4), OPTIONAL, &
                 INTENT(IN) :: PSI  ! Horiz.angle [deg] of near end of path
!
! FUNCTION TYPE
    REAL(R8) :: SRCBFX(SIZE(WNOFIN)) ! Output array size of WNOFIN
!
! LOCAL VARIABLES
    INTEGER(I4) :: ILEVF  ! Atmospheric profile level# farthest from observer
    INTEGER(I4) :: ILEVN  ! Atmospheric profile level# nearest to observer
    REAL(R4)    :: HGTN   ! Altitude [km] at near level
    REAL(R4)    :: TANHGT ! Projected tangent altitude [km]
    REAL(R4)    :: TEMN   ! Temperature [K] at near profile level
    REAL(R4)    :: XTANF  ! Distance [km] from farthest level to tangent point
    REAL(R4)    :: XTANN  ! Distance [km] from nearest level to tangent point
    REAL(R8)    :: BBFN(SIZE(WNOFIN)) ! Planck fn evaluated at near level
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Set index of profile level nearest the observer
  IF ( IDIR .EQ. 1 ) THEN
    ILEVN = IATM
    ILEVF = IATM + 1
  ELSE
    ILEVN = IATM + 1
    ILEVF = IATM
  END IF
!
  IF ( LIMTAN ) THEN
! Limb viewing: calculate distances from nearest/farthest level to tangent point
! In tangent layer, need to interpolate for tangent point temperature
! as nearest point on upward path
    TANHGT = MIN ( TAN(ITAN)%HGT, TAN(ITAN)%GEO )
    XTANN = SQRT ( 2.0 * ( SNGL(RADCRV) + TANHGT ) & 
                      * MAX ( ( HGTATM(ILEVN) - TANHGT ), 0.0 ) ) 
    XTANF = SQRT ( 2.0 * ( SNGL(RADCRV) + TANHGT ) &
                      * MAX ( ( HGTATM(ILEVF) - TANHGT ), 0.0 ) ) 
    HGTN = MAX ( TAN(ITAN)%HGT, HGTATM(ILEVN) ) 
  ELSE              ! plane-parallel or homogeneous atmosphere
    XTANN  = 0.0
    XTANF  = 0.0
    HGTN = HGTATM(ILEVN) 
  END IF
!
  IF ( PRESENT ( PSI ) ) THEN
    TEMN = TEMVAL ( HGTN, PSI ) 
  ELSE
    TEMN = TEMVAL ( HGTN ) 
  END IF
  BBFN = PLANCK ( TEMN, WNOFIN )
  SRCBFX = BBFN + BBFWGT ( XTANN, XTANF, OPTLAY ) * ( SRCLAY - BBFN ) 
!
END FUNCTION SRCBFX
END MODULE SRCBFX_FNC

