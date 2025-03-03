MODULE BBFWGT_FNC
CONTAINS
FUNCTION BBFWGT ( XTAN0, XTAN1, CHI )
!
! VERSION
!   03NOV23 AD Checked.
!   01MAY17 AD F90 conversion of subroutine bbfopt.for. Checked.
!
! DESCRIPTION
!   Interpolation weights for BFX flag
!   Called by SRCBFX
!   If XTAN0 = 0: interpreted as Zen/Nad viewing (=no layer curvature).
!   These weights (W) are to be applied in the equation:
!       DRad = Tau * ( B0 + W * ( BCG - B0 ) ) * ( 1 - EXP(-Opt) )
!   Where Drad = radiance contribution from layer, 
!          Tau = transmission from observer to near end of layer
!           B0 = Planck Fn at near end of layer
!          BCG = Absorber-weighted (CG) Planck Function
!          Opt = Optical depth of layer (CHI)
!            W = Weight (BBFWGT) returned by routine, 0 .LE. W .LE. 1
!     
!   XTAN0 is distance from edge of layer nearest to the observer to tangent pt
!   XTAN1 is distance from edge of layer remote from the observer to tan.pt.
!   If both XTAN0 = XTAN1 = 0, interpreted as zero curvature.
!   XTAN1,XTAN0 are distances in arbitrary units,
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: XTAN0  ! Distance from near edge of layer to tan. pt
    REAL(R4), INTENT(IN) :: XTAN1  ! Distance from  far edge of layer to tan. pt
    REAL(R8), INTENT(IN) :: CHI(:) ! Spectrum of layer optical depths
!
! FUNCTION TYPE
    REAL(R8) :: BBFWGT ( SIZE(CHI) ) ! Function array same size as CHI
!
! LOCAL CONSTANTS
    REAL(R8), PARAMETER :: CHIMIN = 0.01D0 ! Minimum optical depth for full calc
    REAL(R4), PARAMETER :: XLIMIT = 0.01   ! Ignore curvature if layer thickness
                                           ! < XLIMIT * dist to tan pt
! LOCAL VARIABLES
    REAL(R8) :: CFACT  !  Curv. param (2:infty =tangent layer:no curv).
    REAL(R8) :: FACTOR !  1 / ( 2 * CFACT - 1 )
    REAL(R8) :: RDTAU(SIZE(CHI)) ! Recip. of layer transm 1 / ( 1 - exp(-CHI) )
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! Both XTAN0 and XTAN1 should be either +ve or zero
  IF ( XTAN0 .LT. 0.0 ) STOP 'F-BBFOPT: Logical error#0'
  IF ( XTAN1 .LT. 0.0 ) STOP 'F-BBFOPT: Logical error#1'
!
! Curvature is ignored if both XTAN1=XTAN0=0, or if the difference between them
! is small compared to distance to tangent point (ie thin layer and/or remote 
! from tangent point).
  IF ( ABS(XTAN0-XTAN1) .LE. XLIMIT * MIN(XTAN0,XTAN1) ) THEN ! No Curve
    WHERE ( CHI .LT. CHIMIN )            ! Use limiting value
      BBFWGT = 1.0D0
    ELSEWHERE                            ! Explicit calculation
      RDTAU = 1.0D0 / ( 1.0D0 - EXP(-CHI) )
      BBFWGT = 2.0D0 + 2.0D0 / CHI - 2.0D0 * RDTAU
    END WHERE
!
  ELSE                                           ! curvature
! For near-side tangent layer XTAN1=0 so CFACT=1, FACTOR=1
! For far-side tangent layer, XTAN0=0 so CFACT=0, FACTOR=-1
! For other near-side layers, CFACT > 1 so +1 < FACTOR < +1/XLIMIT (approx)
! For other far-side layers,  CFACT < 0 so -1 > FACTOR > -1/XLIMIT
    CFACT = DBLE ( XTAN0 / ( XTAN0 - XTAN1 ) )
    FACTOR = 1.0D0 / ( 2.0D0 * CFACT - 1.0D0 )
    WHERE ( CHI .LT. CHIMIN ) 
      BBFWGT = 1.0D0                      ! Use limiting value
    ELSEWHERE                                       ! Explicit calculation
      RDTAU = 1.0D0 / ( 1.0D0 - EXP(-CHI) )
      BBFWGT = 6.0D0 / ( 3.0D0 + FACTOR ) &
               * ( 1.0D0 + 1.0D0/CHI - RDTAU + ( FACTOR / CHI**2 ) & 
                   * ( 2.0D0 * CHI * RDTAU - CHI - 2.0D0 ) )
    ENDWHERE
  END IF
!
END FUNCTION BBFWGT
END MODULE BBFWGT_FNC
