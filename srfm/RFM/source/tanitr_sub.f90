MODULE TANITR_SUB
CONTAINS
SUBROUTINE TANITR ( HGTINI, HGTGEO, HGTTAN, SZNTAN )
!
! VERSION
!   27DEC23 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Iterative determination of tangent point
!   Called by TANCNV
!   This uses the relationships: n r sin(z) = n t = constant along a 
!   refracted path, where 
!     n = refractive index, n = 1 + R where R is refractivity.
!     r = radial distance from the centre of curvature 
!     z = zenith angle of ray (=angle of incidence/refraction) =90 at Tan.Pt.
!     t = radial distance of projected tangent point.
!   The function F is defined as 
!             F = n(z) z - n_ini t_ini      where F = 0 if z = z_tan
!   This is solved iteratively: 
!             z' = z - F/(dF/dz)  where dF/dz =  1 + N(z) + z (dN/dz)
!   Since refractivity is proportional to air density, it is assumed that 
!   the refractivity gradient is given by
!             dN/dz = -N/H      where H(z) is the density scale height
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE ATMCOM_DAT, ONLY: HGTSFC ! Altitude of surface/base of atmosphere
    USE FLGCOM_DAT, ONLY: GEOFLG ! T = geometric (ie non-refrac) ray paths
    USE PHYADJ_DAT, ONLY: RADCRV ! Local radius of curvature [km]
!
! SUBROUTINES
    USE DSHVAL_FNC ! Interpolate Density Scale Height
    USE RFRVAL_FNC ! Interpolate Refractivity profile
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN)  :: HGTINI ! Initial altitude [km]
    REAL(R4), INTENT(IN)  :: HGTGEO ! Initial geometric tangent height [km]
    REAL(R4), INTENT(OUT) :: HGTTAN ! Tan.ht or surface height [km]
    REAL(R8), INTENT(OUT) :: SZNTAN ! Zenith angle at tan.pt or surface [deg]
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: MAXITR = 10 ! Max. no. of iterations for tan pt
    REAL(R4),    PARAMETER :: ALTTOL = 0.0003 ! Tol. [km] for refr.tan.ht.calcs
!
! LOCAL VARIABLES
    LOGICAL     :: CNVRGD ! T=converged
    INTEGER(I4) :: ITR    ! Iteration counter
    REAL(R4)    :: HGTPRV ! Value of ALTTAN on previous iteration
    REAL(R4)    :: HGTTST ! (Iterated) refracted tangent point altitude
    REAL(R8)    :: DFBYDZ ! dF/dz
    REAL(R8)    :: DRBYDZ ! d(Refractivity)/dz, assuming Rf.has dens.sc.ht.
    REAL(R8)    :: F      ! Iterated function for tangent point (see above)
    REAL(R8)    :: FCONST ! Constant term in calculation of F
    REAL(R8)    :: RADGEO ! Radial distance [km] of geometric tangent point
    REAL(R8)    :: RADSFC ! Radial distance [km] of surface
    REAL(R8)    :: RADTST ! (Iterated) Radial distance [km] of refr.t.pt.
    REAL(R8)    :: RFRINI ! Refractivity at obs.location (or zero if space)
    REAL(R8)    :: RFRSFC ! Refractivity at base of atmosphere
    REAL(R8)    :: RFRTST ! (Iterated) refractivity at tangent point
    REAL(R8)    :: SZNSFC ! Sin(Zen.angle) at surface
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  RADSFC = DBLE ( HGTSFC ) + RADCRV
  RADGEO = DBLE ( HGTGEO ) + RADCRV
!
  IF ( GEOFLG ) THEN
    RFRSFC = 0.0D0
    RFRINI = 0.0D0
  ELSE
    RFRSFC = RFRVAL ( HGTSFC )
    RFRINI = RFRVAL ( HGTINI )
  END IF
!
  SZNSFC = ( 1.0D0 + RFRINI ) * RADGEO / ( 1.0D0 + RFRSFC ) / RADSFC 
!
  IF ( SZNSFC .LE. 1.0D0 ) THEN      ! Ray intersects surface
    SZNTAN = SZNSFC
    HGTTAN = HGTSFC
  ELSE IF ( GEOFLG ) THEN
    SZNTAN = 1.0D0
    HGTTAN = HGTGEO
  ELSE
    SZNTAN = 1.0D0
    HGTTST = HGTGEO                        ! First guess
    FCONST  = ( 1.D0 + RFRINI ) * RADGEO 
    DO ITR = 1, MAXITR
      RFRTST = RFRVAL ( HGTTST ) 
      RADTST = DBLE ( HGTTST ) + RADCRV
      F = ( 1.0D0 + RFRTST ) * RADTST - FCONST
      DRBYDZ = - RFRTST / DSHVAL ( HGTTST )
      DFBYDZ = 1.D0 + RFRTST + ( RADTST * DRBYDZ )
      HGTPRV = HGTTST
      HGTTST = MAX ( HGTTST - SNGL ( F / DFBYDZ ), HGTSFC )
      CNVRGD = ABS ( HGTPRV - HGTTST ) .LE. ALTTOL 
      IF ( CNVRGD ) THEN
        HGTTAN = HGTTST
        EXIT
      END IF
    END DO
! If routine fails to converge, assume this is due to small oscillations 
! about the true value 
    IF ( .NOT. CNVRGD ) THEN 
      WRITE ( *, * ) 'W-TANITR: Failed to converge, residual=', HGTTST - HGTPRV
      HGTTAN = 0.5 * ( HGTTST + HGTPRV )
    END IF
  END IF
!
END SUBROUTINE TANITR
END MODULE TANITR_SUB
