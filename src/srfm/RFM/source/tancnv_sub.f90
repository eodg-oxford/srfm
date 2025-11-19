MODULE TANCNV_SUB
CONTAINS
SUBROUTINE TANCNV ( MODE, ITAN )
!
! VERSION
!   23SEP24 AD Bug#48 Ignore surface intercept if upwards-viewing
!   05AUG19 AD Change ALTOBS to HGTOBS. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Convert tangent point specifications
!   Called by CHKLIM.
!   This uses the relationships: n r sin(z) = n t = constant along a 
!   refracted path, where 
!     n = refractive index, n = 1 + R where R is refractivity.
!     r = radial distance from the centre of curvature 
!     z = zenith angle of ray (=angle of incidence/refraction) =90 at Tan.Pt.
!     t = radial distance of projected tangent point.
!   The function F is defined as 
!             F = n(z) z - n_obs t_obs      where F = 0 if z = z_tan
!   This is solved iteratively: 
!             z' = z - F/(dF/dz)  where dF/dz =  1 + N(z) + z (dN/dz)
!   Since refractivity is proportional to air density, it is assumed that 
!   the refractivity gradient is given by
!             dN/dz = -N/H      where H(z) is the density scale height
!
!   Note that single precision is adequate for 1 metre resolution
!   when adding Earth's radius to altitude terms, but not for adding 
!   refractivity to unity to get refractive index.
!
!   The ELETAN variable is only set if OBSFLG =.T. (Observer alt.specified)
!   In the case of rays intersecting the surface, HGTTAN is set = surface
!   altitude and SZNTAN is set to the sine(Zenith Angle) at that point.
!   For +ve ELETAN values, HGTTAN is set to observer altitude, and GEOTAN set
!   to projected tangent altitude in opposite direction (-ELETAN).
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE FLGCOM_DAT ! Option flags
    USE OBSCOM_DAT ! Observer location data
    USE TANCOM_DAT ! Tangent path data
    USE PHYADJ_DAT, ONLY: RADCRV ! Local radius of curvature [km]
    USE PHYCON_DAT, ONLY: DG2RAD ! [rad]/[deg] conversion factor (pi/180).
!
! SUBROUTINES
    USE RFRVAL_FNC ! Interpolate Refractivity profile
    USE TANITR_SUB ! Iterative determination of tangent point
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: MODE ! #1=ELETAN,#2=GEOTAN,#3=HGTTAN defined 
    INTEGER(I4), INTENT(IN) :: ITAN ! Index of tangent path
!
! LOCAL VARIABLES
    REAL(R8) :: COSELE ! Cosine(Elevation Angle)
    REAL(R8) :: RADGEO ! Radial distance [km] of geometric tangent point
    REAL(R8) :: RADOBS ! Radial distance [km] of observer
    REAL(R8) :: RADTAN ! Radial distance [km] of refracted tangent point
    REAL(R8) :: RFROBS ! Refractivity at obs.location (or zero if space)
    REAL(R8) :: RFRTAN ! (Iterated) refractivity at tangent point
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  SELECT CASE ( MODE )
!
  CASE ( 1 ) ! Mode#1: TAN(ITAN)%ELE specified (and OBSFLG is TRUE).
    IF ( .NOT. OBSFLG ) STOP 'F-TANCNV: Logical Error#1'
    COSELE = COS ( DBLE ( TAN(ITAN)%ELE * DG2RAD ) )
    RADOBS = DBLE ( HGTOBS ) + RADCRV
    RADGEO = RADOBS * COSELE            ! Either up or down viewing
    TAN(ITAN)%GEO = SNGL ( RADGEO - RADCRV ) 
    IF ( TAN(ITAN)%ELE .GT. 0.0 ) THEN         ! Upward viewing
      TAN(ITAN)%HGT = HGTOBS
      TAN(ITAN)%SZN = COSELE
    ELSE                                ! Downward viewing
      CALL TANITR ( HGTOBS, TAN(ITAN)%GEO, TAN(ITAN)%HGT, TAN(ITAN)%SZN )
    END IF
!
  CASE ( 2 ) ! Mode#2: TAN(ITAN)%GEO specified 
    IF ( OBSFLG ) THEN
      RADGEO = DBLE ( TAN(ITAN)%GEO ) + RADCRV
      RADOBS = DBLE ( HGTOBS ) + RADCRV
      COSELE = RADGEO / RADOBS
      TAN(ITAN)%ELE = - SNGL ( ACOS ( COSELE ) ) / DG2RAD
      CALL TANITR ( HGTOBS, TAN(ITAN)%GEO, TAN(ITAN)%HGT, TAN(ITAN)%SZN )
    ELSE
      CALL TANITR ( HGTTOA, TAN(ITAN)%GEO, TAN(ITAN)%HGT, TAN(ITAN)%SZN )
    END IF
!          
  CASE ( 3 ) ! Mode#3: TAN(ITAN)%HGT specified 
    RADTAN = DBLE ( TAN(ITAN)%HGT ) + RADCRV
    IF ( OBSFLG ) THEN
      RADOBS = DBLE ( HGTOBS ) + RADCRV
      IF ( GEOFLG ) THEN
        RADGEO = RADTAN
      ELSE
        RFRTAN = RFRVAL ( TAN(ITAN)%HGT )
        RFROBS = DBLE ( RFRATM(IATOBS) ) 
        RADGEO = ( 1.0D0 + RFRTAN ) * RADTAN / ( 1.0D0 + RFROBS )
      END IF
      TAN(ITAN)%GEO = SNGL ( RADGEO - RADCRV )
      COSELE = RADGEO / RADOBS
      TAN(ITAN)%ELE = - SNGL ( ACOS ( COSELE ) ) / DG2RAD
      TAN(ITAN)%SZN = 1.0D0 
    ELSE
      IF ( GEOFLG ) THEN
        RADGEO = RADTAN
      ELSE
        RFRTAN = RFRVAL ( TAN(ITAN)%HGT )
        RADGEO = ( 1.0D0 + RFRTAN ) * RADTAN 
      END IF
      TAN(ITAN)%GEO = SNGL ( RADGEO - RADCRV )
      TAN(ITAN)%SZN = 1.0D0 
    END IF
!
  CASE DEFAULT 
    STOP 'F-TANCNV: Logical Error#2'
  END SELECT 
!
  TAN(ITAN)%SFC = SFCFLG .AND. TAN(ITAN)%HGT .EQ. HGTSFC
! Special case of observer at the surface looking upward  Bug#48
  IF ( OBSFLG .AND. TAN(ITAN)%SFC ) &
    TAN(ITAN)%SFC = TAN(ITAN)%ELE .LT. 0.0 
!
END SUBROUTINE TANCNV
END MODULE TANCNV_SUB

