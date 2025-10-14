MODULE GRACNV_SUB
CONTAINS
SUBROUTINE GRACNV ( MODE, ITAN ) 
!
! VERSION
!   26SEP24 AD Checked.
!   05AUG19 AD Change ALTOBS to HGTOBS. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Convert tangent point specifications for 2D atmosphere
!   Called by CHKLIM.
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
    USE OBSCOM_DAT ! Observer location data
    USE TANCOM_DAT ! Tangent path data
    USE ATMCOM_DAT, ONLY: HGTSFC, & ! Surface altitude [km]
                          HGTTOA    ! Top of atmosphere altitude [km]
    USE FLGCOM_DAT, ONLY: OBSFLG, & ! T = observer located within atmosphere
                          SFCFLG    ! T = T = Allow for surface 
    USE PHYADJ_DAT, ONLY: RADCRV ! Local radius of curvature [km]
    USE PHYCON_DAT, ONLY: DG2RAD ! [rad]/[deg] conversion factor
!
! SUBROUTINES
    USE RAYGRA_SUB ! Ray-tracing in a 2-D atmosphere (z,psi)
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: MODE ! #1=ELETAN,#2=GEOTAN,#3=HGTTAN defined 
    INTEGER(I4), INTENT(IN) :: ITAN ! Index of tangent path
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: MAXITR = 10 ! Max. no. of iterations for tan point
    REAL(R4),    PARAMETER :: ALTTOL = 0.0003 ! Tolerance [km] for refrac.tan.ht
    REAL(R4),    PARAMETER :: DPATH  = 1.0    ! Path increment for ray-tracing
!                                             ! 1km seems to give accuracy limit
! LOCAL VARIABLES
    LOGICAL     :: CNVRGD ! T = converged
    INTEGER(I4) :: ITER   ! Iteration counter
    REAL(R4)    :: DS     ! Max Path length [km] for each ray-trace segment
    REAL(R4)    :: DZENDH ! Rate of change of ZENOBS with HGTTST
    REAL(R4)    :: HGTPRV ! Previous value of HGTTST
    REAL(R4)    :: HGTTST ! Test value iterated towards HGTTAN
    REAL(R4)    :: PSITOA ! Horizontal angle [deg] of top of atmosphere (+ve)
    REAL(R4)    :: ZENOBS ! Zenith angle [deg] at observer
    REAL(R4)    :: ZENPRV ! Previous estimate of ZENOBS
    REAL(R4)    :: ZENTAN ! Zenith angle [deg] at actual tan.ht (nominally 270)
    REAL(R4)    :: ZENTOA ! Zenith angle [deg] at top of atmosphere
    REAL(R8)    :: COSELE ! Cosine(Elevation Angle)
    REAL(R8)    :: COSTOA ! Cosine(Elevation Angle at top of atmosphere)
    REAL(R8)    :: RADGEO ! Radial distance [km] of geometric tangent point
    REAL(R8)    :: RADOBS ! Radial distance [km] of observer
    REAL(R8)    :: RADTAN ! Radial distance [km] of tangent point
    REAL(R8)    :: RADTOA ! Radial distance [km] of top of atmosphere
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DS = DPATH       ! 1.0 = 1km  max path length for ray tracing`
!
  SELECT CASE ( MODE ) 
!
  CASE ( 1 )                 ! Mode#1: TAN(ITAN)%ELE specified (and OBS = TRUE)
    IF ( .NOT. OBSFLG ) STOP 'F-GRACNV: Logical Error#1'
    COSELE = COS ( DBLE ( TAN(ITAN)%ELE * DG2RAD ) )
    RADOBS = DBLE ( HGTOBS ) + RADCRV 
    RADGEO = RADOBS * COSELE
    TAN(ITAN)%GEO = SNGL ( RADGEO - RADCRV ) ! For both +/- TAN(ITAN)%ELE
    IF ( TAN(ITAN)%ELE .GT. 0.0 ) THEN         ! Upward viewing
      TAN(ITAN)%HGT = HGTOBS
      TAN(ITAN)%SZN = COSELE
    ELSE
      ZENOBS = 270.0 + TAN(ITAN)%ELE
      ZENTAN = 270.0 
      CALL RAYGRA ( HGTOBS, PSIOBS, ZENOBS, 3, DS, &
                    TAN(ITAN)%HGT, TAN(ITAN)%PSI, ZENTAN )
      IF ( ZENTAN .EQ. 270.0 ) THEN ! Do this to avoid precision problems
        TAN(ITAN)%SZN = -1.0D0      
      ELSE
        TAN(ITAN)%SZN = SIN ( DBLE ( ZENTAN * DG2RAD ) )
      END IF
    END IF
!
  CASE ( 2 )                 ! Mode#2: TAN(ITAN)%GEO specified 
    IF ( OBSFLG ) THEN                      ! Observer location specified
      RADGEO = DBLE ( TAN(ITAN)%GEO ) + RADCRV 
      RADOBS = DBLE ( HGTOBS ) + RADCRV 
      COSELE = RADGEO / RADOBS
      TAN(ITAN)%ELE = -SNGL ( ACOS ( COSELE ) ) / DG2RAD  ! -ve for down viewing
      ZENOBS = 270.0 + TAN(ITAN)%ELE
      TAN(ITAN)%SZN = 270.0
      CALL RAYGRA ( HGTOBS, PSIOBS, ZENOBS, 3, DS, &
                    TAN(ITAN)%HGT, TAN(ITAN)%PSI, ZENTAN )
      IF ( ZENTAN .EQ. 270.0 ) THEN ! Do this to avoid precision problems
        TAN(ITAN)%SZN = -1.0D0      
      ELSE
        TAN(ITAN)%SZN = SIN ( DBLE ( TAN(ITAN)%SZN * DG2RAD ) )
      END IF
    ELSE
      RADTOA = DBLE ( HGTTOA ) + RADCRV
      RADGEO = DBLE ( TAN(ITAN)%GEO ) + RADCRV
      PSITOA = SNGL ( ACOS ( RADGEO / RADTOA ) ) / DG2RAD
      ZENTOA = 270.0 - PSITOA                  ! Assume PSITOA is positive
      ZENTAN = 270.0
      CALL RAYGRA ( HGTTOA, PSITOA, ZENTOA, 3, DS, &
                    TAN(ITAN)%HGT, TAN(ITAN)%PSI, ZENTAN )
      IF ( ZENTAN .EQ. 270.0 ) THEN ! Do this to avoid precision problems
        TAN(ITAN)%SZN = -1.0D0      
      ELSE
        TAN(ITAN)%SZN = SIN ( DBLE ( ZENTAN * DG2RAD ) )
      END IF
    END IF
!
  CASE ( 3 )                 ! Mode#3: TAN(ITAN)%HGT specified
    IF ( OBSFLG ) THEN
      RADOBS = DBLE ( HGTOBS ) + RADCRV
      RADTAN = DBLE ( TAN(ITAN)%HGT ) + RADCRV
      COSELE = RADTAN / RADOBS         ! Geometric tan.hgt for first guess
      ZENOBS = 270.0 - SNGL ( ACOS ( COSELE ) ) / DG2RAD   !-ve elevation
      ZENPRV = 270.0
      HGTPRV = HGTOBS
      CNVRGD = .FALSE.
      DO ITER = 1, MAXITR
        ZENTAN = 270.0                  
        DS = DPATH
        CALL RAYGRA ( HGTOBS, PSIOBS, ZENOBS, 3, DS, &
                      HGTTST, TAN(ITAN)%PSI, ZENTAN )
        CNVRGD = ( ABS ( HGTTST - TAN(ITAN)%HGT ) .LT. ALTTOL ) .OR. &
                 ( ABS ( HGTTST - HGTPRV ) .LT. ALTTOL ) 
        IF ( CNVRGD )  EXIT 
        DZENDH = ( ZENOBS - ZENPRV ) / ( HGTTST - HGTPRV )
        ZENPRV = ZENOBS
        HGTPRV = HGTTST
        ZENOBS = ZENPRV + ( TAN(ITAN)%HGT - HGTTST ) * DZENDH 
      END DO
      IF ( .NOT. CNVRGD ) THEN
        WRITE ( *, * ) 'W-GRACNV: Max.iterations exceeded, residual=', &
           TAN(ITAN)%HGT - HGTTST
        ZENOBS = 0.5 * ( ZENOBS + ZENPRV ) ! Assuming oscillating about soln
      END IF
      IF ( ZENOBS .GT. 270.0 ) STOP 'F-GRACNV: Logical Error#2'
      TAN(ITAN)%ELE = ZENOBS - 270.0                          ! Should be -ve
      COSELE = COS ( DBLE ( TAN(ITAN)%ELE * DG2RAD ) )
      RADGEO = RADOBS * COSELE 
      TAN(ITAN)%GEO = SNGL ( RADGEO - RADCRV )
      TAN(ITAN)%SZN = 1.0D0
    ELSE
      TAN(ITAN)%PSI = 0.0
      ZENTAN = 90.0                   ! 90 since heading *towards* observer
      CALL RAYGRA ( TAN(ITAN)%HGT, TAN(ITAN)%PSI, ZENTAN, 1, DS, &
                    HGTTOA, PSITOA, ZENTOA )
      RADTOA = DBLE ( HGTTOA ) + RADCRV
      COSTOA = SIN ( DBLE ( ZENTOA * DG2RAD ) )
      TAN(ITAN)%ELE = - SNGL ( ACOS ( COSTOA ) ) / DG2RAD
      RADGEO = RADTOA * COSTOA
      TAN(ITAN)%GEO = SNGL ( RADGEO - RADCRV )
      TAN(ITAN)%SZN = 1.0D0
    END IF
  CASE DEFAULT
    STOP 'F-GRACNV: Logical Error#3'      ! Mode .NE. 1, 2 or 3
  END SELECT
!
  TAN(ITAN)%SFC = SFCFLG .AND. TAN(ITAN)%HGT .EQ. HGTSFC
!        
END SUBROUTINE GRACNV
END MODULE GRACNV_SUB
