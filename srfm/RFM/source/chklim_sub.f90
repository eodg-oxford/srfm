MODULE CHKLIM_SUB
CONTAINS
SUBROUTINE CHKLIM ( ITAN, ELETAN, GEOTAN, FAIL, ERRMSG )
!
! VERSION
!   23SEP34 AD Bug#48 Ignore surface intercept if upwards-viewing
!   24MAY24 AD Checked.
!   05AUG19 AD Change ALTOBS to HGTOBS   
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Check Limb-viewing tangent paths
!   Called by CHKTAN, CHKFOV (if FOV flag)
!   Also fills TAN%ELE, TAN%GEO and TAN%HGT
!   Requires inputs from both *TAN and *OBS sections.
!
! VARIABLE KINDS
    USE KIND_DAT
! 
! GLOBAL DATA
    USE TANCOM_DAT ! Tangent path data
    USE FLGCOM_DAT, ONLY: GRAFLG, SFCFLG  ! Option flags
    USE ATMCOM_DAT, ONLY: HGTATM, HGTTOA, HGTSFC ! Atmospheric profile data
    USE OBSCOM_DAT, ONLY: HGTOBS ! Observer altitude [km]
    USE PHYADJ_DAT, ONLY: RADCRV ! Local radius of curvature [km]
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
    USE GRACNV_SUB ! Convert tangent point specifications for 2D atmosphere
    USE IBRAKT_GEN ! Lower index of array interpolation
    USE TANCNV_SUB ! Convert tangent point specifications
!
  IMPLICIT NONE 
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: ITAN   ! Index of tan path to be checked
    LOGICAL,       INTENT(IN)  :: ELETAN ! T= specified as elev.angle [deg]
    LOGICAL,       INTENT(IN)  :: GEOTAN ! T= specified as geom.tan.hgt [km]
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)  :: IATM = 1 ! Guessed value for atm.layer
    REAL(R4)     :: TANTST   ! User-specified value to be checked
    CHARACTER(9) :: TANSTR   ! String containing TANTST for error messages
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  TANTST = TAN(ITAN)%USR
  TANSTR = C9REAL ( TANTST )
!
  FAIL = .TRUE.
  IF ( GEOTAN ) THEN ! TAN(ITAN)%USR contains geometric tangent heights
    TAN(ITAN)%GEO = TANTST
    IF ( TANTST .GT. HGTTOA ) THEN
      ERRMSG = 'F-CHKLIM: Tangent Height=' // TANSTR // &
               '[km] > Top of Atmosphere=' // C9REAL ( HGTTOA ) // '[km]'
    ELSE IF ( SFCFLG .AND. TANTST .LT. -RADCRV ) THEN
      ERRMSG = 'F-CHKLIM: Tangent Height=' // TANSTR // &
               '[km] < Centre of Earth=' // C9REAL ( SNGL ( -RADCRV ) )//'[km]'
    ELSE IF ( .NOT. SFCFLG .AND. TANTST .LT. HGTSFC ) THEN
      ERRMSG = 'F-CHKLIM: Tangent Height=' // TANSTR // &
               '[km] < Base of Atmosphere=' // C9REAL ( HGTSFC ) // '[km]'
    ELSE
      IF ( GRAFLG ) THEN
        CALL GRACNV ( 2, ITAN )
      ELSE
        CALL TANCNV ( 2, ITAN )
      END IF
      IF ( ABS ( TAN(ITAN)%SZN ) .LT. 1.0D0 .AND. .NOT. SFCFLG ) THEN
        ERRMSG = 'F-CHKLIM: Tan.Ht.=' // TANSTR // &
                 ' [km] < atmos. after refraction, alt=' // &
                 C9REAL ( TAN(ITAN)%HGT ) // ' [km]' 
      ELSE
        FAIL = .FALSE.
      END IF
    END IF
!
  ELSE IF ( ELETAN ) THEN ! TAN(ITAN)%USR contains elevation angles.
    TAN(ITAN)%ELE = TANTST
    IF ( ABS ( TANTST ) .GT. 90.0 ) THEN
      ERRMSG = 'F-CHKLIM: Specified Elev.  angle outside range' // &
               ' +/-90deg, value=' // TANSTR
    ELSE IF ( TANTST .GT. 0.0 .AND. HGTOBS .GE. HGTTOA ) THEN       
      ERRMSG = 'F-CHKLIM: Positive Elev.Angles not allowed '// &
               'if OBServer above atmosphere'
    ELSE
      IF ( GRAFLG ) THEN
        CALL GRACNV ( 1, ITAN )
      ELSE
        CALL TANCNV ( 1, ITAN ) 
      END IF
      IF ( TAN(ITAN)%GEO .GT. HGTTOA ) THEN
        TAN(ITAN)%GEO = MIN ( TAN(ITAN)%GEO, 999999. )
        ERRMSG = 'F-CHKLIM: Elev.Ang=' // TANSTR // &
                 ' [deg] has geom.tan.pt. > atmos., alt=' //  &
                 C9REAL ( TAN(ITAN)%GEO ) // ' [km]'
      ELSE IF ( TAN(ITAN)%ELE .LT. 0.0 ) THEN ! Downward-viewing only Bug#48
        IF ( .NOT. SFCFLG .AND. TAN(ITAN)%GEO .LE. HGTSFC ) THEN
          ERRMSG = 'F-CHKLIM: Elev.Ang=' // TANSTR // &
                   ' [deg] has geom.tan.pt. < atmos., alt=' //  &
                   C9REAL ( TAN(ITAN)%GEO ) // ' [km]'
        ELSE IF ( .NOT. SFCFLG .AND. TAN(ITAN)%HGT .LE. HGTSFC ) THEN
          ERRMSG = 'F-CHKLIM: Elev.Ang=' // TANSTR // &
                   ' [deg] has refr.tan.pt. < atmos., alt=' // &
                   C9REAL ( TAN(ITAN)%HGT ) // ' [km]'
        ELSE
          FAIL = .FALSE.
        END IF
      ELSE
        FAIL = .FALSE.
      END IF
    END IF
!
  ELSE        ! TAN(ITAN)%USR contains refracted tangent heights
    TAN(ITAN)%HGT = TANTST
    IF ( TANTST .GE. HGTTOA ) THEN
      ERRMSG = 'F-CHKLIM: Tangent Height=' // TANSTR // &
               '[km] >= Top of Atmosphere=' // C9REAL ( HGTTOA ) // '[km]'
    ELSE IF ( TANTST .LT. HGTSFC ) THEN
      ERRMSG = 'F-CHKLIM: Tangent Height=' // TANSTR // &
               '[km] < Base of Atmosphere=' // C9REAL ( HGTSFC ) // '[km]'
    ELSE
      IF ( GRAFLG ) THEN
        CALL GRACNV ( 3, ITAN )
      ELSE
        CALL TANCNV ( 3, ITAN ) 
      END IF
      FAIL = .FALSE.
    END IF
  END IF
  IF ( FAIL ) RETURN
!
  TAN(ITAN)%IAT = IBRAKT ( HGTATM, TAN(ITAN)%HGT, IATM, .TRUE. )
! Could return IAT for level below if TANHGT is exactly same as atm level
  IATM = TAN(ITAN)%IAT
  IF ( TAN(ITAN)%HGT .EQ. HGTATM(IATM+1) ) TAN(ITAN)%IAT = TAN(ITAN)%IAT + 1
!
END SUBROUTINE CHKLIM
END MODULE CHKLIM_SUB

