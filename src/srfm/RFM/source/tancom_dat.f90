MODULE TANCOM_DAT
!
! VERSION
!   25MAY24 AD Checked
!   25MAR19 AD Add UNITAN, USRUNI. Checked.
!   05MAR19 AD Add TAN%SEC, TAN%SKY, TAN%ISK
!   30MAY18 AD Make TAN a TARGET
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Tangent path data.
!   Output tangent heights are 1:NTAN (1=lowest, NTAN=highest) limb view
!   also sorted for NAD, ZEN or HOM views.
!   Extra tangent heights for FOV calculation are in NTAN+1:MTAN
!   More may be added for Jacobian calculations
!   Reset via TANCOM_RESET for repeated library invocations.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
  PUBLIC :: TANCOM_RESET
!
! GLOBAL CONSTANTS
    INTEGER, PARAMETER :: LENTAN = 6 ! Max length of tan.ht string in filenames
!
  TYPE :: TANTYP
    LOGICAL     :: CLC ! T = Radiance.calc reqd for Tan.Hgt
    LOGICAL     :: SFC ! T = limb ray intersects surface
    LOGICAL     :: SKY ! T = use for diffuse sky radiance
    INTEGER(I4) :: IAT ! Atmospheric layer containing tangent pt 
    INTEGER(I4) :: ISK ! Index of tan for diffuse sky radiance component, or 0
    INTEGER(I4) :: ITN ! Index of original tan (or 0 if this is original).
    INTEGER(I4) :: JDX ! Target species for Jacobian perturbations, or 0
    REAL(R4)    :: ELE ! Elevation angles of paths [deg]
    REAL(R4)    :: GEO ! Projected (=Geom) tan.alt. [km]
    REAL(R4)    :: HGT ! Lowest altitude along paths [km] 
    REAL(R4)    :: PSI ! LOS angle of tan.pt. [deg] (with GRA flag)
    REAL(R4)    :: SEC ! Sec(zenith angle) for p/p atmospheres
    REAL(R4)    :: USR ! Tan.path values from Drv.Table
    REAL(R8)    :: SZN ! Sin(Zenith angle) at 'tangent point'/obs/surface
    CHARACTER(LENTAN) :: STR ! Strings for filenames
  END TYPE TANTYP
!
! GLOBAL VARIABLES
    TYPE(TANTYP), TARGET, ALLOCATABLE :: TAN(:) 
!
    LOGICAL      :: USRELE = .FALSE. ! T=USRTAN are elevation angles
    LOGICAL      :: USRGEO = .FALSE. ! T=USRTAN are Geom.Tan.points.
    LOGICAL      :: LIMTAN        ! T=limb-viewing geom., F=plane or homog.path
    INTEGER(I4)  :: MTAN          ! No. ray paths reqd for calculation
    INTEGER(I4)  :: NTAN          ! No. nominal tan pts for output
    REAL(R4)     :: UNITAN = 1.0  ! Scale factor if units not km
    CHARACTER(2) :: USRUNI = 'km' ! Units for scale factor
!
CONTAINS

  SUBROUTINE TANCOM_RESET()
    IF ( ALLOCATED ( TAN ) ) DEALLOCATE ( TAN )
    USRELE = .FALSE.
    USRGEO = .FALSE.
    LIMTAN = .FALSE.
    MTAN   = 0
    NTAN   = 0
    UNITAN = 1.0_R4
    USRUNI = 'km'
  END SUBROUTINE TANCOM_RESET

END MODULE TANCOM_DAT
