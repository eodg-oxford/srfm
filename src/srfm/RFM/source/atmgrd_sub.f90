MODULE ATMGRD_SUB
CONTAINS
SUBROUTINE ATMGRD ( LABEL, NLEV, PROFIL, PREGRD, FAIL, ERRMSG )
!
! VERSION
!   02APR21 AD Bug#27 Check that Hgt grid used with GRA flag. Checked.
!   21JUN17 AD Original. Checked. 
!
! DESCRIPTION
!   Set atmospheric profile vertical grid
!   Called by ATMFIL for first profile of first file.
!   First profile should either be '*HGT' or '*PRE', unless HOM flag is enabled
!   (homogeneous path calculation - no vertical grid required).
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
!
! SUBROUTINES
    USE ATMINI_SUB ! Initialise atmospheric profile data in ATMCOM
    USE ATMLAY_SUB ! Subdivide atmospheric profile layers
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: LABEL  ! Profile label
    INTEGER(I4),   INTENT(IN)  :: NLEV   ! No. levels in profile
    REAL(R4),      INTENT(IN)  :: PROFIL(:) ! Profile
    LOGICAL,       INTENT(OUT) :: PREGRD ! T=pressure grid, F=altitude grid
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    REAL(R4), ALLOCATABLE :: LEVGRD(:) ! Profile grid levels
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( NLEV .EQ. 1 .AND. .NOT. HOMFLG ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-ATMGRD: Single atmos.level specified but ' // &
             'HOM flag (homog. path calc) not set'
    RETURN
  ELSE IF ( HOMFLG .AND. NLEV .NE. 1 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-ATMGRD: HOM flag set (homog. path calc) but' // &
             ' multi-layer atmos.profile'
    RETURN
  END IF
!
  ALLOCATE ( LEVGRD(NLEV) ) 
  LEVGRD = PROFIL
!
  SELECT CASE ( LABEL ) 
  CASE ( '*HGT' ) 
    PREGRD = .FALSE.
  CASE ( '*PRE' ) 
    PREGRD = .TRUE.
    IF ( GRAFLG ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-ATMGRD: GRA flag requires all atmos profiles on Height grid'
      RETURN
    END IF
  CASE DEFAULT
    IF ( HOMFLG ) THEN   
      PREGRD = .FALSE.
      LEVGRD(1) = 0.0  ! Arbitrarily assign HGT=0.0km for HOM paths
    ELSE
      FAIL = .TRUE.
      ERRMSG = 'F-ATMGRD: First profile must be *HGT or *PRE, not: ' // LABEL
      RETURN
    END IF
  END SELECT
  IF ( LAYFLG ) CALL ATMLAY ( LEVGRD, PREGRD ) 
  CALL ATMINI ( LEVGRD, PREGRD )   
!
  FAIL = .FALSE.
!
END SUBROUTINE ATMGRD
END MODULE ATMGRD_SUB
