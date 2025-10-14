MODULE CHKPRF_SUB
CONTAINS
SUBROUTINE CHKPRF ( TYP, PROFIL, FAIL, ERRMSG )
!
! VERSION
!   05SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Check atmospheric profile on input
!   Called by ATMFIL and ATMPRF for various profiles on input
!   Requirements for each profile type are as follows
!   'HGT' all values must be monotonically increasing
!   'TEM' all values must be positive
!   'PRE' all values must be positive and monotonically decreasing
!   'VMR' all values must be zero or positive, multiply by 1.0e-6
!   'EXT' no checks.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string 
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: TYP    ! Prof type: 'HGT','TEM','PRE' etc
    REAL(R4),      INTENT(IN)  :: PROFIL(:) ! Profile
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error occurs
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: ILEV ! Profile level counter
    INTEGER(I4) :: NLEV ! No. profile levels
    REAL(R4)    :: PREV ! Previous (lower alt) profile value
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .TRUE.
  NLEV = SIZE ( PROFIL ) 
!
  SELECT CASE ( TYP ) 
  CASE ( 'HGT' ) ! Check that altitudes increase monotonically
    PREV = PROFIL(1)
    DO ILEV = 2, NLEV
      IF ( PROFIL(ILEV) .LE. PREV ) THEN
        ERRMSG = 'F-CHKPRF: HGT profile contains non-increasing ' // &
                 'altitude value=' // C9REAL ( PROFIL(ILEV) )
        RETURN
      END IF
      PREV = PROFIL(ILEV)
    END DO
  CASE ( 'TEM' ) ! Check that temperatures are all positive
    DO ILEV = 1, NLEV
      IF ( PROFIL(ILEV) .LE. 0.0 ) THEN
        ERRMSG = 'F-CHKPRF: TEM profile contains non-positive ' // &
                 'temperature value=' // C9REAL ( PROFIL(ILEV) )
        RETURN
      END IF
    END DO
  CASE ( 'PRE' ) ! Check that pressures are all +ve and decrease monotonically
    PREV = PROFIL(1) 
    DO ILEV = 2, NLEV
      IF ( PROFIL(ILEV) .LT. 0.0 ) THEN
        ERRMSG = 'F-CHKPRF: PRE profile contains negative value,' // &
                 'pressure value=' // C9REAL ( PROFIL(ILEV) )
        RETURN
      ELSE IF ( PROFIL(ILEV) .GE. PREV ) THEN
        ERRMSG = 'F-CHKPRF: PRE profile contains non-decreasing ' // &
                 'pressure value=' // C9REAL ( PROFIL(ILEV) )
        RETURN
      END IF
      PREV = PROFIL(ILEV)
    END DO
  CASE ( 'VMR' ) ! Apart from aerosol, assume VMR and check all values 0 or +ve
    DO ILEV = 1, NLEV
      IF ( PROFIL(ILEV) .LT. 0.0 ) THEN
        ERRMSG = 'F-CHKPRF: ' // TYP // ' profile contains negative VMR ' // &
                 'value=' // C9REAL ( PROFIL(ILEV) )
        RETURN
      END IF
    END DO
  CASE ( 'EXT' ) ! No checks on extinction profiles
    CONTINUE
  CASE DEFAULT
    WRITE (*,*) 'F-CHKPRF: Logical error - unrecognised argument TYP=' // TYP
    STOP
  END SELECT
!
  FAIL = .FALSE.
!
END SUBROUTINE CHKPRF
END MODULE CHKPRF_SUB
