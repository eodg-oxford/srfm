MODULE ATMFIL_SUB
CONTAINS
SUBROUTINE ATMFIL ( NAMATM, FAIL, ERRMSG )
!
! VERSION
!   18MAR21 AD Add more log messages. Checked.
!   16AUG19 AD Add ATMQAL. Checked.
!   21JUN17 AD Read in entire profile, allow for PRE or HGT grids.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Read file containing atmospheric profiles
!   Called by DRVATM for each file listed in *ATM section of driver table.
!   The first profile sets the atmospheric grid for subsequent profiles
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT, ONLY: LENATM ! Max length of .atm profile label
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE ATMGRD_SUB ! Set atmospheric profile vertical grid
    USE ATMPRF_SUB ! Load profile into ATMCOM
    USE ATMQAL_SUB ! Extract qualifiers from .atm filename
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE C9REAL_GEN ! Write real number as C*9 string
    USE CHKPRF_SUB ! Check atmospheric profile on input
    USE NXTPRF_SUB ! Get next profile from .atm file
    USE OPNFIL_SUB ! Open ASCII input file, and skip/Log any initial comments.
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(INOUT) :: NAMATM ! Name of file 
    LOGICAL,       INTENT(OUT)   :: FAIL   ! TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT)   :: ERRMSG ! Error message written if FAIL=t
!
! LOCAL VARIABLES
    LOGICAL       :: ANYUSE         ! T = at least one profile used
    LOGICAL       :: FIRST = .TRUE. ! T=first profile
    LOGICAL       :: PREGRD = .FALSE. ! T = pressure, F = height grid
    LOGICAL       :: USEPRF         ! T = current profile was used 
    INTEGER(I4)   :: IDXHGT         ! Index of *HGT profile in PRFDAT
    INTEGER(I4)   :: IDXPRE         ! Index of *PRE profile in PRFDAT
    INTEGER(I4)   :: IOS            ! Saved value IOSTAT for error
    INTEGER(I4)   :: IPRF           ! Counter for profiles in file
    INTEGER(I4)   :: NLEV           ! No. levels in supplied profiles
    INTEGER(I4)   :: NPRF           ! No. profiles in file
    INTEGER(I4)   :: NQAL           ! No. qualifers attached to filename
    CHARACTER(20) :: LABEL          ! Label identifying profile
    REAL(R4),          ALLOCATABLE :: PRFDAT(:,:) ! All profiles in file
    REAL(R4),          ALLOCATABLE :: PRFLEV(:)   ! Input profile values       
    REAL(R4),          ALLOCATABLE :: PROFIL(:)   ! Input profile values       
    REAL(R4),          ALLOCATABLE :: SAVDAT(:,:) ! Saved PRFDAT during realloc
    CHARACTER(20),     ALLOCATABLE :: PRFLAB(:)   ! Profile labels
    CHARACTER(20),     ALLOCATABLE :: SAVLAB(:)   ! Saved PRFLAB during realloc
    CHARACTER(LENATM), ALLOCATABLE :: QALPRF(:)   ! List of qualifiers
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL ATMQAL ( NAMATM, NQAL, QALPRF, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
  CALL OPNFIL ( LUNTMP, NAMATM, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN
  ANYUSE = FIRST    ! First file will always be used (for vertical grid)
!
! 1st non-comment record in file must be no. levels (integer)
  READ ( LUNTMP, *, IOSTAT=IOS ) NLEV
  IF ( IOS .NE. 0 ) THEN
    FAIL = .TRUE.
    WRITE ( ERRMSG, * ) &
      'F-ATMFIL: Failed to read No.profile levels from file. IOSTAT=', IOS
    RETURN
  END IF
!
  ALLOCATE ( PRFLEV(NLEV) )
  ALLOCATE ( PROFIL(NLEV) )
  ALLOCATE ( PRFDAT(NLEV,1) )
  ALLOCATE ( PRFLAB(1) )
!
! Read in all profiles from this file
  NPRF = 0
  IDXPRE = 0
  IDXHGT = 0
  DO
    CALL NXTPRF ( LUNTMP, LABEL, PROFIL, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( LABEL .EQ. '*END' ) EXIT
!
! List of specific profiles for inclusion/exclusion appended to filename
    IF ( NQAL .NE. 0 ) THEN
      IF ( NQAL .LT. 0 ) THEN                      ! exclusion list
        USEPRF = .NOT. ANY ( QALPRF .EQ. LABEL )         
      ELSE                                         ! inclusion list
        USEPRF = ANY ( QALPRF .EQ. LABEL )         
! Logic is complicated by the requirement that once a grid is established
! (.not. first) the user shouldn't have to specify the grid profile in the
! inclusion list
        IF ( .NOT. FIRST ) THEN
          IF ( PREGRD ) THEN 
            IF ( LABEL .EQ. '*PRE' ) USEPRF = .TRUE.
          ELSE
            IF ( LABEL .EQ. '*HGT' ) USEPRF = .TRUE.
          END IF
        END IF
      END IF
      IF ( .NOT. USEPRF ) CYCLE
    END IF
!
    NPRF = NPRF + 1
    IF ( NPRF .GT. 1 ) THEN
      CALL MOVE_ALLOC ( PRFLAB, SAVLAB )
      ALLOCATE ( PRFLAB(NPRF) ) 
      PRFLAB(1:NPRF-1) = SAVLAB
      CALL MOVE_ALLOC ( PRFDAT, SAVDAT )
      ALLOCATE ( PRFDAT(NLEV,NPRF) )
      PRFDAT(:,1:NPRF-1) = SAVDAT
    END IF
    SELECT CASE ( LABEL )
    CASE ( '*PRE' ) 
      IDXPRE = NPRF
      CALL CHKPRF ( 'PRE', PROFIL, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    CASE ( '*HGT' )
      IDXHGT = NPRF
      CALL CHKPRF ( 'HGT', PROFIL, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    CASE DEFAULT
      CONTINUE
    END SELECT
!
! First profile of first .atm file sets vertical grid
    IF ( FIRST ) THEN
      CALL ATMGRD ( LABEL, NLEV, PROFIL, PREGRD, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      IF ( PREGRD ) THEN
        CALL WRTLOG ( 'I-ATMFIL: Set Pressure grid: ' // &
          TRIM ( C11INT(NLEV) ) // ' profile levels from ' // & 
          C9REAL(PROFIL(1)) // ':' // &
          C9REAL(PROFIL(NLEV)) // ' [hPa]' )
      ELSE 
        CALL WRTLOG ( 'I-ATMFIL: Set Altitude grid: ' // &
          TRIM ( C11INT(NLEV) ) // ' profile levels from ' // & 
          C9REAL(PROFIL(1)) // ':' // &
          C9REAL(PROFIL(NLEV)) // ' [km]' )
      END IF
      FIRST = .FALSE.
    END IF
!
    PRFLAB(NPRF) = LABEL
    PRFDAT(:,NPRF) = PROFIL
  END DO
!
  IF ( PREGRD ) THEN
    IF ( IDXPRE .EQ. 0 ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-ATMFIL: No *PRE profile in .atm file' 
      RETURN
    END IF
    PRFLEV = LOG ( PRFDAT(:,IDXPRE) )      
  ELSE
    IF ( IDXHGT .EQ. 0 ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-ATMFIL: No *HGT profile in .atm file' 
      RETURN
    END IF
    PRFLEV = PRFDAT(:,IDXHGT)      
  END IF
!
  DO IPRF = 1, NPRF
    IF ( IPRF .EQ. IDXPRE .AND. PREGRD ) CYCLE
    IF ( IPRF .EQ. IDXHGT .AND. .NOT. PREGRD ) CYCLE
    CALL ATMPRF ( PRFLAB(IPRF), PRFLEV, PRFDAT(:,IPRF), USEPRF, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    ANYUSE = ANYUSE .OR. USEPRF
  END DO
!
  CLOSE ( LUNTMP )
!
  IF ( .NOT. ANYUSE ) &
    CALL WRTLOG ( 'W-ATMFIL: No useful profiles found in file' )
!
  FAIL = .FALSE.
!
END SUBROUTINE ATMFIL
END MODULE ATMFIL_SUB

