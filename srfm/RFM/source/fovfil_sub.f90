MODULE FOVFIL_SUB
CONTAINS
SUBROUTINE FOVFIL ( NAMFOV, FAIL, ERRMSG )
!
! VERSION
!   09NOV23 AD Checked.
!   04MAY18 AD Bug#4 - read into local array ALT rather than FOV%ALT directly
!   26FEB18 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Read FOV data file
!   Called once by DRVFOV if FOV flag enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FOVCOM_DAT ! Field of View data
    USE FLGCOM_DAT, ONLY: OBSFLG ! T = observer within atmosphere
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE OPNFIL_SUB ! Open input file

  IMPLICIT NONE 
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMFOV ! Name of FOV data file
    LOGICAL,       INTENT(OUT) :: FAIL   ! T=A fatal error was detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES      
    INTEGER(I4) :: IFOV   ! Counter for FOV tabulation points
    INTEGER(I4) :: IOS    ! Saved value of IOSTAT for error message
    REAL(R4), ALLOCATABLE :: ALT(:)    ! Temporary array for input
    REAL(R8), ALLOCATABLE :: TABFOV(:) ! Tabulated FOV response function
    REAL(R8), ALLOCATABLE :: DALT(:)   ! Differences between altitudes
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  CALL OPNFIL ( LUNTMP, NAMFOV, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!      
  READ ( LUNTMP, *, IOSTAT=IOS ) NFOV
  IF ( IOS .NE. 0 ) THEN
    FAIL = .TRUE.
    WRITE ( ERRMSG, * ) 'F-FOVFIL: Error reading No.FOV pts. IOSTAT=', IOS
    RETURN
  END IF
  ELEFOV = NFOV .LT. 0
  NFOV = ABS ( NFOV )
!
  FAIL = .TRUE.
  IF ( NFOV .LE. 2 ) THEN 
    WRITE ( ERRMSG, '(A,I2,A)' ) 'F-FOVFIL: No.tabulated FOV pts=', NFOV, &
      ' - should be at least 3'
  ELSE IF ( ELEFOV .AND. .NOT. OBSFLG ) THEN 
    ERRMSG = 'F-FOVFIL: Require Observation altitude (OBS Flag)' // &
             ' to use angular FOV tabulation'
  ELSE IF ( .NOT. ELEFOV .AND. OBSFLG ) THEN 
    ERRMSG ='F-FOVFIL: Require angular FOV tabulation with OBS flag'
  ELSE
    FAIL = .FALSE.
  END IF
  IF ( FAIL ) RETURN
!
! Read list of FOV altitudes
  ALLOCATE ( FOV(NFOV), TABFOV(NFOV), DALT(NFOV-1), ALT(NFOV) )
  READ ( LUNTMP, *, IOSTAT=IOS ) ALT
  IF ( IOS .NE. 0 ) THEN 
    FAIL = .TRUE.
    WRITE ( ERRMSG, * ) 'F-FOVFIL: Error reading FOV altitudes. IOSTAT=', IOS
    RETURN
  END IF
  FOV%ALT = ALT
!
! Read tabulated FOV function 
  READ ( LUNTMP, *, IOSTAT=IOS ) TABFOV
  IF ( IOS .NE. 0 ) THEN 
    FAIL = .TRUE.
    WRITE ( ERRMSG, * ) 'F-FOVFIL: Error reading FOV function. IOSTAT=', IOS
    RETURN
  END IF
!
  CLOSE ( LUNTMP )
!
  IF ( TABFOV(1) .NE. 0.0 .OR. TABFOV(NFOV) .NE. 0.0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-FOVFIL: FOV Function non-zero at one or both ends'
    RETURN
  END IF
!
! Convert tabulated FOV points to weights for numerical integration
  DO IFOV = 1, NFOV-1
    DALT(IFOV) = DBLE ( FOV(IFOV+1)%ALT - FOV(IFOV)%ALT )
  END DO
!
  FOV(1)%FNC = DALT(1) * TABFOV(2) 
  DO IFOV = 2, NFOV-1
    FOV(IFOV)%FNC = DALT(IFOV-1) * ( TABFOV(IFOV-1) + 2 * TABFOV(IFOV) ) &
                    + DALT(IFOV) * ( TABFOV(IFOV+1) + 2 * TABFOV(IFOV) ) 
  END DO
  FOV(NFOV)%FNC = DALT(NFOV-1) * TABFOV(NFOV-1) 
!
! Normalise
  FOV%FNC = FOV%FNC / SUM ( FOV%FNC ) 
!
  FAIL = .FALSE.
!
END SUBROUTINE FOVFIL
END MODULE FOVFIL_SUB
