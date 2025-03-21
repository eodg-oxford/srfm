MODULE SPCRNG_SUB
CONTAINS
SUBROUTINE SPCRNG ( NAMGRD, GHZGRD, TYPGRD, NPT, PT1, PT2, PTD, FAIL, ERRMSG )
!
! VERSION
!   28JUN24 AD Checked.
!   30DEC17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Extract spectral range/resln from file
!   Called by SPCFIL, GRDFIL
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE PHYCON_DAT, ONLY: GHZ2CM ! GHz to cm-1 conv.fac (~1/30)
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMGRD ! Name of GRD file 
    LOGICAL,       INTENT(IN)  :: GHZGRD ! T=range in GHz, F=range in cm-1
    CHARACTER(3),  INTENT(OUT) :: TYPGRD ! GRD file type ('GRD','SPC','BIN')
    INTEGER(I4),   INTENT(OUT) :: NPT    ! No. spectral points in file
    REAL(R8),      INTENT(OUT) :: PT1    ! First spectral point in file
    REAL(R8),      INTENT(OUT) :: PTD    ! Spectral increment (or 0 if irreg.)
    REAL(R8),      INTENT(OUT) :: PT2    ! Last spectral point in file
    LOGICAL,       INTENT(OUT) :: FAIL   ! T=Fatal error detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: NTYP = 3   ! No different GRD file types
!
! LOCAL VARIABLES
    LOGICAL       :: GHZRNG       ! T=file grid in GHz, F=grid in cm-1.
    INTEGER(I4)   :: IOS          ! Saved value of IOSTAT for error messages
    INTEGER(I4)   :: ITYP         ! Counter for different grid file types
    INTEGER(I4)   :: NTOT         ! Total no.full grid points in .grd file
    CHARACTER(80) :: REC80        ! Record from GRD file
    CHARACTER(3)  :: TYPLST(NTYP) ! List of different GRD types
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Guess probable type and put first in list to test
  IF ( INDEX ( NAMGRD, '.grd' ) .GT. 0 ) THEN
    TYPLST = (/ 'GRD', 'SPC', 'BIN' /)
  ELSE IF ( INDEX ( NAMGRD, '.bin' ) .GT. 0 ) THEN
    TYPLST = (/ 'BIN', 'SPC', 'GRD' /)
  ELSE
    TYPLST = (/ 'SPC', 'GRD', 'BIN' /)
  END IF
!
  DO ITYP = 1, NTYP
    TYPGRD = TYPLST(ITYP)
    SELECT CASE ( TYPGRD ) 
    CASE ( 'GRD' ) 
      OPEN ( UNIT=LUNTMP, FILE=NAMGRD, STATUS='OLD', ACTION='READ', &
             IOSTAT=IOS, ERR=100 ) 
      READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=100 ) REC80
      REC80(1:1) = '!'    ! 1st record may not start with '!' character
      DO WHILE ( REC80(1:1) .EQ. '!' )
        READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=100 ) REC80 
      END DO
      IF ( REC80(1:3) .NE. 'lin' ) GOTO 100
      READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) NTOT, NPT, PT1, PTD
      GHZRNG = NTOT .LT. 0
      PT2 = PT1 + ( ABS(NTOT) - 1 ) * PTD
      PTD = 0.0       ! flag to indicate irregular grid
    CASE ( 'SPC' ) 
      OPEN ( UNIT=LUNTMP, FILE=NAMGRD, STATUS='OLD', ACTION='READ', &
             IOSTAT=IOS, ERR=100 ) 
      READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=100 ) REC80
      DO WHILE ( REC80(1:1) .EQ. '!' )  
        READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=100 ) REC80
      END DO
      READ ( REC80, *, IOSTAT=IOS, ERR=100 ) NPT, PT1, PTD, PT2
      GHZRNG = NPT .LT. 0
      NPT = ABS ( NPT ) 
    CASE ( 'BIN' ) 
      OPEN ( UNIT=LUNTMP, FILE=NAMGRD, STATUS='OLD', ACTION='READ', &
             FORM='UNFORMATTED', IOSTAT=IOS, ERR=100 ) 
      READ ( LUNTMP, IOSTAT=IOS, ERR=100 ) REC80
      DO WHILE ( REC80(1:1) .EQ. '!' )  
        READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=100 ) REC80
      END DO
      BACKSPACE ( LUNTMP, IOSTAT=IOS, ERR=100 ) 
      READ ( LUNTMP, IOSTAT=IOS, ERR=100 ) NPT, PT1, PTD, PT2
      GHZRNG = NPT .LT. 0
      NPT = ABS ( NPT ) 
    CASE DEFAULT
      STOP 'F-SPCRNG: Logical error'   ! unknown TYP
    END SELECT
!
! If this point reached, spectral range/resl succesfull read
    IF ( GHZRNG .AND. .NOT. GHZGRD ) THEN  ! convert GHz to cm-1
      PT1 = PT1 * GHZ2CM
      PT2 = PT2 * GHZ2CM
      PTD = PTD * GHZ2CM
    ELSE IF ( .NOT. GHZRNG .AND. GHZGRD ) THEN ! convert cm-1 to GHz
      PT1 = PT1 / GHZ2CM
      PT2 = PT2 / GHZ2CM
      PTD = PTD / GHZ2CM
    END IF
    CLOSE ( LUNTMP ) 
    FAIL = .FALSE. 
    RETURN           ! Normal exit
!
! Jump here if file read error, assumed due to wrong type
100 CONTINUE   
    CLOSE ( LUNTMP )
  END DO
!
! If this point is reached, unable to successfully read any type of file
900 CONTINUE
  TYPGRD = ' '
  FAIL = .TRUE.
  WRITE ( ERRMSG, * ) &
    'F-SPCRNG: Error reading spectral range file, IOSTAT=', IOS
!
END SUBROUTINE SPCRNG
END MODULE SPCRNG_SUB
