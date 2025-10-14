MODULE SPCFIL_SUB
CONTAINS
SUBROUTINE SPCFIL ( NAMSPC, FAIL, ERRMSG )
!
! VERSION
!   27JUN24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Load spectral range/resln from file
!   Called by DRVSPC if filename specified in *SPC section
!   If regular grid, just copies range/resl to SPCCOM.
!   If irregular grid, saves filename in GFLCOM.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE SPCCOM_DAT ! Spectral range data
    USE FLGCOM_DAT, ONLY: GHZFLG ! T = Use GHz spectral axis, F = Wavenumber
    USE GRDCOM_DAT, ONLY: ISPGRD ! Spc.Rng index of currently loaded grid
    USE PHYCON_DAT, ONLY: GHZ2CM ! GHz to cm-1 conv.factor (1e7/c)
!
! SUBROUTINES
    USE ADDGFL_SUB ! Add irreg. grid file to GFLCOM
    USE REAGRD_SUB ! Read irregular grid data from .grd file into GRDCOM
    USE REASPC_SUB ! Read irregular grid data from .spc file into GRDCOM
    USE SPCRNG_SUB ! Extract spectral range/resln from file
!
  IMPLICIT NONE 
!
! ARGUMENTS 
    CHARACTER(*),  INTENT(IN)  :: NAMSPC ! Name of spectral/grid file
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)  :: IGFL ! Index of file loaded into GFLCOM
    INTEGER(I4)  :: NPT  ! No. spectral points in file
    REAL(R8)     :: PT1  ! First spectral point in file
    REAL(R8)     :: PT2  ! Last spectral point in file
    REAL(R8)     :: PTD  ! Spectral increment (or 0 if irreg.)
    CHARACTER(3) :: TYP  ! Type of spc.range data file
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! Extract spectral range/resl from file header
  CALL SPCRNG ( NAMSPC, GHZFLG, TYP, NPT, PT1, PT2, PTD, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN   
!           
  SPC(NSPC)%NPT = NPT
  SPC(NSPC)%WNL = PT1
  SPC(NSPC)%WNU = PT2
  SPC(NSPC)%WNR = PTD
! From this point onwards PT1,PT2,PTD are always wavenumber
  IF ( GHZFLG ) THEN
    PT1 = PT1 * GHZ2CM
    PT2 = PT2 * GHZ2CM
    PTD = PTD * GHZ2CM
  END IF
!
! For irregular grids, save filename in GFLCOM
  IF ( PTD .EQ. 0.0D0 ) THEN
    CALL ADDGFL ( NAMSPC, TYP, PT1, PT2, IGFL ) 
    SPC(NSPC)%IGD = IGFL
! Check that file can be successfully read
    IF ( TYP .EQ. 'GRD' ) THEN
      CALL REAGRD ( NAMSPC, PT1, PT2, FAIL, ERRMSG )
    ELSE
      CALL REASPC ( NAMSPC, (TYP .EQ. 'BIN'), PT1, PT2, FAIL, ERRMSG )
    END IF
    IF ( FAIL ) RETURN
    ISPGRD = NSPC
  ELSE
    SPC(NSPC)%IGD = 0
  END IF
!
! Normal exit
  FAIL = .FALSE.
!
END SUBROUTINE SPCFIL
END MODULE SPCFIL_SUB
