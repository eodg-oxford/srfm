MODULE REAGRD_SUB
CONTAINS
SUBROUTINE REAGRD ( FILGRD, WNOMIN, WNOMAX, FAIL, ERRMSG )
!
! VERSION
!   21NOV18 AD Bug#11: Set upper limit of loop to MIN(IDXMAX,NTOT). Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Read irregular grid data from .grd file into GRDCOM
!   Called by GRDFIL, SPCFIL, SPCFUL
!   This sets up WNOGRD as
!     WNOGRD(1) = WNOMIN
!     WNOGRD(2) = WNOMIN+WNRGRD if 4-pt interpolation function to be used
!     WNOGRD(...) = values from .grd file
!     WNOGRD(N-1) = WNOMAX-WNRGRD if 4-pt interpolation function to be used
!     WNOGRD(N) = WNOMAX
!   This ensures that the WNOGRD array will contain grid points at the 
!   specified wno limits and, if 4-pt interpolation function to be used, at
!   the next points as well. 
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE GRDCOM_DAT ! Irregular grid
    USE PHYCON_DAT, ONLY: GHZ2CM ! GHz to cm-1 conv.fac [~1/30]
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: FILGRD ! Name of .grd file
    REAL(R8),      INTENT(IN)  :: WNOMIN ! Lowest wavenumber [cm-1] to be set
    REAL(R8),      INTENT(IN)  :: WNOMAX ! Highest wavenumber [cm-1] to be set
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL CONSTANTS
    CHARACTER(3), PARAMETER :: FNCLST(10) = &
      (/ 'lin', 'qad', 'cub', '1li', '1qa', '1cu', '1sq', 'lor', 'lnl', 'lnc' /)
!
! LOCAL VARIABLES
    INTEGER(I4)   :: ICHAR  ! Counter for hex characters within record (1:50)
    INTEGER(I4)   :: IDXMAX ! Inde of Highest grid pt to be used from grd file
    INTEGER(I4)   :: IDXMIN ! Inde of Lowest  grid pt to be used from grd file
    INTEGER(I4)   :: IHEX   ! Integer(I4) value of hexadecimal character
    INTEGER(I4)   :: IOFF   ! Offset between spectral grid and irreg grid file
    INTEGER(I4)   :: IOS    ! Saved value of IOSTAT for error messages
    INTEGER(I4)   :: ITOT   ! Counter for regular grid points
    INTEGER(I4)   :: ITWO   ! Factors of two for selecting hex bits
    INTEGER(I4)   :: NTOT   ! Total no. of regular grid points,-ve=GHz
    INTEGER(I4)   :: NUSE   ! Number of irreg.grid points represented
    REAL(R8)      :: WNO1   ! Wavenumber [cm-1/GHz] at start of reg.grid
    CHARACTER(80) :: CDUMMY ! Header record read from GRD file (dummy)
    CHARACTER(25) :: CONFIG ! GRD file configuration info
    CHARACTER(50) :: REC50  ! Grid record read from file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  OPEN ( LUNTMP, FILE=FILGRD, STATUS='OLD', ACTION='READ', IOSTAT=IOS, ERR=900 ) 
  READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) CONFIG
  READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) CDUMMY
  READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) FNCGRD
  IF ( .NOT. ANY ( FNCLST .EQ. FNCGRD ) ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-REAGRD: Unrecognised interpolation method='//FNCGRD
    RETURN
  END IF
!
  READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) NTOT, NUSE, WNO1, WNRGRD
  READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) CDUMMY
!
! Convert spectral axis from GHz to cm-1
  IF ( NTOT .LT. 0 ) THEN        
    NTOT = -NTOT
    WNO1 = WNO1 * GHZ2CM
    WNRGRD = WNRGRD * GHZ2CM
  END IF
!
  TWOGRD = INDEX ( 'lin 1li 1sq lor lnl', FNCGRD ) .NE. 0
!
! Max number of grid points that might need to be stored is number in .grd file 
! (NUSE) plus up to two extra points inserted at each end WNOMIN and WNOMAX
  IF ( ALLOCATED ( GRD ) ) DEALLOCATE ( GRD )
  ALLOCATE ( GRD(NUSE+4) ) 
!
  GRD(1)%WNO = WNOMIN
  GRD(1)%IFN = 1
  NGRD = 1
  IDXMIN = NINT ( ( WNOMIN - WNO1 ) / WNRGRD ) + 2
  IDXMAX = NINT ( ( WNOMAX - WNO1 ) / WNRGRD ) 
  IF ( .NOT. TWOGRD ) THEN
    GRD(2)%WNO = WNOMIN + WNRGRD
    GRD(2)%IFN = 2
    IDXMIN = IDXMIN + 1
    IDXMAX = IDXMAX - 1
    NGRD = 2
  END IF
!
  IOFF = NINT ( ( WNOMIN - WNO1 ) / WNRGRD ) 
  DO ITOT = 1, MIN ( IDXMAX, NTOT )                  ! Bug#11
    IF ( MOD ( ITOT, 200 ) .EQ. 1 ) THEN
      READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) REC50
      ICHAR = 0
    END IF
    IF ( MOD ( ITOT, 4 ) .EQ. 1 ) THEN
      ICHAR = ICHAR + 1
      READ ( REC50(ICHAR:ICHAR), '(Z1)', IOSTAT=IOS, ERR=900 ) IHEX
      ITWO = 16
    ENDIF
    ITWO = ITWO/2
    IF ( IHEX .GE. ITWO ) THEN
      IF ( ITOT .GE. IDXMIN ) THEN
        NGRD = NGRD + 1
        GRD(NGRD)%WNO = WNO1 + ( ITOT - 1 ) * WNRGRD
        GRD(NGRD)%IFN = ITOT - IOFF
      END IF
      IHEX = IHEX - ITWO
    END IF
  END DO
!
  IF ( .NOT. TWOGRD ) THEN
    NGRD = NGRD + 1
    IDXMAX = IDXMAX + 1
    GRD(NGRD)%WNO = WNOMAX - WNRGRD
    GRD(NGRD)%IFN = IDXMAX 
  END IF        
  NGRD = NGRD + 1
  IDXMAX = IDXMAX + 1
  GRD(NGRD)%WNO = WNOMAX
  GRD(NGRD)%IFN = IDXMAX 
!
  CLOSE ( LUNTMP, IOSTAT=IOS, ERR=900 )
!
900  CONTINUE
  FAIL = IOS .NE. 0
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-REAGRD: I/O error reading file, IOSTAT=', IOS
!
END SUBROUTINE REAGRD
END MODULE REAGRD_SUB
