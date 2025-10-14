MODULE REASPC_SUB
CONTAINS
SUBROUTINE REASPC ( FILSPC, BINFIL, WNOMIN, WNOMAX, FAIL, ERRMSG )
!
! VERSION
!   30JAN24 AD Checked.
!   01SEP17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Read irregular grid data from .spc file into GRDCOM
!   Called by GRDFIL, SPCFIL, SPCFUL
!   Assumes that FILSPC is an irregular spectral file which has some
!   overlap with the WNOMIN:WNOMAX range.
!   Ensures that WNOGRD(1) = WNOMIN, WNOGRD(NGRD) = WNOMAX
!
! VARIABLE KINDS
    USE KIND_DAT  
!
! GLOBAL DATA
    USE GRDCOM_DAT ! Irregular grid
    USE FINCOM_DAT, ONLY: NOMFIN ! Nominal No.pts/cm-1 for fine mesh
    USE PHYCON_DAT, ONLY: GHZ2CM ! GHz to cm-1 conv.fac
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: FILSPC ! Name of .spc file
    LOGICAL,       INTENT(IN)  :: BINFIL ! T=binary file, F=text file
    REAL(R8),      INTENT(IN)  :: WNOMIN ! Lowest wavenumber [cm-1] to be set
    REAL(R8),      INTENT(IN)  :: WNOMAX ! Highest wavenumber [cm-1] to be set
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL       :: GHZFIL ! T=file spectal axis in GHz, F=cm-1
    INTEGER(I4)   :: IOS    ! Saved value of IOSTAT for error messages
    INTEGER(I4)   :: NWNO   ! No. spectral points in file
    REAL(R8)      :: WNO    ! Wavenumber [cm-1/GHz] read from file
    REAL(R8)      :: WNO1   ! Wavenumber [cm-1/GHz] at start of grid
    REAL(R8)      :: WNOTOL ! Tolerance [cm-1] for matching Waveno.
    CHARACTER(80) :: REC80  ! Header record from file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  REC80(1:1) = '!'
  IF ( BINFIL ) THEN
    OPEN ( LUNTMP, FILE=FILSPC, STATUS='OLD', ACTION='READ', &
           FORM='UNFORMATTED', IOSTAT=IOS, ERR=900 ) 
    DO WHILE ( REC80(1:1) .EQ. '!' ) 
      READ ( LUNTMP, IOSTAT=IOS, ERR=900 ) REC80
    END DO
    BACKSPACE ( LUNTMP ) 
    READ ( LUNTMP, IOSTAT=IOS, ERR=900 ) NWNO, WNO1
  ELSE 
    OPEN ( LUNTMP, FILE=FILSPC, STATUS='OLD', ACTION='READ', &
           IOSTAT=IOS, ERR=900 ) 
    DO WHILE ( REC80(1:1) .EQ. '!' ) 
      READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) REC80
    END DO
    READ ( REC80, *, IOSTAT=IOS, ERR=900 ) NWNO, WNO1
  END IF
!
  GHZFIL = NWNO .LT. 0
  IF ( GHZFIL ) THEN
    NWNO = ABS ( NWNO ) 
    WNO1 = WNO1 * GHZ2CM  ! Convert GHz to cm-1 
  END IF
!
! Max no grid points required will be number in file plus WNOMIN,WNOMAX
  IF ( ALLOCATED ( GRD ) ) DEALLOCATE ( GRD ) 
  ALLOCATE ( GRD(NWNO+2) ) 
!
  WNRGRD = 1.0D0 / DBLE ( NOMFIN ) 
  WNOTOL = 0.1 * WNRGRD
  GRD(1)%WNO = WNOMIN
  GRD(1)%IFN = 1
  NGRD = 1
! Skip any grid points .le. WNOMIN
  WNO = WNO1
  DO WHILE ( WNO .LT. WNOMIN+WNOTOL ) 
    IF ( BINFIL ) THEN
      READ ( LUNTMP, IOSTAT=IOS, ERR=900 ) WNO
    ELSE
      READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) WNO
    END IF
    IF ( GHZFIL ) WNO = WNO * GHZ2CM ! Convert GHz to cm-1 
  END DO
! WNO now contains second irreg grid point within reqd range
  DO WHILE ( WNO .LT. WNOMAX-WNOTOL ) 
    NGRD = NGRD + 1
    GRD(NGRD)%WNO = WNO
    GRD(NGRD)%IFN = NINT ( ( WNO - WNOMIN ) / WNRGRD ) + 1
    IF ( BINFIL ) THEN
      READ ( LUNTMP, IOSTAT=IOS, ERR=900, END=100 ) WNO
    ELSE
      READ ( LUNTMP, *, IOSTAT=IOS, ERR=900, END=100 ) WNO
    END IF
    IF ( GHZFIL ) WNO = WNO * GHZ2CM   ! Convert GHz to cm-1 
  END DO
! Add final point at upper limit of range
100  CONTINUE
  NGRD = NGRD + 1
  GRD(NGRD)%WNO = WNOMAX
  GRD(NGRD)%IFN = NINT ( ( WNOMAX - WNOMIN ) / WNRGRD ) + 1
!
  CLOSE ( LUNTMP, IOSTAT=IOS, ERR=900 )
!
900  CONTINUE
  FAIL = IOS .NE. 0
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-REASPC: I/O error reading file, IOSTAT=', IOS
!
END SUBROUTINE REASPC
END MODULE REASPC_SUB
