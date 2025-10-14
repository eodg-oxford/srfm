MODULE LUTINF_SUB
CONTAINS
SUBROUTINE LUTINF ( LUNLUT, NAMLUT, BINFIL, IDXMOL, IDXISO, V1, V2, &
                    FAIL, ERRMSG )
! VERSION
!   18JUN24 AD Checked.
!   20DEC17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Extract info from TAB-LUT header record
!   Called by LUTFIL for each LUT file.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL CONSTANTS
    USE PHYCON_DAT, ONLY: GHZ2CM ! GHz to cm-1 conv.fac
!
! SUBROUTINES
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNLUT ! LUN for reading TAB-LUT file
    CHARACTER(*),  INTENT(IN)  :: NAMLUT ! Name of TAB-LUT file
    LOGICAL,       INTENT(OUT) :: BINFIL ! T=Read as binary file, F=ASCII
    INTEGER(I4),   INTENT(OUT) :: IDXMOL ! RFM/HITRAN ID of gas
    INTEGER(I4),   INTENT(OUT) :: IDXISO ! Isotopic ID of gas, or -1
    REAL(R8),      INTENT(OUT) :: V1     ! Lower Wno limit [cm-1] of data
    REAL(R8),      INTENT(OUT) :: V2     ! Upper Wno limit [cm-1] of data
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IOS    ! Saved value of IOSTAT for error messages
    INTEGER(I4)   :: ITRY   ! Counter for different attempts to open/read file
    INTEGER(I4)   :: NV     ! No.spectral points in TAB-LUT file
    REAL(R4)      :: RFMT   ! File format identifier
    REAL(R4)      :: RGAS   ! Molec.isotope encoded as real number
    CHARACTER(5)  :: GASSTR ! Part of RECHDR containing molec ID
    CHARACTER(80) :: RECORD ! Record read from TAB-LUT file
    CHARACTER(80) :: RECHDR ! Record containing header info
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IOS = 0
!
! Initially guess BINFIL=T if '.bin' string present in filename
  BINFIL = INDEX ( NAMLUT, '.bin' ) .GT. 0 
!
  DO ITRY = 1, 2 
    IF ( ITRY .EQ. 2 ) BINFIL = .NOT. BINFIL   ! 2nd time around, reverse
    IF ( BINFIL ) THEN
      OPEN ( UNIT=LUNLUT, FILE=NAMLUT, STATUS='OLD', ACTION='READ', &
             FORM='UNFORMATTED', IOSTAT=IOS, ERR=800 ) 
      READ ( LUNLUT, IOSTAT=IOS, ERR=800 ) RECORD
      IF ( RECORD(1:1) .NE. '!' ) GOTO 800     ! Not binary LUT file
      CALL WRTLOG ( RECORD ) 
      DO WHILE ( RECORD(1:1) .EQ. '!' ) 
        READ ( LUNLUT, IOSTAT=IOS, ERR=800 ) RECORD
      END DO
      READ ( LUNLUT, IOSTAT=IOS, ERR=800 ) RECHDR
    ELSE
      OPEN ( UNIT=LUNLUT, FILE=NAMLUT, STATUS='OLD', ACTION='READ', &
             IOSTAT=IOS, ERR=800 )
      READ ( LUNLUT, '(A)', IOSTAT=IOS, ERR=800 ) RECORD
      IF ( RECORD(1:1) .NE. '!' ) GOTO 800    ! Not ASC LUT file
      CALL WRTLOG ( RECORD ) 
      DO WHILE ( RECORD(1:1) .EQ. '!' ) 
        READ ( LUNLUT, '(A)', IOSTAT=IOS, ERR=800 ) RECORD
      END DO
      READ ( LUNLUT, '(A)', IOSTAT=IOS, ERR=800 ) RECHDR
    END IF
800 CONTINUE
    CLOSE ( LUNLUT, IOSTAT=IOS, ERR=900 )
    IF ( IOS .EQ. 0 ) EXIT
  END DO
  IF ( IOS .NE. 0 ) GOTO 900
!
  READ ( RECORD, *, IOSTAT=IOS, ERR=900 ) RFMT   
  IF ( RFMT .NE. 1.0 ) THEN
    WRITE ( ERRMSG, '(A,F10.2)' ) &
       'F-LUTINF: unrecognised LUT file format, version=', RFMT
    FAIL = .TRUE.
    RETURN
  END IF
!
  GASSTR = RECHDR(1:5)
  IF ( INDEX ( GASSTR, '.' ) .EQ. 0 ) THEN
    READ ( GASSTR, *, IOSTAT=IOS, ERR=900 ) IDXMOL
    IDXISO = -1
  ELSE
    READ ( GASSTR, *, IOSTAT=IOS, ERR=900 ) RGAS
    IDXMOL = INT ( RGAS )
    IDXISO = NINT( 10*(RGAS-IDXMOL) )
  END IF
!
  READ ( RECHDR(6:), *, IOSTAT=IOS, ERR=900 ) NV, V1, V2
  IF ( NV .LT. 0 ) THEN         ! Convert from GHz to cm-1
    NV = -NV 
    V1 = V1 * GHZ2CM
    V2 = V2 * GHZ2CM
  END IF
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-LUTINF: Error reading LUT file header. IOSTAT=', IOS
!
END SUBROUTINE LUTINF
END MODULE LUTINF_SUB
