MODULE XSCFIL_SUB
CONTAINS
SUBROUTINE XSCFIL ( NAMXSC, FAIL, ERRMSG )
!
! VERSION
!   19FEB20 AD Simplified logic. Checked.
!   17JAN18 AD BACKSPACE rather than REWIND in case RFM format has comment recs
!   20NOV17 AD should really check that .xsc file has correct format at first
!   23JUN17 AD Allow for original HITRAN format files
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Check if .xsc file required
!   Called by DRVXSC for each filename in *XSC section, and XSCDEF if wildcard.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE SPCCOM_DAT ! Spectral range data
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
    USE MOLIDX_SUB ! Give molecule name for HITRAN/RFM index, or vice-versa
    USE OPNFIL_SUB ! Open input file
    USE REAXSC_SUB ! Read data from .xsc file
    USE SAVXSC_SUB ! Load X/S data into XSCCOM
    USE WRTLOG_SUB ! Write text message to log file
!   
  IMPLICIT NONE 
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMXSC ! Name of .xsc file 
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL       :: NODATA ! T=no data within reqd spectral range
    LOGICAL       :: RFMFMT ! T=file in RFM format, F=original HITRAN format
    INTEGER(I4)   :: IDXMOL ! HITRAN/RFM Molecule ID
    INTEGER(I4)   :: IGAS   ! Index of molecule in GASCOM
    INTEGER(I4)   :: IOS    ! Saved value of IOSTAT for error message
    INTEGER(I4)   :: NXST   ! No. (p,T) datasets in spectral range
    REAL(R8)      :: WNL    ! Lower wavenumber [cm-1] of spectral range
    REAL(R8)      :: WNU    ! Upper wavenumber [cm-1] of spectral range
    CHARACTER(20) :: C20    ! First 20 characters of first record in file
    CHARACTER(10) :: FMT    ! .xsc file format string
    CHARACTER(20) :: MOLEC  ! Right adjusted name of molec in orig.fmt file 
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  CALL OPNFIL ( LUNTMP, NAMXSC, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
! Check if RFM format or original format
  READ ( LUNTMP, '(A20)', IOSTAT=IOS, ERR=900 ) C20
  RFMFMT = C20(1:3) .EQ. '***'
  BACKSPACE ( LUNTMP )
!
! Identify molecule from first tabulation header
  IF ( RFMFMT ) THEN
    READ ( LUNTMP, '(3X,I7,I10,A10,2F10.4)', IOSTAT=IOS, ERR=900 ) &
      IDXMOL, NXST, FMT, WNL, WNU
    IF ( INDEX ( FMT, 'HITRAN2K' ) .EQ. 0 ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-XSCFIL: .xsc file is not HITRAN2k format'
      RETURN
    END IF
    BACKSPACE ( LUNTMP ) 
    CALL MOLIDX ( IDXMOL, MOLEC )
  ELSE
    MOLEC = ADJUSTL(C20)
    IDXMOL = 0
    CALL MOLIDX ( IDXMOL, MOLEC )  
  END IF

! Check if this molecule is required
  IGAS = 0
  IF ( IDXMOL .GT. 0 ) IGAS = IDXGAS ( IDXMOL )
  IF ( IGAS .EQ. 0 ) THEN
    CLOSE ( LUNTMP )
    CALL WRTLOG ('W-XSCFIL: File ignored - Molecule not required')
    RETURN
  END IF
  GAS(IGAS)%XSC = .TRUE.     ! At least one x/s file identified for molec.
!
  CALL REAXSC ( LUNTMP, RFMFMT, MOLEC, SPC%WXL, SPC%WXU, NODATA, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
  CLOSE ( LUNTMP ) 
!
  IF ( NODATA ) THEN
    CALL WRTLOG ( 'W-XSCFIL: File ignored - No data in required spectral range')
    RETURN
  END IF
!
  CALL SAVXSC ( IGAS )
!
900 CONTINUE
  FAIL = IOS .GT. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-XSCFIL: Error reading X/S file. IOSTAT=', IOS
!
END SUBROUTINE XSCFIL
END MODULE XSCFIL_SUB
