MODULE NTEFIL_SUB
CONTAINS
SUBROUTINE NTEFIL ( NAMNTE, FAIL, ERRMSG )
!
! VERSION
!   24JUN19 AD Use LOG(PREATM) since LNPATM may not be set yet. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Read NTE Vibrational Temperature data from file
!   Called by DRVNTE or NTEDEF for each file specified in *NTE section.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE FLGCOM_DAT, ONLY: GRAFLG ! T=2d atmosphere
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE ADDNTE_SUB ! Add vibrational temperature to NTE data
    USE ADDQFN_SUB ! Add vibrational partition function
    USE IDXNTE_FNC ! Index in NTECOM of encoded Vib.index
    USE INTERP_GEN ! Interpolate array
    USE NXTREC_SUB ! Load next record from input file
    USE OPNFIL_SUB ! Open input file
    USE PRFGRA_SUB ! Copy atm profile from 1D to 2D field
    USE USEMOL_FNC ! T = molecule is required
    USE WRTLOG_SUB ! Write text message to log file
!   
  IMPLICIT NONE 
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMNTE ! Name of .nte file 
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL       :: ENDSEC ! Set TRUE when '*' marker or EOF reached
    LOGICAL       :: NEWNTE ! T=new vib.temp, F=existing vib.temp
    LOGICAL       :: USEFIL ! T=file contains at least one reqd molecule
    INTEGER(I4)   :: IDUMMY ! Dummy integer 
    INTEGER(I4)   :: IDXISO ! HITRAN/RFM Isotopic ID
    INTEGER(I4)   :: IDXLEV ! Vibrational level ID
    INTEGER(I4)   :: IDXMOL ! HITRAN/RFM Molecule ID
    INTEGER(I4)   :: IDXVIB ! Encoded Vib Temp info
    INTEGER(I4)   :: IMOL   ! Counter for different molecules 
    INTEGER(I4)   :: IOS    ! IOSTAT saved for error message
    INTEGER(I4)   :: IQFN   ! Index of loaded Partition Fn profiles
    INTEGER(I4)   :: IQPR   ! Counter for Partition Fn profiles 
    INTEGER(I4)   :: ISET   ! Counter for VT profiles
    INTEGER(I4)   :: IVIB   ! Index of Vib.Temp profile
    INTEGER(I4)   :: NLEV   ! No.of profile levels 
    INTEGER(I4)   :: NMOL   ! No.of different molecules in file
    INTEGER(I4)   :: NQPR   ! Total No.of different partition functions 
    INTEGER(I4)   :: NSET   ! No.of Vib.Temp profiles in current file
    REAL(R4)      :: ENERGY ! Energy of vibrational state
    CHARACTER(80) :: LOGMSG ! Text message sent to log file
    CHARACTER(80) :: RECORD ! Record read from file
    REAL(R4), ALLOCATABLE :: LNPLEV(:)   ! ln(p) levels of profile
    REAL(R4), ALLOCATABLE :: QFNLEV(:)   ! Vib. Part. Fn profile from file
    REAL(R4), ALLOCATABLE :: TKNLEV(:)   ! Kinetic temp [K] profile from file
    REAL(R4), ALLOCATABLE :: VIBLEV(:)   ! Vib. temp [K] profile from file
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  CALL OPNFIL ( LUNTMP, NAMNTE, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
! First three characters of RECORD should be '***'
  READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) RECORD
  IF ( RECORD(1:3) .NE. '***' ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-NTEFIL: first non-comment record does not start with "***"'
    RETURN
  END IF
! IDUMMY is actually MODEL#
  READ ( RECORD(4:), *, IOSTAT=IOS, ERR=900 ) IDUMMY, NSET, NMOL, NLEV, NQPR 
!
! Read in gas, isotope pairs for which there is vib. tem. data
  USEFIL = .FALSE.
  DO IMOL = 1, NMOL 
    READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) IDXMOL, IDUMMY ! IDUMMY=NISO
    IF ( USEMOL(IDXMOL) ) USEFIL = .TRUE.
  END DO
  IF ( .NOT. USEFIL ) THEN
    LOGMSG = 'W-NTEFIL: file not used - no required molecules'
    CALL WRTLOG ( LOGMSG )
    CLOSE ( LUNTMP )
    RETURN
  END IF
!
  IF ( NLEV .LT. 1 ) THEN
    FAIL = .TRUE.
    WRITE ( ERRMSG, * ) &
      'F-NTEFIL: Value of NLEV must be GE 1, actual value=', NLEV
    RETURN
  END IF
  ALLOCATE ( LNPLEV(NLEV) )
!
! Read in pressure profile for this model and convert from p/mb to ln(p/mb)
  READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) LNPLEV
  LNPLEV = LOG ( LNPLEV ) 
!
  ALLOCATE ( TKNLEV(NLEV) )
  READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) ISET     ! should be ISET=0 
  READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) TKNLEV   ! kinetic temp.profile
!
! Read vibrational partition function profiles
  IF ( NQPR .GT. 0 ) THEN
    ALLOCATE ( QFNLEV(NLEV) ) 
    CALL NXTREC ( LUNTMP, RECORD, ENDSEC, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    BACKSPACE ( LUNTMP, IOSTAT=IOS, ERR=900 )
    DO IQPR = 1, NQPR
      READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) IDXMOL, IDXISO
      READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) QFNLEV
      IF ( .NOT. USEMOL(IDXMOL) ) CYCLE
      CALL ADDQFN ( IDXMOL, IDXISO, IQFN )
      QFNATM(:,IQFN) = INTERP ( LNPLEV, LOG(PREATM), QFNLEV ) 
      IF ( GRAFLG ) CALL PRFGRA ( 'QFN', IQFN ) 
    END DO
  END IF
!
  CALL NXTREC ( LUNTMP, RECORD, ENDSEC, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
  BACKSPACE ( LUNTMP, IOSTAT=IOS, ERR=900 )
!
! Read vibrational temperature profiles 
  ALLOCATE ( VIBLEV(NLEV) )
  DO ISET = 1, NSET                             ! IDUMMY = ISET
    READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) &
      IDUMMY, IDXMOL, IDXISO, IDXLEV, ENERGY
    IDXVIB = IDXMOL * 1000000 + IDXISO * 1000 + IDXLEV
    READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) VIBLEV
    IF ( .NOT. USEMOL(IDXMOL) ) CYCLE   ! not a required molecule
    CALL ADDNTE ( IDXVIB, ENERGY, NEWNTE ) 
    IF ( .NOT. NEWNTE ) CYCLE      ! Don't replace loaded Tvib profiles
    VIBLEV = VIBLEV - TKNLEV
    IVIB = IDXNTE ( IDXVIB ) 
    VIBATM(:,IVIB) = INTERP ( LNPLEV, LOG(PREATM), VIBLEV )
    IF ( GRAFLG ) CALL PRFGRA ( 'VIB', IVIB ) 
  END DO
!
  CLOSE ( LUNTMP )
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-NTEFIL: I/O error on NTE data. IOSTAT=', IOS
!
END SUBROUTINE NTEFIL
END MODULE NTEFIL_SUB
