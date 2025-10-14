MODULE ATMPRF_SUB
CONTAINS
SUBROUTINE ATMPRF ( LABEL, LEVPRF, VALPRF, USEPRF, FAIL, ERRMSG )
!
! VERSION
!   11APR23 AD Checked.
!   19MAR21 AD Add log message if VMR profiles updated.
!   17JAN18 AD Original. Checked.
!
! DESCRIPTION
!   Load profile into ATMCOM
!   Called by ATMFIL for each profile in each .atm file, 
!   and by ATMPAR for single-valued profiles.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Climatology data
    USE FLGCOM_DAT, ONLY: GRAFLG, NTEFLG ! T=horiz.grads, non-LTE calc
    USE IDXCON_DAT, ONLY: IDXAER ! RFM molec# for 'aerosol'
    USE SFCCOM_DAT, ONLY: TEMSFC ! Surface temperature
!
! SUBROUTINES
    USE ADDGAS_SUB ! Add new molecule/isotope to list of absorbers
    USE ADDNTE_SUB ! Add vibrational temperature to NTE data
    USE ADDVMR_SUB ! Add extra VMR profile to ATMCOM
    USE CHKGAS_SUB ! Check for valid molecule name, isotope, Vib.Level
    USE CHKPRF_SUB ! Check atmospheric profile on input
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
    USE IDXNTE_FNC ! Index in NTECOM of encoded Vib.index
    USE INTERP_GEN ! Interpolate array
    USE PRFGRA_SUB ! Copy atm profile from 1D to 2D field
    USE USEMOL_FNC ! T = molecule is required
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: LABEL     ! Profile *LABEL
    REAL(R4),      INTENT(IN)  :: LEVPRF(:) ! Profile grid
    REAL(R4),      INTENT(IN)  :: VALPRF(:) ! Profile values
    LOGICAL,       INTENT(OUT) :: USEPRF    ! T=profile used, F=profile ignored
    LOGICAL,       INTENT(OUT) :: FAIL      ! Set True if a fatal error detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG    ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL     :: NEWGAS ! T=new molecule/isotope added
    INTEGER(I4) :: IDXISO ! Isotope# (or 0 if no isotope specified)
    INTEGER(I4) :: IDXMOL ! RFM/HITRAN ID for molecule
    INTEGER(I4) :: IDXVIB ! Encoded molec/iso/vib.level
    INTEGER(I4) :: IVIB   ! Index of Vib Tem profile
    INTEGER(I4) :: IVMR   ! Index of molec/iso in GASCOM/ATMCOM arrays
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( .NOT. ASSOCIATED ( LEVATM ) ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-ATMPRF: Atmospheric profile levels not yet established'
    RETURN
  END IF
!
  USEPRF = .FALSE.
  SELECT CASE ( LABEL )
  CASE ( '*TEM' )                             ! Found temperature profile
    CALL WRTLOG ( '  '//LABEL ) 
    CALL CHKPRF ( 'TEM', VALPRF, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    TEMATM = INTERP ( LEVPRF, LEVATM, VALPRF, .FALSE., .FALSE. )
    TEMSFC = TEMATM(1)
    IF ( SETTEM ) CALL WRTLOG ( 'W-ATMFIL: Temperature profile replaced' )
    SETTEM = .TRUE.
    IF ( GRAFLG ) CALL PRFGRA ( 'TEM' ) 
  CASE ( '*PRE' )                             ! Found pressure profile
    CALL WRTLOG ( '  '//LABEL ) 
    CALL CHKPRF ( 'PRE', VALPRF, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    PREATM = INTERP ( LEVPRF, LEVATM, VALPRF, .TRUE., .TRUE. )
    IF ( SETPRE ) CALL WRTLOG ( 'W-ATMFIL: Pressure profile replaced' )
    SETPRE = .TRUE.
    IF ( GRAFLG ) CALL PRFGRA ( 'PRE' ) 
  CASE ( '*HGT' )                             ! Found altitude profile
    CALL WRTLOG ( '  '//LABEL ) 
    CALL CHKPRF ( 'HGT', VALPRF, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    HGTATM = INTERP ( LEVPRF, LEVATM, VALPRF, .FALSE., .TRUE. )
    HGTSFC = HGTATM(1)
    HGTTOA = HGTATM(NATM)
    IF ( SETHGT ) CALL WRTLOG ( 'W-ATMFIL: Altitude profile replaced' )
    SETHGT = .TRUE.
    IF ( GRAFLG ) CALL PRFGRA ( 'HGT' ) 
  CASE DEFAULT
    CALL CHKGAS ( LABEL(2:), IDXMOL, IDXISO, IDXVIB )
    IF ( .NOT. USEMOL(IDXMOL) ) RETURN   ! Molecule not recognised or not reqd
    CALL WRTLOG ( '  '//LABEL )
    IF ( IDXVIB .EQ. 0 ) THEN                 ! VMR/aerosol profiles
! atm files won't add new molecules but might add specific isotopes
      CALL ADDGAS ( IDXMOL, IDXISO, NEWGAS )
      IF ( NEWGAS ) CALL ADDVMR             ! Isotope added as a new gas
! load/update VMR profile
      IVMR = IDXGAS ( IDXMOL, IDXISO )
      IF ( IDXMOL .EQ. IDXAER ) THEN
        CALL CHKPRF ( 'EXT', VALPRF, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        EXTATM = INTERP ( LEVPRF, LEVATM, VALPRF )
        IF ( GRAFLG ) CALL PRFGRA ( 'EXT', IVMR ) 
      ELSE
        CALL CHKPRF ( 'VMR', VALPRF, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        VMRATM(:,IVMR) = INTERP ( LEVPRF, LEVATM, VALPRF, .NOT. LINVMR(IVMR) )
        IF ( GRAFLG ) CALL PRFGRA ( 'VMR', IVMR ) 
      END IF
      IF ( SETVMR(IVMR) ) &
        CALL WRTLOG ( 'W-ATMFIL: ' // TRIM(LABEL(2:)) // ' profile replaced' ) 
      SETVMR(IVMR) = .TRUE.
! Vib Temp profile
    ELSE IF ( NTEFLG ) THEN                        
      CALL ADDNTE ( IDXVIB ) 
      IVIB = IDXNTE ( IDXVIB ) 
      VIBATM(:,IVIB) = INTERP ( LEVPRF, LEVATM, VALPRF, .FALSE., .FALSE. )
    END IF
  END SELECT
!
  USEPRF = .TRUE.                 ! Normal exit with profile used
!
END SUBROUTINE ATMPRF
END MODULE ATMPRF_SUB
