MODULE LUTFIL_SUB
CONTAINS
SUBROUTINE LUTFIL ( NAMLUT, FAIL, ERRMSG )
!
! VERSION
!   17JUN24 AD Checked.
!   20DEC17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Check if LUT file is required
!   Called by DRVLUT for each file specified in *LUT section of driver table,
!   and by LUTDEF if filename template used.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LFLCOM_DAT ! Look-Up Table files
    USE SPCCOM_DAT ! Spectral range data
    USE RFMLUN_DAT, ONLY: LUNTMP, LUNNXT ! Logical Unit Numbers of RFM files
!
! SUBROUTINES
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
    USE LUTINF_SUB ! Extract info from TAB-LUT header record
    USE NAMGAS_FNC ! Return molecule name + (iso) associated with GASCOM index
    USE REAQAL_SUB ! Decode Qualifier strings '(IQAL)' or '(IQAL:JQAL)'
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMLUT ! Name of LUT file (+qualifier)
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL       :: ANYQAL ! T=qualifiers appended to end of file name
    LOGICAL       :: BINFIL ! T=binary file, F=ASCII file
    LOGICAL       :: TWOQAL ! T=pair of values in qualifier string, F=one value
    LOGICAL       :: USEFIL ! T=file contains usable data
    INTEGER(I4)   :: IDXISO ! Isotope of tabulated data (0=no differentation)
    INTEGER(I4)   :: IDXMOL ! HITRAN/RFM index of tabulated data
    INTEGER(I4)   :: IGAS   ! Index of molec,iso  in GAS 
    INTEGER(I4)   :: IPT    ! Pointer to '(' marking start of qualifier data
    INTEGER(I4)   :: IQAL   ! First value in QAL string
    INTEGER(I4)   :: ISPC   ! Counter for RFM spectral ranges
    INTEGER(I4)   :: JQAL   ! Second value in QAL string (if more than 1)
    INTEGER(I4)   :: LPT    ! Length of filename part of NAMLUT
    REAL(R8)      :: V1     ! Lower wavenumber limit of LUT [/cm]
    REAL(R8)      :: V2     ! Upper wavenumber limit of LUT [/cm]
    REAL(R8)      :: WNOTOL ! Tolerance for matching wavenumber [cm-1]
    CHARACTER(14) :: QALSTR ! Qualifier string (if any) appended to filename
    TYPE(LFLTYP), ALLOCATABLE :: LFLSAV(:) ! Saved LFL during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL WRTLOG ( 'I-LUTFIL: Checking LUT file:' // NAMLUT )
!
! Separate any qualifier string from the filename
  LPT = LEN_TRIM ( NAMLUT ) 
  ANYQAL = INDEX ( NAMLUT, '(' ) .GT. 1
  IF ( ANYQAL ) THEN
    IPT = INDEX ( NAMLUT, '(' )
    QALSTR = NAMLUT(IPT:LPT)
    LPT = IPT - 1
  END IF
!
! Identify LUT file type, molecule and spectral range
  CALL LUTINF ( LUNTMP, NAMLUT(1:LPT), BINFIL, IDXMOL, IDXISO, V1, V2, &
                FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
  IGAS = IDXGAS ( IDXMOL, IDXISO ) 
  IF ( IGAS .EQ. 0 ) THEN
    CALL WRTLOG ( 'W-LUTFIL: file ignored - not required species' )
    RETURN
  END IF
!
  USEFIL = .FALSE.
  DO ISPC = 1, NSPC
    WNOTOL = SPC(ISPC)%WNR * 0.1
    IF ( V1 .LE. SPC(ISPC)%WNL+WNOTOL .AND. &
         V2 .GE. SPC(ISPC)%WNU-WNOTOL           ) THEN
      IF ( IDXLFL(IGAS,ISPC) .NE. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-LUTFIL: ' // TRIM ( NAMGAS(IGAS) ) // ' has 2 possible' &
                 // ' LUT files for spectral range ' // SPC(ISPC)%LAB
        RETURN
      END IF
      IF ( .NOT. USEFIL ) THEN
        USEFIL = .TRUE.
        IF ( ALLOCATED ( LFL ) ) CALL MOVE_ALLOC ( LFL, LFLSAV )
        NLFL = NLFL + 1
        ALLOCATE ( LFL(NLFL) ) 
        IF ( ALLOCATED ( LFLSAV ) ) LFL(1:NLFL-1) = LFLSAV
        LFL(NLFL)%NAM = NAMLUT(1:LPT)
        LFL(NLFL)%BIN = BINFIL
        LFL(NLFL)%IGS = IGAS
        LFL(NLFL)%NDP = 1
        LFL(NLFL)%NDT = 1
        LFL(NLFL)%LUN = LUNNXT   ! Separate LUNs assigned to each TAB-LUT file
        LUNNXT = LUNNXT + 1
      END IF
      IDXLFL(IGAS,ISPC) = NLFL
    END IF
  END DO
  IF ( .NOT. USEFIL ) THEN
    CALL WRTLOG ( 'W-LUTFIL: file ignored - does not match any spectral range' )
    RETURN
  END IF
!
! Check and save any qualifiers limiting contents
  IF ( ANYQAL ) THEN
    CALL REAQAL ( QALSTR, IQAL, FAIL, ERRMSG, TWOQAL, JQAL ) 
    IF ( FAIL ) RETURN
    IF ( .NOT. TWOQAL ) THEN
      ERRMSG = 'F-LUTQAL: Two qualifiers expected for TAB LUT'
      FAIL = .TRUE.
      RETURN
    END IF 
    LFL(NLFL)%NDP = IQAL
    LFL(NLFL)%NDT = JQAL
  END IF
!
END SUBROUTINE LUTFIL
END MODULE LUTFIL_SUB
