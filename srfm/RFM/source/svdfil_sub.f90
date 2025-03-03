MODULE SVDFIL_SUB
CONTAINS
SUBROUTINE SVDFIL ( NAMSVD, FAIL, ERRMSG )
!
! VERSION
!   19JUL24 AD Checked.
!   26JAN18 AD Adapted from LUTFIL. Checked.
!
! DESCRIPTION
!   Check if SVD file is required
!   Called by DRVSVD for each file specified in *SVD section of driver table,
!   and by SVDDEF if filename template used.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE SFLCOM_DAT ! SVD-compressed LUT files
    USE SPCCOM_DAT ! Spectral range data
!
! SUBROUTINES
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
    USE SVDINF_SUB ! Extract info from SVD-LUT header
    USE NAMGAS_FNC ! Return molecule name + (iso) associated with GASCOM index
    USE REAQAL_SUB ! Decode Qualifier strings '(IQAL)' or '(IQAL:JQAL)'
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMSVD ! Name of SVD-LUT file (+qualifier)
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL       :: ANYQAL ! T=qualifiers appended to end of file name
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
    CHARACTER(3)  :: TAB    ! Tabulated function
    TYPE(SFLTYP), ALLOCATABLE :: SFLSAV(:) ! Saved SFL during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL WRTLOG ( 'I-SVDFIL: Checking SVD file:' // NAMSVD )
!
! Separate any qualifier string from the filename
  LPT = LEN_TRIM ( NAMSVD ) 
  ANYQAL = INDEX ( NAMSVD, '(' ) .GT. 1
  IF ( ANYQAL ) THEN
    IPT = INDEX ( NAMSVD, '(' )
    QALSTR = NAMSVD(IPT:LPT)
    LPT = IPT - 1
  END IF
!
! Identify SVD molecule and spectral range
  CALL SVDINF ( NAMSVD(1:LPT), TAB, IDXMOL, IDXISO, V1, V2, &
                FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
  IGAS = IDXGAS ( IDXMOL, IDXISO ) 
  IF ( IGAS .EQ. 0 ) THEN
    CALL WRTLOG ( 'W-SVDFIL: file ignored - not required species' )
    RETURN
  END IF
!
  USEFIL = .FALSE.
  DO ISPC = 1, NSPC
    WNOTOL = SPC(ISPC)%WNR * 0.1
    IF ( V1 .LE. SPC(ISPC)%WNL+WNOTOL .AND. &
         V2 .GE. SPC(ISPC)%WNU-WNOTOL           ) THEN
      IF ( IDXSFL(IGAS,ISPC) .NE. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-SVDFIL: ' // TRIM ( NAMGAS(IGAS) ) // ' has 2 possible' &
                 // ' SVD-LUT files for spectral range ' // SPC(ISPC)%LAB
        RETURN
      END IF
      IF ( .NOT. USEFIL ) THEN
        USEFIL = .TRUE.
        IF ( ALLOCATED ( SFL ) ) CALL MOVE_ALLOC ( SFL, SFLSAV )
        NSFL = NSFL + 1
        ALLOCATE ( SFL(NSFL) ) 
        IF ( ALLOCATED ( SFLSAV ) ) SFL(1:NSFL-1) = SFLSAV
        SFL(NSFL)%NAM = NAMSVD(1:LPT)
        SFL(NSFL)%IGS = IGAS
        SFL(NSFL)%NL = 0 
      END IF
      IDXSFL(IGAS,ISPC) = NSFL
    END IF
  END DO
  IF ( .NOT. USEFIL ) THEN
    CALL WRTLOG ( 'W-SVDFIL: file ignored - does not match any spectral range' )
    RETURN
  END IF
!
! Check and save any qualifiers limiting contents
  IF ( ANYQAL ) THEN
    CALL REAQAL ( QALSTR, IQAL, FAIL, ERRMSG, TWOQAL, JQAL ) 
    IF ( FAIL ) RETURN
    IF ( TWOQAL ) THEN
      ERRMSG = 'F-LUTQAL: Only one qualifier expected for SVD-LUT'
      FAIL = .TRUE.
      RETURN
    END IF
    SFL(NSFL)%NL = IQAL
  END IF
!
END SUBROUTINE SVDFIL
END MODULE SVDFIL_SUB
