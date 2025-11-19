MODULE GASTQL_SUB
CONTAINS
SUBROUTINE GASTQL ( QALSTR, FAIL, ERRMSG )
!
! VERSION
!   28AUG24 AD Checked.
!   24JUN19 AD Remove CIA flag. Checked.
!   01MAY17 AD F90 conversion of gasctm.for. Checked.
!
! DESCRIPTION
!   Read molecule text qualifier strings
!   Called by GASCHK for any gas with a qualifier string.
!   Checks for any of the following, and removes them from QALSTR.
!   (ctm) - set molecule for continuum-only component
!         Requires 1) CTM flag enabled, and
!                  2) Molecule is one of CO2,H2O,N2 or O2, and
!                  3) No other qualifiers attached to molecule name
!         Set GAS%SHP = SHPCTM for this molecule 
!   (cia) - set molecule for CIA-only component
!         Requires 1) No other qualifiers attached to molecule name
!         Set GAS%SHP = SHPCIA for this molecule
!   (noctm) - switch off continuum for this molecule
!          Warn if a) CTM flag not enabled, or
!                  b) Molecule is one of CO2,H2O,N2 or O2
!          For H2O: sets SUBH2O = TABFLG
!   (sub) - subtract value at 25cm-1 from line centre from H2O lines
!         Requires 1) Molecule is H2O
!                  2) (nosub) qualifier not also present 
!         Warn if  a) SUBH2O already set by default (from CTM or TAB)
!   (nosub) - normal treatment of H2O lines (ie no subtraction of 25cm-1 value)
!         Requires 1) Molecule is H2O
!         Warn if  a) SUBH2O already set F by default (no CTM or TAB)
!
! VARIABLE KINDS
    USE KIND_DAT   ! Variable kinds
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE SHPCON_DAT ! Line-shape indices
    USE FLGCOM_DAT, ONLY: CTMFLG, TABFLG ! Option flags
    USE IDXCON_DAT, ONLY: IDXH2O ! RFM/HITRAN index for H2O
!
! SUBROUTINES
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(INOUT) :: QALSTR ! String containing qualifiers
    LOGICAL,       INTENT(OUT)   :: FAIL   ! T = fatal error is detected    
    CHARACTER(80), INTENT(OUT)   :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)       :: IPT    ! Pointer to '(' element of each qualifier 
    CHARACTER(LENGAS) :: MOLEC  ! Molecule name
    CHARACTER(80)     :: LOGMSG ! Message sent to RFM log file
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  FAIL = .FALSE.
!
  MOLEC = GAS(NGAS)%COD
!
! First check for presence of sole (ctm) qualifier in QALSTR
  IF ( INDEX ( QALSTR, '(ctm)' ) .GT. 0 ) THEN
    FAIL = .TRUE.
    IF ( .NOT. CTMFLG ) THEN
      ERRMSG = 'F-GASTQL: (ctm) qualifier set for gas=' // &
               MOLEC // ' but CTM flag not enabled'
    ELSE IF ( .NOT. GAS(NGAS)%CTM ) THEN
      ERRMSG = 'F-GASTQL: No continuum data available for gas=' // MOLEC
    ELSE IF ( QALSTR .NE. '(ctm)' ) THEN
      ERRMSG = 'F-GASTQL: No additional qualifiers allowed with (CTM)'
    ELSE 
      FAIL = .FALSE.
      GAS(NGAS)%SHP = SHPCTM
      LOGMSG = 'I-GASTQL: setting gas=' // MOLEC // &
               ' to continuum-only calculation'
      CALL WRTLOG ( LOGMSG )
      QALSTR = ''
    END IF
    RETURN                 ! Exit with '(CTM)' qualifier
  END IF
!
! Then check for presence of sole (cia) qualifier in QALSTR
  IF ( INDEX ( QALSTR, '(cia)' ) .GT. 0 ) THEN
    FAIL = .TRUE.
    IF ( QALSTR .NE. '(cia)' ) THEN
      ERRMSG = 'F-GASTQL: No additional qualifiers allowed with (CIA)'
    ELSE 
      FAIL = .FALSE.
      GAS(NGAS)%SHP = SHPCIA
      LOGMSG = 'I-GASTQL: setting gas=' // MOLEC // ' to CIA-only calculation'
      CALL WRTLOG ( LOGMSG )
      QALSTR = ''
    END IF
    RETURN                     ! Exit with (cia) only qualifier
  END IF
!
! Other qualifiers may exist in combination and with (ISO)(LGQ:UGQ) qualifiers
!
! Check for (noctm) qualifier
  IPT = INDEX ( QALSTR, '(noctm)' ) 
  IF ( IPT .NE. 0 ) THEN
    IF ( .NOT. CTMFLG ) THEN
      LOGMSG = 'W-GASTQL: CTM Flag not set so (noctm) qualifier for gas=' // &
               MOLEC // ' is redundant'
    ELSE IF ( .NOT. GAS(NGAS)%CTM ) THEN        
      LOGMSG = 'W-GASTQL: No continuum data for gas=' // &
               MOLEC // ' so (noctm) qualifier is redundant'
    ELSE
      GAS(NGAS)%CTM = .FALSE.
      LOGMSG = 'I-GASTQL: Excluding continuum for gas=' // MOLEC
    END IF
    CALL WRTLOG ( LOGMSG )
! If no continuum for H2O, subtraction of 25cm-1 offset depends only on TAB flag
    IF ( GAS(NGAS)%IDM .EQ. IDXH2O ) SUBH2O = TABFLG
    QALSTR = QALSTR(1:IPT-1) // QALSTR(IPT+7:) ! remove '(noctm)' from QALSTR
  END IF
!
! Report current state of H2O subtraction
  IF ( GAS(NGAS)%IDM .EQ. IDXH2O ) THEN
    IF ( SUBH2O ) THEN
      LOGMSG = 'I-GASTQL: Default setting: ' // &
               'subtract 25cm-1 offset from H2O lines'
    ELSE
      LOGMSG = 'I-GASTQL: Default setting: '// &
               'do not subtract 25cm-1 offset from H2O lines'
    END IF
    CALL WRTLOG ( LOGMSG )
  END IF
!
! Check for (sub) qualifier
  IPT = INDEX ( QALSTR, '(sub)' )
  IF ( IPT .NE. 0 ) THEN
    FAIL = .TRUE.
    IF ( GAS(NGAS)%IDM .NE. IDXH2O ) THEN
      ERRMSG = 'F-GASTQL: (SUB) qualifier only applicable to H2O'
    ELSE IF ( INDEX ( QALSTR, '(nosub)' ) .GT. 0 ) THEN
      ERRMSG = 'F-GASTQL: both (SUB) and (NOSUB) qualifiers not allowed'
    ELSE
      FAIL = .FALSE.
    END IF
    IF ( FAIL ) RETURN
    QALSTR = QALSTR(1:IPT-1) // QALSTR(IPT+5:)  ! remove '(sub)' from QALSTR
    IF ( SUBH2O ) THEN
      LOGMSG = 'W-GASTQL: (SUB) already set by default for H2O'
    ELSE
      SUBH2O = .TRUE.
      LOGMSG = 'I-GASTQL: (SUB): change to subtract 25cm-1 ' // &
               'offset from H2O lines'
    END IF
    CALL WRTLOG ( LOGMSG )
    RETURN           ! Normal exit with (sub) qualifier
  END IF
!
! Check for (nosub) qualifier
  IPT = INDEX ( QALSTR, '(nosub)' )
  IF ( IPT .NE. 0 ) THEN
    IF ( GAS(NGAS)%IDM .NE. IDXH2O ) THEN
      ERRMSG = 'F-GASTQL: (NOSUB) qualifier only applicable to H2O'
      FAIL = .TRUE.
      RETURN
    END IF
    QALSTR = QALSTR(1:IPT-1) // QALSTR(IPT+7:)
    IF ( .NOT. SUBH2O ) THEN
      LOGMSG = 'W-GASTQL: (NOSUB) already set by default for H2O'
    ELSE
      SUBH2O = .FALSE.
      LOGMSG = 'I-GASTQL: (NOSUB): change to not subtract ' // & 
               '25cm-1 offset from H2O lines'
    END IF
    CALL WRTLOG ( LOGMSG )
  END IF
!
END SUBROUTINE GASTQL
END MODULE GASTQL_SUB

