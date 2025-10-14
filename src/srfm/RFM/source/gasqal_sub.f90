MODULE GASQAL_SUB
CONTAINS
SUBROUTINE GASQAL ( QALSTR, FAIL, ERRMSG )
!
! VERSION
!   23FEB23 AD Checked.
!   03OCT17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Set isotope/band qualifiers for each molecule 
!   Called by GASCHK for each molecule with appended qualifier string.
!   The layout of qualifiers is, eg for CO2,
!     CO2(3:4)(1:2)(2)(1)(2:3)(3)(1:3)(2:4)(4)(5)
!   which will (a) select bands (3:4), (1:2) for all isotopes
!              (b) select band (2:3) for isotopes 1,2
!              (c) select bands (1:3), (2:3) for isotope 3
!              (d) select all bands for isotopes 4,5
!   Implicitly this is CO2(*)(user-qualifiers)(*:*), ie starts with all
!   isotopes selected, ends with all bands selected. Decoded as pairs of 
!   (isotopes)x(applicable bands)
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
!
! SUBROUTINES
    USE ADDQAL_SUB ! Add to list of molec iso/band qualifiers
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE REAQAL_SUB ! Decode Qualifier strings '(IQAL)' or '(IQAL:JQAL)'
    USE USEQAL_FNC ! Set TRUE if listed in line molecule qualifiers
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(INOUT) :: QALSTR ! String containing qualifiers
    LOGICAL,       INTENT(OUT)   :: FAIL   ! Set TRUE if fatal error is detected
    CHARACTER(80), INTENT(OUT)   :: ERRMSG ! Error message written if FAIL
!
! LOCAL VARIABLES
    LOGICAL       :: LASBND  ! T=last qualifier read was a band (F=an isotope)
    LOGICAL       :: TWOQAL  ! T=pair of qualifiers, F=single qualifier
    INTEGER(I4)   :: IDXMOL  ! HITRAN/RFM index of molecule
    INTEGER(I4)   :: IISO    ! Isotope counter
    INTEGER(I4)   :: ILS     ! Selected lower state Global Quantum number
    INTEGER(I4)   :: IQAL    ! Qualifier counter
    INTEGER(I4)   :: ISO     ! Isotope counter/index
    INTEGER(I4)   :: IUS     ! Selected upper state Global Quantum number
    INTEGER(I4)   :: JISO    ! Secondary isotope counter
    INTEGER(I4)   :: NISMAX  ! Max no.isotopes for this molecule
    INTEGER(I4)   :: NISO    ! No.of isotopes in ISOLST
    INTEGER(I4)   :: ISOLST(MAXISO) ! List of isotopes using current band
    CHARACTER(30) :: ERRSTR  ! Substring for error messages
    CHARACTER(80) :: LOGMSG  ! Message sent to RFM log file
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  IF ( LEN_TRIM ( QALSTR ) .EQ. 0 ) RETURN
!
  LOGMSG = 'I-GASQAL: Gas=' // GAS(NGAS)%COD // &
           ' List of selected isotopes and/or bands ... '
  CALL WRTLOG ( LOGMSG )
!
! Initialise for each new gas
  GAS(NGAS)%QAL = .TRUE.
  IQAL = 0
  LASBND = .TRUE. 
  IDXMOL = GAS(NGAS)%IDM
  NISMAX = GAS(NGAS)%NIS
!
  IF ( NISMAX .EQ. 0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-GASQAL: No isotopes defined for ' // GAS(NGAS)%COD
    RETURN
  END IF
!
! Repeat until all qualifiers decoded from QALSTR
  DO WHILE ( QALSTR(1:1) .EQ. '(' )
    IQAL = IQAL + 1
!
! Set up substring for any error messages applicable to this gas/qualifier
    WRITE ( ERRSTR, '(A,I2,A)' ) &
      'F-GASQAL: Gas='//GAS(NGAS)%COD//' Qual#', IQAL, ':'
!
! Read qualifier contents
    CALL REAQAL ( QALSTR, ISO, FAIL, ERRMSG, TWOQAL, IUS )
    IF ( FAIL ) RETURN
! Single qualifier, interpreted as isotope#
    IF ( .NOT. TWOQAL ) THEN 
!
! Check valid isotope#
      IF ( ISO .EQ. 0 .OR. ABS ( ISO ) .GT. NISMAX ) THEN
        FAIL = .TRUE.
        ERRMSG = ERRSTR // 'bad isotope#=' // TRIM(C11INT(ISO)) // &
                 '  (range=1:' // TRIM(C11INT(NISMAX)) // ')'
        RETURN
      END IF
!
! If starting new isotope list then initialise ISOLST depending on whether
! isotopes are specified as being added or removed
      IF ( LASBND ) THEN  
        LASBND = .FALSE.              ! modifying isotope list
        IF ( ISO .GT. 0 ) THEN        ! specifying isotope to be added
          NISO = 0  ! start new list of isotopes  
        ELSE                          ! specifying isotope to be removed
          ISOLST = (/ ( IISO, IISO = 1, NISMAX ) /)   ! 1,2... etc
          NISO = NISMAX
        END IF
      END IF
!
      IF ( ISO .GT. 0 ) THEN  ! add isotope to list
        CALL WRTLOG ( 'I-GASQAL:   Select isotope#' // C11INT(ISO) ) 
        IF ( ANY ( ISOLST(1:NISO) .EQ. ISO ) ) THEN ! check iso is new in list
          FAIL = .TRUE.
          ERRMSG = ERRSTR // 'repeated isotope#' // C11INT(ISO)
          RETURN
        END IF
        NISO = NISO + 1
        ISOLST(NISO) = ISO
!
      ELSE                ! remove isotope from list
        CALL WRTLOG ( 'I-GASQAL: Deselect isotope#' // C11INT(ISO) )
        JISO = 0
        DO IISO = 1, NISO   ! Copy all remaining isotopes back to ISOLST
          IF ( ISOLST(IISO) .NE. ABS(ISO) ) THEN
            JISO = JISO + 1
            ISOLST(JISO) = ISOLST(IISO) 
          END IF
        END DO
        IF ( JISO .NE. NISO-1 ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I2,A)' ) &
            ERRSTR // 'isotope#', ABS(ISO), ' not previously selected'
          RETURN
        END IF
        NISO = NISO - 1
      END IF  

! Two qualifiers, interpreted as low:high vib.level numbers
    ELSE
      ILS = ISO
      IF ( ILS .LT. 0 .OR. ILS .GT. 999 ) THEN
        FAIL = .TRUE.
        ERRMSG = ERRSTR // 'ILSGQ out of range 0:999, value=' // C11INT(ILS)
        RETURN
      END IF
      IF ( IUS .LT. 0 .OR. IUS .GT. 999 ) THEN
        FAIL = .TRUE.
        ERRMSG = ERRSTR // 'IUSGQ out of range 0:999, value=' // C11INT(IUS)
        RETURN
      END IF
      LASBND = .TRUE.
! If using all isotopes just use 'ISO=0' as wildcard
      IF ( NISO .EQ. NISMAX ) THEN
        ISOLST(1) = 0
        NISO = 1
      END IF
! Check that qualifier not already listed
      DO IISO = 1, NISO
        IF ( USEQAL ( IDXMOL, ISOLST(IISO), ILS, IUS ) ) THEN 
          FAIL = .TRUE.
          ERRMSG = ERRSTR // 'isotope#' // TRIM(C11INT(ISO)) // & 
                   'repeated band (ILS:IUS)'
          RETURN
        END IF
      END DO
      CALL ADDQAL ( IDXMOL, ISOLST(1:NISO), ILS, IUS )
!
      LOGMSG = 'I-GASQAL:     Select lines from vib.band low:high=' // &
               TRIM ( C11INT(ILS) ) // ':' // TRIM ( C11INT(IUS) ) 
       CALL WRTLOG ( LOGMSG )
    END IF
!
  END DO                           
!
! Successfully decoded entire QALSTR 
!
! Check if finished with an isotope list, in which case set ILS,IUS wildcards
  IF ( .NOT. LASBND ) THEN
    IF ( NISO .EQ. NISMAX ) THEN
      ISOLST(1) = 0
      NISO = 1
      CALL WRTLOG ( 'W-GASQAL: all isotopes selected' )
    END IF
    CALL WRTLOG ( 'I-GASQAL:     Select all bands' )
    DO IISO = 1, NISO
      IF ( USEQAL ( IDXMOL, ISOLST(IISO), 0, 0 ) ) THEN
        FAIL = .TRUE.
        ERRMSG = ERRSTR // 'bands already selected for isotope#' // &
                 C11INT( ISOLST(ISO) ) 
        RETURN
      END IF
    END DO
    CALL ADDQAL ( IDXMOL, ISOLST(1:NISO), 0, 0 )
  END IF
!
! Normal exit
  FAIL = .FALSE.
!
END SUBROUTINE GASQAL
END MODULE GASQAL_SUB
