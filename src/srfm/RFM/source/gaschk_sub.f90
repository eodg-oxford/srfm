MODULE GASCHK_SUB
CONTAINS
SUBROUTINE GASCHK ( GASSTR, NOTGAS, FAIL, ERRMSG )
!
! VERSION
!   03SEP24 AD Checked.
!   01JUL19 AD Bug#22 Reset IDXISO=0 if interpreting hdo or ch3d. Checked.
!   23JUN17 AD Add ADDNEW argument to CHKGAS. Checked.
!   01MAY17 AD F90 conversion
!
! DESCRIPTION
!   Check Gas name and set indices
!   Called by DRVGAS for each field in *GAS section.
!   Note that all strings are converted to lower case for internal use.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! SUBROUTINES
    USE ADDGAS_SUB ! Add new molecule/isotope to list of absorbers
    USE CHKGAS_SUB ! Check for valid molecule name, isotope, Vib.Level
    USE GASQAL_SUB ! Set isotope/band qualifiers for each molecule
    USE GASTQL_SUB ! Read molecule text qualifier strings
    USE LOCASE_FNC ! Convert text string to lower case
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: GASSTR ! Gas name to be tested 
    LOGICAL,       INTENT(OUT) :: NOTGAS ! Set TRUE if gas not recognised
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL     :: ANYQAL ! Set TRUE if any qualifiers added
    LOGICAL     :: NEWGAS ! Set TRUE if new molec/isotope identified
    INTEGER(I4) :: IDUMMY ! Dummy argument for CHKGAS
    INTEGER(I4) :: IDXISO ! Isotope# 
    INTEGER(I4) :: IDXMOL ! HITRAN index of molecule
    INTEGER(I4) :: LGS    ! Length of GAS part of GASSTR without qualifiers
    INTEGER(I4) :: LPT    ! Length of non-blank part of GASSTR
    INTEGER(I4) :: LQS    ! Length of QALSTR written
    CHARACTER(LEN(GASSTR)) :: QALSTR  ! Qualifier data associated with molec.
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL   = .FALSE.
  NOTGAS = .TRUE.
!
  LPT = LEN_TRIM ( GASSTR ) ! Set LPT to last non-blank character in GASSTR
!
! Check for any qualifier data associated with molecule -look for '(' in string
  ANYQAL = INDEX ( GASSTR, '(' ) .GT. 1   ! GASSTR Cannot start with '('
  IF ( ANYQAL ) THEN
    LGS = INDEX ( GASSTR, '(' ) - 1       
    LQS = LPT - LGS
    QALSTR = LOCASE ( GASSTR(LGS+1:LPT) )
  ELSE
    LGS = LPT 
    LQS = 0
  END IF
!
! Check if GASSTR can be decoded as a molecule (+iso)
! T=allow new molecules to be defined
  CALL CHKGAS ( GASSTR(1:LGS), IDXMOL, IDXISO, IDUMMY, FAIL, ERRMSG, .TRUE. )
  IF ( FAIL ) RETURN
  NOTGAS = IDXMOL .EQ. 0 
  IF ( NOTGAS ) RETURN     ! Exit with no molecule identified
!
! Normally CHKGAS will also read any '(...)' and return Isotope#, but in the
! *GAS section of the driver table the qualifiers can be more complicated so
! these are removed earlier. However if gas is a recognised isotopomeric 
! abbreviation ('hdo' or 'ch3d') CHKGAS could still return isotope numbers, 
! in which case pretend these were part of the original qualifier string
  IF ( IDXISO .GT. 0 ) THEN   ! Identified isotopomeric abbreviation
    IF ( ANYQAL ) QALSTR(4:LQS+3) = QALSTR(1:LQS)
    WRITE ( QALSTR(1:3), '(A,I1,A)' ) '(', IDXISO, ')'
    LQS = LQS + 3
    ANYQAL = .TRUE.
    IDXISO = 0   ! Isotope info now back in QALSTR
  END IF
!
  CALL ADDGAS ( IDXMOL, IDXISO, NEWGAS ) 
! Check molecule isn't already loaded
  IF ( .NOT. NEWGAS ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-GASCHK: Repeated Gas=' // GASSTR(1:LGS)
    RETURN
  END IF
!
! Identify and remove any text qualifiers, '(ctm)', '(cia)', '(nosub)' etc
  IF ( ANYQAL ) THEN
    CALL GASTQL ( QALSTR, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    LQS = LEN_TRIM ( QALSTR )  
    ANYQAL = LQS .GT. 0
  END IF
!
! Set isotope/band selection qualifiers
  IF ( ANYQAL ) THEN
    CALL GASQAL ( QALSTR(1:LQS), FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
  END IF
!
  FAIL = .FALSE.                      ! Normal exit with new molecule added 
!
END SUBROUTINE GASCHK
END MODULE GASCHK_SUB
