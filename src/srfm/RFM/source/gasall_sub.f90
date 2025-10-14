MODULE GASALL_SUB
CONTAINS
SUBROUTINE GASALL ( IMIN, FAIL, ERRMSG )
!
! VERSION
!   10JUL24 AD Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Add wildcard gases to list of absorbers
!   Called once by DRVGAS if a wildcard symbol '*' is present in *GAS section
!   of the RFM driver table
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE SPCCOM_DAT ! Spectral range data
!
! SUBROUTINES
    USE ADDGAS_SUB ! Add new molecule/isotope to list of absorbers
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE LSTABS_SUB ! Construct ordered list of absorbers for given spec.range
    USE MOLIDX_SUB ! Give molecule name for HITRAN/RFM index, or vice-versa
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS 
    INTEGER(I4),   INTENT(IN)  :: IMIN   ! Limit for including absorbers
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL       :: NEWGAS ! T=new molecule added
    INTEGER(I4)   :: IABS   ! Absorber counter
    INTEGER(I4)   :: IDXMOL ! HITRAN index of each molecule from LSTABS
    INTEGER(I4)   :: ISPC   ! Spectral range counter
    INTEGER(I4)   :: NABS   ! No. of molecules required/found by LSTABS
    CHARACTER(7)  :: CODMOL ! Molecules name or index
    CHARACTER(80) :: LOGMSG ! Information message for log file
    INTEGER(I4), ALLOCATABLE :: MOLLST(:) ! List of molecule IDs
    REAL(R4),    ALLOCATABLE :: OPTLST(:) ! List of "optical strengths"
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  SELECT CASE ( IMIN ) 
  CASE ( :-1 )      
    ERRMSG = 'F-GASALL: -ve qualifier not allowed for wildcard in *GAS section'
    FAIL = .TRUE.
    RETURN
  CASE ( 99: ) 
    ERRMSG = 'F-GASALL: qualifier >99 not allowed for wildcard in *GAS section'
    FAIL = .TRUE.
    RETURN
  CASE ( 0 ) 
    CALL WRTLOG ( 'I-GASALL: Adding extra absorbers ' ) 
  CASE DEFAULT
    LOGMSG = 'I-GASALL: Adding extra absorbers with ''optical strength'' >=' &
             // C11INT(IMIN) 
    CALL WRTLOG ( LOGMSG )
  END SELECT
!
  DO ISPC = 1, NSPC
    IF ( NSPC .EQ. 1 ) THEN
      LOGMSG = 'I-GASALL: Add: '
    ELSE
      LOGMSG = 'I-GASALL: Spec.Range' // SPC(ISPC)%LAB // ' add: '
    END IF
    CALL WRTLOG ( LOGMSG, .TRUE. )
    CALL LSTABS ( SPC(ISPC)%WNL, SPC(ISPC)%WNU, IMIN, NABS, MOLLST, OPTLST )
!
    DO IABS = 1, NABS
      IDXMOL = MOLLST(IABS)
      CALL ADDGAS ( IDXMOL, 0, NEWGAS )
      IF ( NEWGAS ) THEN
        CODMOL = ''
        CALL MOLIDX ( IDXMOL, CODMOL )
        CALL WRTLOG ( CODMOL, .TRUE. )
      END IF
    END DO           ! End loop over potential new absorbers
    CALL WRTLOG ( '', .FALSE. )
  END DO             ! End loop over spectral ranges
!
END SUBROUTINE GASALL
END MODULE GASALL_SUB

