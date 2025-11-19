MODULE OPNOUT_SUB
CONTAINS
SUBROUTINE OPNOUT ( LUN, NAMTMP, FAIL, ERRMSG, IGAS, IJAC, ILEV, ISPC, ITAN )
!
! VERSION
!   27MAR23 AD Checked.
!   08NOV17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Open spectral output files for current spectral range
!   Called by SPCTAB, SPCWRT, WRTSTT for each spectral range
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE NAMCOM_DAT, ONLY: DIRNAM, LENNAM ! output dir & Max length of filename
!
! SUBROUTINES
    USE MAKNAM_SUB ! Construct filename for RFM output files
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),           INTENT(IN)  :: LUN    ! Next available LUN
    CHARACTER(*),          INTENT(IN)  :: NAMTMP ! Filename template
    LOGICAL,               INTENT(OUT) :: FAIL   ! TRUE if a fatal error occurs
    CHARACTER(80),         INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
    INTEGER(I4), OPTIONAL, INTENT(IN)  :: IGAS   ! Absorber index
    INTEGER(I4), OPTIONAL, INTENT(IN)  :: IJAC   ! Jacobian index 
    INTEGER(I4), OPTIONAL, INTENT(IN)  :: ILEV   ! Level index 
    INTEGER(I4), OPTIONAL, INTENT(IN)  :: ISPC   ! Spectral range index 
    INTEGER(I4), OPTIONAL, INTENT(IN)  :: ITAN   ! Tangent path index
!
! LOCAL VARIABLES
    INTEGER(I4)         :: IOS    ! Value of IOSTAT on OPEN
    CHARACTER(11)       :: FORSTR ! Value for FORM keyword in OPEN statement
    CHARACTER(2*LENNAM) :: NAMOUT ! Name of file actually opened, incl.RUNID
    CHARACTER(7)        :: STASTR ! Value for STATUS keyword in OPEN statement
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( PRESENT ( IGAS ) ) THEN            ! .tab file
    CALL MAKNAM ( NAMTMP, NAMOUT, IGAS=IGAS, ISPC=ISPC )
  ELSE                                    ! spectral file
    CALL MAKNAM ( NAMTMP, NAMOUT, ISPC=ISPC, ITAN=ITAN, IJAC=IJAC, ILEV=ILEV )
  END IF
!
  NAMOUT = TRIM ( DIRNAM ) // TRIM ( NAMOUT ) 
  CALL WRTLOG ( 'I-OPNOUT: Opening output file: ' // NAMOUT ) 
!
! Open file with appropriate STATUS, FORM according to NEWFLG, BINFLG
  IF ( NEWFLG ) THEN
    STASTR = 'NEW' 
  ELSE 
    STASTR = 'UNKNOWN'
  END IF
  IF ( BINFLG ) THEN
    FORSTR = 'UNFORMATTED' 
  ELSE 
    FORSTR = 'FORMATTED'
  END IF
  OPEN ( UNIT=LUN, FILE=NAMOUT, STATUS=STASTR, FORM=FORSTR, ACTION='WRITE', &
         IOSTAT=IOS )
!
  FAIL = IOS .NE. 0  
  IF ( FAIL ) WRITE ( ERRMSG, * ) 'F-OPNOUT: Error opening file. IOSTAT=', IOS
!
END SUBROUTINE OPNOUT
END MODULE OPNOUT_SUB

