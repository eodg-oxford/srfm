MODULE ATMPSI_SUB
CONTAINS
SUBROUTINE ATMPSI ( NAMATM, FAIL, ERRMSG )
!
! VERSION
!   16AUG19 AD Allow for multiple qualifiers. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Extract Psi angle from brackets following .atm filename
!   Called by DRVATM, DRVNTE (also for .nte filename)
!   If called with GRA flag this returns the index IGRA for the Psi angle.
!   If called without GRA flag this returns IGRA=0
!   Also strips off '(...)' from NAMATM leaving just the filename.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT, ONLY: GRAFLG ! T = use horizontal gradients
    USE GRACOM_DAT, ONLY: PSIATM ! Horiz angle for next profiles
!
! SUBROUTINES
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(INOUT) :: NAMATM ! Name of file 
    LOGICAL,       INTENT(OUT)   :: FAIL   ! Set TRUE if fatal error detected
    CHARACTER(80), INTENT(OUT)   :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: PSIMAX = 90.0 ! Max horizontal angle [deg]
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IEND   ! Location of end of text field in RECORD 
    INTEGER(I4)   :: ILEN   ! Length of qualifier substring
    INTEGER(I4)   :: IOS    ! Saved value of IOSTAT for error message
    INTEGER(I4)   :: ISTA   ! Location of start of text field in RECORD 
    CHARACTER(80) :: MESSGE ! Text message for LOG file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
!
! Search for each pair of '(...)' in filename
  IEND = 0
  ISTA = INDEX ( NAMATM, '(' ) 
  DO WHILE ( ISTA .GT. IEND )                ! found another '('
    ILEN = INDEX ( NAMATM(ISTA:), ')' )   ! relative index of next ')'
    IF ( ILEN .LE. 2 ) THEN               ! ILEN = 2 means '()', ie no contents
      FAIL = .TRUE.
      ERRMSG = 'F-ATMPSI: Failed to identify (.) qualifier in filename'
      RETURN
    END IF
    IEND = ISTA + ILEN - 1                ! absolute index of next ')'
!
! Try reading a number from '(...)' contents
    READ ( NAMATM(ISTA+1:IEND-1), *, IOSTAT=IOS ) PSIATM
! 
    IF ( IOS .EQ. 0 ) THEN      ! A readable real number
! At this point assume that the psi angle qualifier has been identified
      IF ( GRAFLG ) THEN
        IF ( ABS(PSIATM) .LE. PSIMAX ) THEN
          WRITE ( MESSGE, '(A,F7.3)' ) &
            'I-ATMPSI: Loading horizontal gradient profiles at PSI=', PSIATM
          CALL WRTLOG ( MESSGE ) 
        ELSE
          WRITE ( ERRMSG, '(A,F7.3,A)' ) &
          'F-ATMPSI: |Psi angle| > allowed maximum magnitude=', PSIMAX, ' [deg]'
          FAIL = .TRUE.
          RETURN
        END IF
      ELSE
        CALL WRTLOG ( 'W-ATMPSI: GRA flag disabled ' &
                    // 'so profile Horiz.angles ignored.' )
      END IF
      NAMATM = NAMATM(1:ISTA-1) // NAMATM(IEND+1:)
      RETURN               ! Normal exit with Psi angle read and removed
    END IF
    ISTA = IEND + INDEX ( NAMATM(IEND+1:), '(' )  ! index of next '(', or 0.
  END DO
!      
! If it gets this far, no Psi angle has been successfully read
  IF ( GRAFLG ) THEN
    CALL WRTLOG ( 'W-ATMPSI: PSI angle not specified - assuming PSI=0' )
    PSIATM = 0.0
  END IF      
!
END SUBROUTINE ATMPSI
END MODULE ATMPSI_SUB
