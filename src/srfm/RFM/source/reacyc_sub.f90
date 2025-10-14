MODULE REACYC_SUB 
CONTAINS
SUBROUTINE REACYC ( WNOMIN, WNOMAX, FAIL, ERRMSG )
!
! VERSION
!   22AUG24 AD Checked.
!   11AUG23 AD Remove BINCOM_DAT.
!   23MAY23 AD Remove WNOHFL. Reposition REAHIT at end of loop.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Read HITRAN data into cyclic buffers
!   Called by SPCFIN at start of each fine mesh calculation.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE HITCOM_DAT ! HITRAN line data
!
! SUBROUTINES
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE REAHIT_SUB ! Read record from HITRAN line data file
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8),      INTENT(IN)  :: WNOMIN ! Min. Wno to store [cm-1]
    REAL(R8),      INTENT(IN)  :: WNOMAX ! Max. Wno to store [cm-1]
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: NCYC1  = 1000   ! Initial size of cyclic buffer
!
! LOCAL VARIABLES
    LOGICAL       :: EOF    ! Set TRUE if end of HITRAN file reached
    INTEGER(I4)   :: ICYC   ! Index of line cyclic buffer
    INTEGER(I4)   :: IERR   ! Error status allocating CYC
    INTEGER(I4)   :: ILIN   ! Counter for lines in Cyclic buffer
    INTEGER(I4)   :: MCYC   ! Old value of NCYC
    CHARACTER(80) :: LOGMSG ! Message sent to log file
    TYPE(HITTYP), ALLOCATABLE :: CYCSAV(:) ! Saved CYC during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! On first call allocate CYC to initial size
  IF ( .NOT. ALLOCATED ( CYC ) ) THEN
    ALLOCATE ( CYC(NCYC1) ) 
    NCYC = NCYC1
  END IF
!
! Count the no. of useful lines already in buffer and label the first 
! ILIN ends up storing 1 + (no.lines which can be dropped)
!
  ILIN = 1
  ICYC = MOD ( ILIN + ICYC1 - 2, NCYC ) + 1  ! Lowest WNO line in buffer
  DO WHILE ( CYC(ICYC)%WNO .LT. WNOMIN .AND. NLIN .GT. 0 )
    NLIN = NLIN - 1
    ILIN = ILIN + 1
    ICYC = MOD ( ILIN + ICYC1 - 2, NCYC ) + 1
  END DO
  ICYC1 = ICYC
!
! Note that to avoid back-stepping and rereading a line it is necessary to load
! one line beyond the required range into the cyclic buffer from each file
!
  HIT%WNO = -1.0D0
  DO WHILE ( HIT%WNO .LT. WNOMAX ) 
    CALL REAHIT ( EOF, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( EOF ) EXIT
    IF ( NLIN .GE. NCYC ) THEN      ! increase size of CYC buffer
      CALL MOVE_ALLOC ( CYC, CYCSAV )
      MCYC = NCYC
      NCYC = NINT ( NCYC * 1.5 ) 
      LOGMSG = 'W-REACYC: NCYC being increased to NCYC=' // C11INT(NCYC)
      CALL WRTLOG ( LOGMSG ) 
      ALLOCATE ( CYC(NCYC), STAT=IERR )
      IF ( IERR .NE. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-REACYC: error resizing CYC buffer'
        RETURN
      END IF
! Reorder useful lines from lowest wavenumber at start of cyclic buffer
      DO ILIN = 1, NLIN
        ICYC = MOD ( ILIN + ICYC1 - 2, MCYC ) + 1  ! original index of line
        CYC(ILIN) = CYCSAV(ICYC)                   ! new index from 1...
      END DO
      ICYC1 = 1             ! new index of lowest wavenumber line
      DEALLOCATE ( CYCSAV ) 
    END IF
! Save new line data from HITCOM
    NLIN = NLIN + 1
    ICYC = MOD ( ICYC1 + NLIN - 2, NCYC ) + 1
    CYC(ICYC) = HIT
  END DO
!
END SUBROUTINE REACYC
END MODULE REACYC_SUB
