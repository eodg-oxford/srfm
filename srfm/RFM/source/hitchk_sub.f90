MODULE HITCHK_SUB
CONTAINS
SUBROUTINE HITCHK 
!
! VERSION
!   12AUG24 AD Checked.
!   11AUG23 AD Rewritten and redefined.
!   30MAY23 AD Renamed from OPNHIT and simplified structure.
!   29JAN20 AD Open all forms of HITRAN data file. Checked.
!   12JUN17 AD Allow for ASCII files. Checked.
!   O1MAY17 AD F90 conversion. 
!
! DESCRIPTION
!   Check molecule assignments for HITRAN files
!   Called once by DRVHIT after all HITRAN files specified.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE HFLCOM_DAT ! HITRAN file data
    USE HITCOM_DAT ! HITRAN line data
!
! SUBROUTINES
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IDM    ! HITRAN index of molecule
    INTEGER(I4)   :: IFLUSE ! Secondary counter for HITRAN files
    INTEGER(I4)   :: IGAS   ! Counter for required molecules
    INTEGER(I4)   :: IHFL   ! Counter for HITRAN files
    INTEGER(I4)   :: NFLUSE ! No. of files used
    CHARACTER(80) :: LOGMSG ! Message sent to log file
    LOGICAL       :: SETIDM(MAXIDM) ! T=file identified for molec#
    TYPE(HFLTYP), ALLOCATABLE :: HFLSAV(:) ! Saved HFL during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  SETIDM = .FALSE.
  NFLUSE = 0
  IFLIDM = 0 
! Go through in reverse order since last file supersedes earlier files
  DO IHFL = NHFL, 1, -1
    WHERE ( SETIDM ) HFL(IHFL)%IDM(:) = .FALSE.
    WHERE ( HFL(IHFL)%IDM ) IFLIDM = IHFL
    IF ( .NOT. ANY ( HFL(IHFL)%IDM ) ) THEN
      WRITE ( LOGMSG, * ) 'W-HITCHK: No molecules used from File#', IHFL
      CALL WRTLOG ( LOGMSG ) 
    ELSE
      NFLUSE = NFLUSE + 1      
    END IF
    SETIDM = SETIDM .OR. HFL(IHFL)%IDM
  END DO
!
! Remove any redundant files from list
  IF ( NFLUSE .LT. NHFL ) THEN 
    IFLIDM = 0
    IFLUSE = 0
    CALL MOVE_ALLOC ( HFL, HFLSAV ) 
    ALLOCATE ( HFL(NFLUSE) ) 
    DO IHFL = 1, NHFL 
      IF ( ANY ( HFLSAV(IHFL)%IDM ) ) THEN
        IFLUSE = IFLUSE + 1
        HFL(IFLUSE) = HFLSAV(IHFL) 
        WHERE ( HFL(IFLUSE)%IDM ) IFLIDM = IFLUSE
      END IF
    END DO
    NHFL = NFLUSE
  END IF
!
! If only one file make sure HIT%IST is set appropriately
! since RFM will read directly into HIT rather than HFL%HIT
  IF ( NHFL .EQ. 1 ) HIT = HFL(1)%HIT
!
! Flag molecules using this file for spectroscopic data 
! (even if no actual lines within spectral range of calculation) 
  DO IGAS = 1, NGAS
    IDM = GAS(IGAS)%IDM
    IF ( IDM .GT. 0 .AND. IDM .LE. MAXIDM ) GAS(IGAS)%HIT = SETIDM(IDM)
  END DO
!
END SUBROUTINE HITCHK
END MODULE HITCHK_SUB
