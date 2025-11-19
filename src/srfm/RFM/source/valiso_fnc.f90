MODULE VALISO_FNC 
CONTAINS
LOGICAL FUNCTION VALISO ( IGAS, ISO ) 
!
! VERSION
!   08NOV23 AD Checked.
!   03MAY18 AD Bug#3 SAVE WRNLST array. Checked.
!   01MAY17 AD Original.
!
! DESCRIPTION
!   Check recognised isotope
!   Called by REAHIT
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
!
! SUBROUTINES
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE WRTLOG_SUB ! Write text message to log file
!    
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IGAS ! Index of Gas in GASCOM
    INTEGER(I4), INTENT(IN) :: ISO  ! Isotope#
!
! LOCAL VARIABLES
    LOGICAL       :: WARN      ! T=send warning message to log file
    INTEGER(I4)   :: IDXWRN    ! Encoded IGAS,ISO value
    INTEGER(I4)   :: NWRN = 0  ! No.different warnings so far
    CHARACTER(80) :: LOGMSG    ! Message sent to log file
    INTEGER(I4), SAVE, ALLOCATABLE :: WRNLST(:) ! List of IDXWRN values
    INTEGER(I4), ALLOCATABLE :: WRNSAV(:) ! Saved WRNLST during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  VALISO = ( ISO .LE. GAS(IGAS)%NIS ) 
!
  IF ( .NOT. VALISO ) THEN
    IDXWRN = IGSMOL(IGAS)*1000 + ISO
    IF ( NWRN .EQ. 0 ) THEN
      WARN = .TRUE.
    ELSE 
      WARN = .NOT. ANY ( WRNLST .EQ. IDXWRN ) 
    END IF
    IF ( WARN ) THEN
      LOGMSG = 'W-VALISO: ignoring ' // TRIM ( GAS(IGAS)%COD ) // &
                 ' lines for unlisted Isotope#=' // C11INT(ISO)
      CALL WRTLOG ( LOGMSG )
!
      IF ( ALLOCATED ( WRNLST ) ) CALL MOVE_ALLOC ( WRNLST, WRNSAV )
      NWRN = NWRN + 1
      ALLOCATE ( WRNLST(NWRN ) )
      IF ( ALLOCATED ( WRNSAV ) ) WRNLST(1:NWRN-1) = WRNSAV
      WRNLST(NWRN) = IDXWRN
    END IF
  END IF
!
END FUNCTION VALISO
END MODULE VALISO_FNC
