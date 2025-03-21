MODULE IDXSET_FNC
CONTAINS
INTEGER(I4) FUNCTION IDXSET ( GOTPAR )
!
! VERSION
!   16AUG24 AD Checked.
!   11AUG23 AD Original
!
! DESCRIPTION
!   Return index of line parameter set
!   General purpose module.
!   Creates a new set if not already listed
!   If called without GOTPAR argument, creates/returns set assuming
!   mandatory line parameters only.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL VARIABLES
    USE SETCOM_DAT ! HITRAN line parameter sets
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL, OPTIONAL, INTENT(IN) :: GOTPAR(:) ! T=parameter expected for file
!
! LOCAL VARIABLES
    INTEGER(I4) :: ISET           ! Counter for datasets
    LOGICAL     :: FILPAR(MAXPAR) ! Flags for new dataset 
    LOGICAL, ALLOCATABLE :: SETSAV(:,:) ! Saved copy during reallocation
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  IF ( PRESENT ( GOTPAR ) ) THEN
    FILPAR = GOTPAR
  ELSE
    FILPAR(1:NRQPAR) = .TRUE.
    FILPAR(NRQPAR+1:) = .FALSE.
  END IF
!
  DO ISET = 1, NSET
    IF ( ALL ( FILPAR .EQV. SET(:,ISET) ) ) THEN
      IDXSET = ISET       ! Found match with existing dataset
      RETURN
    END IF
  END DO
!
! Add new dataset
  NSET = NSET + 1
  IF ( NSET .GT. 1 ) CALL MOVE_ALLOC ( SET, SETSAV ) 
  ALLOCATE ( SET(MAXPAR,NSET) ) 
  IF ( NSET .GT. 1 ) SET(:,1:NSET-1) = SETSAV
  SET(:,NSET) = FILPAR
!
  IDXSET = NSET
!
END FUNCTION IDXSET
END MODULE IDXSET_FNC


