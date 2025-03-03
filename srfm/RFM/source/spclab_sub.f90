MODULE SPCLAB_SUB
CONTAINS
SUBROUTINE SPCLAB ( LABEL, FAIL, ERRMSG )
!
! VERSION
!   17NOV23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Add Spc.range label to SPCCOM
!   Called by DRVSPC for each spectral range
!   Also increases NSPC
!
! VARIABLE KINDS
    USE KIND_DAT  
!
! GLOBAL DATA
    USE SPCCOM_DAT ! Spectral range data
!
  IMPLICIT NONE 
!
! ARGUMENTS 
    CHARACTER(*),  INTENT(IN)  :: LABEL  ! Spectral range label.
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL :: ONESPC = .FALSE.            ! T=only one spc.range allowed
    TYPE(SPCTYP), ALLOCATABLE :: SPCSAV(:) ! Copy of SPC during reallocation
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  FAIL = .FALSE.
!
! If unlabelled range, only permitted if this is the first and only range
  IF ( NSPC .EQ. 0 ) THEN
    ONESPC = LABEL .EQ. ''
  ELSE IF ( ONESPC .OR. LABEL .EQ. '' ) THEN  
    FAIL = .TRUE.
    ERRMSG = 'F-SPCLAB: Unlabelled spectral range only allowed ' // &
             'if it is the *only* range'
    RETURN
  END IF
!
! Check Label within permitted length (LENSPC, in SPCCOM_DAT)
  IF ( LEN ( LABEL ) .GT. LENSPC ) THEN
    FAIL = .TRUE.
    WRITE ( ERRMSG, '(A,I2)' ) 'F-SPCLAB: Spc.range Label: ' // &
      TRIM ( LABEL(1:24) ) // ' exceeds max.length LENSPC=', LENSPC
    RETURN
  END IF
!  
! Find if label already listed 
  IF ( NSPC .GT. 0 ) THEN
    IF ( ANY ( SPC%LAB .EQ. LABEL ) ) THEN 
      FAIL = .TRUE.
      ERRMSG = 'F-SPCLAB: Repeated spectral range label: ' // LABEL
      RETURN
    END IF
  END IF
!
  IF ( ALLOCATED ( SPC ) ) CALL MOVE_ALLOC ( SPC, SPCSAV )
  NSPC = NSPC + 1
  ALLOCATE ( SPC(NSPC) ) 
  IF ( ALLOCATED ( SPCSAV ) ) SPC(1:NSPC-1) = SPCSAV
!
  SPC(NSPC)%LAB = LABEL
  SPC(NSPC)%IGD = 0
!
END SUBROUTINE SPCLAB
END MODULE SPCLAB_SUB

