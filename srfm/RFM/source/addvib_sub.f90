MODULE ADDVIB_SUB
CONTAINS
SUBROUTINE ADDVIB
!
! VERSION
!   28DEC23 AD Checked.
!   01MAY17 AD F90 Original. Checked.
!
! DESCRIPTION
!   Add extra Vib Temp profile to ATMCOM
!   Called by ADDNTE
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE GRACOM_DAT ! Atmospheric 2-D field
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    LOGICAL,  ALLOCATABLE :: L2DSAV(:,:)   ! Saved value of 2D logical arrays
    REAL(R4), ALLOCATABLE :: R2DSAV(:,:)   ! Saved value of 2D real arrays
    REAL(R4), ALLOCATABLE :: R3DSAV(:,:,:) ! Saved value of 3D real arrays

! EXECUTABLE CODE --------------------------------------------------------------
!
  IF ( NVIB .GT. 0 ) THEN
    CALL MOVE_ALLOC ( VIBATM, R2DSAV )
    IF ( NPSI .GT. 0 ) THEN
      CALL MOVE_ALLOC ( VIBGRA, R3DSAV )
      CALL MOVE_ALLOC ( FLNGRA, L2DSAV ) 
    END IF
  END IF
!
  NVIB = NVIB + 1
  ALLOCATE ( VIBATM(NATM,NVIB) )
  IF ( NPSI .GT. 0 ) THEN
    ALLOCATE ( VIBGRA(NATM,NPSI,NVIB) )
    ALLOCATE ( FLNGRA(NPSI,NVIB) ) 
  END IF
!
  IF ( NVIB .GT. 1 ) THEN
    VIBATM(:,1:NVIB-1) = R2DSAV
    IF ( NPSI .GT. 0 ) THEN
      VIBGRA(:,:,1:NVIB-1) = R3DSAV 
      FLNGRA(:,1:NVIB-1) = L2DSAV
      FLNGRA(:,NVIB) = .FALSE.
    END IF
  END IF
!
END SUBROUTINE ADDVIB
END MODULE ADDVIB_SUB
