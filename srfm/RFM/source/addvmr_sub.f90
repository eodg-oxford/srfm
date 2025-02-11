MODULE ADDVMR_SUB
CONTAINS
SUBROUTINE ADDVMR ( IVMR ) 
!
! VERSION
!   21NOV23 AD Checked.
!   01MAY17 AD Original. Checked.
!
! DESCRIPTION
!   Add extra VMR profile to ATMCOM
!   Called by ATMPRF, JACISO
!   Assumes GAS(NGAS) contains info on new molecule
!   If argument IVMR is present this copies existing profiles from IVMR into
!   new profiles (eg if separating out an isotope, but keep the same VMR 
!   profile for the isotope)
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE GASCOM_DAT ! Molecule and isotope data
    USE GRACOM_DAT ! Atmospheric 2-D field
    USE FLGCOM_DAT, ONLY:LINFLG !  T = Assume VMR varies linearly with altitude
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), OPTIONAL, INTENT(IN) :: IVMR ! Index of profile to be copied
!
! LOCAL VARIABLES
    LOGICAL :: COPY ! T=copy from existing profile#IVMR
    LOGICAL,  ALLOCATABLE :: L1DSAV(:)     ! Saved value of 1D logical arrays
    LOGICAL,  ALLOCATABLE :: L2DSAV(:,:)   ! Saved value of 2D logical arrays
    REAL(R4), ALLOCATABLE :: R2DSAV(:,:)   ! Saved value of 2D real arrays
    REAL(R4), ALLOCATABLE :: R3DSAV(:,:,:) ! Saved value of 3D real arrays

! EXECUTABLE CODE --------------------------------------------------------------
!
  COPY = PRESENT ( IVMR ) 
!
  NVMR = NGAS
  CALL MOVE_ALLOC ( VMRATM, R2DSAV )
  CALL MOVE_ALLOC ( SETVMR, L1DSAV )
  ALLOCATE ( VMRATM(NATM,NVMR), SETVMR(NVMR) )
  VMRATM(:,1:NVMR-1) = R2DSAV
  IF ( COPY ) VMRATM(:,NVMR) = VMRATM(:,IVMR)
  SETVMR(1:NVMR-1) = L1DSAV
  SETVMR(NVMR) = .FALSE.
  IF ( COPY ) SETVMR(NVMR) = SETVMR(IVMR) 
!
  CALL MOVE_ALLOC ( LINVMR, L1DSAV )
  ALLOCATE ( LINVMR(NVMR) )
  LINVMR(1:NVMR-1) = L1DSAV
  LINVMR(NVMR) = LINFLG
  IF ( COPY ) LINVMR(NVMR) = LINVMR(IVMR)
!
  CALL MOVE_ALLOC ( NTEVMR, L1DSAV )
  ALLOCATE ( NTEVMR(NVMR) )
  NTEVMR(1:NVMR-1) = L1DSAV
  NTEVMR(NVMR) = GAS(NGAS)%NTE
  IF ( COPY ) NTEVMR(NVMR) = NTEVMR(IVMR)
!
  IF ( NPSI .GT. 0 ) THEN
    CALL MOVE_ALLOC ( VMRGRA, R3DSAV )
    ALLOCATE ( VMRGRA(NATM,NPSI,NVMR) )
    VMRGRA(:,:,1:NVMR-1) = R3DSAV 
    IF ( COPY ) VMRGRA(:,:,NVMR) = VMRGRA(:,:,IVMR) 
    CALL MOVE_ALLOC ( FLVGRA, L2DSAV ) 
    ALLOCATE ( FLVGRA(NPSI,NVMR) ) 
    FLVGRA(:,1:NVMR-1) = L2DSAV
    FLVGRA(:,NVMR) = .FALSE.
    IF ( COPY ) FLVGRA(:,NVMR) = FLVGRA(:,IVMR)
  END IF
!
END SUBROUTINE ADDVMR
END MODULE ADDVMR_SUB
