MODULE ADDGRA_SUB
CONTAINS
SUBROUTINE ADDGRA
!
! VERSION
!   16SEP24 AD Checked.
!   01MAY17 AD F90 conversion of newpsi.for. Checked.
!
! DESCRIPTION
!   Set horiz location in 2D field
!   Called by PRFGRA
!   This module inserts location into GRACOM arrays in order of increasing
!   PSI value.
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
    LOGICAL,  ALLOCATABLE :: L1DSAV(:)     ! Saved 1D logical array 
    LOGICAL,  ALLOCATABLE :: L2DSAV(:,:)   ! Saved 2D logical array 
    REAL(R4), ALLOCATABLE :: R1DSAV(:)     ! Saved 1D real array 
    REAL(R4), ALLOCATABLE :: R2DSAV(:,:)   ! Saved 2D real array 
    REAL(R4), ALLOCATABLE :: R3DSAV(:,:,:) ! Saved 3D real array 
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! On first call add reference profile at PSI=0.0
  IF ( NPSI .EQ. 0 ) THEN
    NPSI = 1
    NOMGRA = 1
    ALLOCATE ( PSIGRA(1) ) ; PSIGRA = 0.0
    ALLOCATE ( FLPGRA(1) ) ; FLPGRA = .FALSE.
    ALLOCATE ( FLTGRA(1) ) ; FLTGRA = .FALSE.
    ALLOCATE ( FLVGRA(1,NVMR) ) ; FLVGRA = .FALSE.
    ALLOCATE ( TEMGRA(NATM,1) ) 
    ALLOCATE ( PREGRA(NATM,1) ) 
    ALLOCATE ( VMRGRA(NATM,1,NVMR) ) 
    ALLOCATE ( DNSGRA(NATM,1) )
    ALLOCATE ( DSHGRA(NATM,1) )
    ALLOCATE ( RFRGRA(NATM,1) )
    IF ( IAXVMR .GT. 0 ) ALLOCATE ( EXTGRA(NATM,1) )
  END IF
!
! Subsequent calls may be for an existing or a new PSI angle
  IF ( MINVAL ( ABS ( PSIATM - PSIGRA ) ) .LE. PSITOL ) THEN
    MPSI = MINLOC ( ABS ( PSIATM - PSIGRA ), 1 ) 
    RETURN   ! existing angle
  END IF
!
  MPSI = COUNT ( PSIGRA .LT. PSIATM - PSITOL ) + 1
  CALL MOVE_ALLOC ( PSIGRA, R1DSAV )
  NPSI = NPSI + 1
  ALLOCATE ( PSIGRA(NPSI) ) 
  PSIGRA = (/ R1DSAV(1:MPSI-1), PSIATM, R1DSAV(MPSI:) /)
!
! Insert flag for new location, marked FALSE.
  CALL MOVE_ALLOC ( FLTGRA, L1DSAV ) 
  ALLOCATE ( FLTGRA(NPSI) )
  FLTGRA = (/ L1DSAV(1:MPSI-1), .FALSE., L1DSAV(MPSI:) /)
  CALL MOVE_ALLOC ( FLPGRA, L1DSAV ) 
  ALLOCATE ( FLPGRA(NPSI) )
  FLPGRA = (/ L1DSAV(1:MPSI-1), .FALSE., L1DSAV(MPSI:) /)
  CALL MOVE_ALLOC ( FLVGRA, L2DSAV ) 
  ALLOCATE ( FLVGRA(NPSI,NVMR) )
  FLVGRA(1:MPSI-1,:) = L2DSAV(1:MPSI-1,:)
  FLVGRA(MPSI,:)     = .FALSE.
  FLVGRA(MPSI+1:,:)  = L2DSAV(MPSI:,:)
!
! Insert profiles for new location
  CALL MOVE_ALLOC ( TEMGRA, R2DSAV ) 
  ALLOCATE ( TEMGRA(NATM,NPSI) )
  TEMGRA(:,1:MPSI-1) = R2DSAV(:,1:MPSI-1)
  TEMGRA(:,MPSI)     = 0.0
  TEMGRA(:,MPSI+1:)  = R2DSAV(:,MPSI:)
!
  CALL MOVE_ALLOC ( PREGRA, R2DSAV ) 
  ALLOCATE ( PREGRA(NATM,NPSI) )
  PREGRA(:,1:MPSI-1) = R2DSAV(:,1:MPSI-1)
  PREGRA(:,MPSI)     = 0.0
  PREGRA(:,MPSI+1:)  = R2DSAV(:,MPSI:)
!
  IF ( IAXVMR .NE. 0 ) THEN
    CALL MOVE_ALLOC ( EXTGRA, R2DSAV ) 
    ALLOCATE ( EXTGRA(NATM,NPSI) )
    EXTGRA(:,1:MPSI-1) = R2DSAV(:,1:MPSI-1)
    EXTGRA(:,MPSI)     = 0.0
    EXTGRA(:,MPSI+1:)  = R2DSAV(:,MPSI:)
  END IF
!
  CALL MOVE_ALLOC ( VMRGRA, R3DSAV ) 
  ALLOCATE ( VMRGRA(NATM,NPSI,NVMR) )
  VMRGRA(:,1:MPSI-1,:) = R3DSAV(:,1:MPSI-1,:)
  VMRGRA(:,MPSI,:)     = 0.0
  VMRGRA(:,MPSI+1:,:)  = R3DSAV(:,MPSI:,:)
!
! The following arrays are filled later by ATMAUX so no need to save contents
  DEALLOCATE ( DNSGRA ) ; ALLOCATE ( DNSGRA(NATM,NPSI) )
  DEALLOCATE ( DSHGRA ) ; ALLOCATE ( DSHGRA(NATM,NPSI) )
  DEALLOCATE ( RFRGRA ) ; ALLOCATE ( RFRGRA(NATM,NPSI) )
!
! If loading a -ve PSI, index of PSI=0 has been incremented
  IF ( PSIATM .LT. 0.0 ) NOMGRA = NOMGRA + 1
!
END SUBROUTINE ADDGRA
END MODULE ADDGRA_SUB
