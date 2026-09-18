MODULE ADDATM_SUB
CONTAINS
SUBROUTINE ADDATM ( LEV, USEHGT, IDXATM ) 
!
! VERSION
!   01MAY26 AD Bug#55 check for IDXATM .GT. NATM for logical error.
!   01AUG25 AD Simplified. Also update LNPATM here if allocated.
!   26JUL24 AD Checked.
!   24JUN19 AD Remove ATMAUX. Checked.
!   21JUN17 AD Add USEHGT argument to allow for interp in lnp as well as hgt.
!              Avoid interpolation except at added level.
!   01JUN17 AD Original. Checked.
!
! DESCRIPTION
!   Add extra level to atm profiles in ATMCOM
!   Called by ATMLEV.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE GRACOM_DAT ! Atmospheric 2-D field
    USE JACCOM_DAT ! Jacobian data
    USE TANCOM_DAT ! Tangent path data
    USE OBSCOM_DAT, ONLY: IATOBS ! Index of atm level of observer
    USE QFNCOM_DAT, ONLY: NQFN   ! No. Vib.Partition functions
!
! SUBROUTINES
    USE IBRAKT_GEN ! Lower index of array interpolation
    USE VAL1DI_GEN ! Interpolate value from 1D array
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4),    INTENT(IN)  :: LEV    ! Alt/Press of level to be inserted
    LOGICAL,     INTENT(IN)  :: USEHGT ! T=insert altitude, F=insert pressure
    INTEGER(I4), INTENT(OUT) :: IDXATM ! Index assigned to new level
!
! LOCAL VARIABLES
    INTEGER(I4) :: IATM   ! Index of level below inserted level
    INTEGER(I4) :: IPSI   ! Counter for horizontal locations
    INTEGER(I4) :: IQFN   ! Counter for Vib Part Fnc profiles
    INTEGER(I4) :: IVIB   ! Counter for Vib Tem profiles
    INTEGER(I4) :: IVMR   ! Counter for VMR profiles
    REAL(R4)    :: GRDLEV ! Vertical coordinate of level to be inserted
    INTEGER(I4), ALLOCATABLE :: IDXOLD(:)  ! Indices of orig profile levels 
    REAL(R4), ALLOCATABLE :: GRDORG(:)     ! Original profile altitudes
    REAL(R4), ALLOCATABLE :: R1DSAV(:)     ! Saved value of 1d profile
    REAL(R4), ALLOCATABLE :: R2DSAV(:,:)   ! Saved value of VMRATM, VIBATM
    REAL(R4), ALLOCATABLE :: R3DSAV(:,:,:) ! Saved value of VMRATM, VIBATM
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( USEHGT ) THEN
    GRDLEV = LEV
    IDXATM = IBRAKT ( HGTATM, GRDLEV ) + 1
  ELSE
    GRDLEV = LOG ( LEV ) 
    IDXATM = IBRAKT ( LNPATM, GRDLEV ) + 1
  END IF
  IF ( IDXATM .LE. 1 .OR. IDXATM .GT. NATM ) STOP 'F-ADDATM: Logical error'
!
  ALLOCATE ( IDXOLD(NATM+1) )
! Indices of orignal profile levels in new profiles (missing out IDXATM)
  IDXOLD = (/ ( IATM, IATM=1, IDXATM-1 ), ( IATM, IATM=IDXATM+1, NATM+1 ) /)

  NATM = NATM + 1
  IF ( USEHGT ) THEN
    CALL MOVE_ALLOC ( HGTATM, GRDORG ) 
    ALLOCATE ( HGTATM(NATM) ) 
    HGTATM(IDXOLD) = GRDORG
    HGTATM(IDXATM) = GRDLEV
!
    IF ( ALLOCATED ( LNPATM ) ) THEN
      CALL MOVE_ALLOC ( LNPATM, R1DSAV ) 
      ALLOCATE ( LNPATM(NATM) ) 
      LNPATM(IDXOLD) = R1DSAV
      LNPATM(IDXATM) = VAL1DI ( GRDORG, GRDLEV, R1DSAV, .FALSE. )
    END IF
  ELSE
    CALL MOVE_ALLOC ( LNPATM, GRDORG ) 
    ALLOCATE ( LNPATM(NATM) ) 
    LNPATM(IDXOLD) = GRDORG
    LNPATM(IDXATM) = GRDLEV
!
    CALL MOVE_ALLOC ( HGTATM, R1DSAV ) 
    ALLOCATE ( HGTATM(NATM) ) 
    HGTATM(IDXOLD) = R1DSAV
    HGTATM(IDXATM) = VAL1DI ( GRDORG, GRDLEV, R1DSAV, .FALSE. )
  END IF
!
  CALL MOVE_ALLOC ( PREATM, R1DSAV ) 
  ALLOCATE ( PREATM(NATM) ) 
  PREATM(IDXOLD) = R1DSAV
  IF ( USEHGT ) THEN
    PREATM(IDXATM) =  VAL1DI ( GRDORG, GRDLEV, R1DSAV, .TRUE. ) 
  ELSE
    PREATM(IDXATM) = LEV
  END IF
!
  CALL MOVE_ALLOC ( TEMATM, R1DSAV ) 
  ALLOCATE ( TEMATM(NATM) ) 
  TEMATM(IDXOLD) = R1DSAV
  TEMATM(IDXATM) = VAL1DI ( GRDORG, GRDLEV, R1DSAV, .FALSE. )
!
  IF ( ALLOCATED ( EXTATM ) ) THEN
    CALL MOVE_ALLOC ( EXTATM, R1DSAV ) 
    ALLOCATE ( EXTATM(NATM) ) 
    EXTATM(IDXOLD) = R1DSAV
    EXTATM(IDXATM) = VAL1DI ( GRDORG, GRDLEV, R1DSAV, .FALSE. )
  END IF
!
  CALL MOVE_ALLOC ( VMRATM, R2DSAV ) 
  ALLOCATE ( VMRATM(NATM,NVMR) )
  DO IVMR = 1, NVMR
    VMRATM(IDXOLD,IVMR) = R2DSAV(:,IVMR)
    VMRATM(IDXATM,IVMR) = &
      VAL1DI ( GRDORG, GRDLEV, R2DSAV(:,IVMR), .NOT. LINVMR(IVMR) )
  END DO
!
  IF ( NQFN .GT. 0 ) THEN
    CALL MOVE_ALLOC ( QFNATM, R2DSAV ) 
    ALLOCATE ( QFNATM(NATM,NQFN) )
    DO IQFN = 1, NQFN
      QFNATM(IDXOLD,IQFN) = R2DSAV(:,IQFN)
      QFNATM(IDXATM,IQFN) = &
        VAL1DI ( GRDORG, GRDLEV, R2DSAV(:,IQFN) )
    END DO
  END IF   
!
  IF ( NVIB .GT. 0 ) THEN
    CALL MOVE_ALLOC ( VIBATM, R2DSAV ) 
    ALLOCATE ( VIBATM(NATM,NVIB) )
    DO IVIB = 1, NVIB
      VIBATM(IDXOLD,IVIB) = R2DSAV(:,IVIB)
      VIBATM(IDXATM,IVIB) = &
        VAL1DI ( GRDORG, GRDLEV, R2DSAV(:,IVIB), .FALSE. )
    END DO
  END IF
!
  IF ( NPSI .GT. 0 ) THEN    ! 2D atmosphere
    CALL MOVE_ALLOC ( PREGRA, R2DSAV ) 
    ALLOCATE ( PREGRA(NATM,NPSI) ) 
    DO IPSI = 1, NPSI
      PREGRA(IDXOLD,IPSI) = R2DSAV(:,IPSI)
      PREGRA(IDXATM,IPSI) = &
        VAL1DI ( GRDORG, GRDLEV, R2DSAV(:,IPSI), .TRUE. )
    END DO
!
    CALL MOVE_ALLOC ( TEMGRA, R2DSAV ) 
    ALLOCATE ( TEMGRA(NATM,NPSI) ) 
    DO IPSI = 1, NPSI
      TEMGRA(IDXOLD,IPSI) = R2DSAV(:,IPSI)
      TEMGRA(IDXATM,IPSI) = &
        VAL1DI ( GRDORG, GRDLEV, R2DSAV(:,IPSI), .FALSE. )
    END DO
!
    IF ( ALLOCATED ( EXTGRA ) ) THEN
      CALL MOVE_ALLOC ( EXTGRA, R2DSAV ) 
      ALLOCATE ( EXTGRA(NATM,NPSI) ) 
      DO IPSI = 1, NPSI
        EXTGRA(IDXOLD,IPSI) = R2DSAV(:,IPSI)
        EXTGRA(IDXATM,IPSI) = &
          VAL1DI ( GRDORG, GRDLEV, R2DSAV(:,IPSI), .FALSE. )
      END DO
    END IF
!
    CALL MOVE_ALLOC ( VMRGRA, R3DSAV ) 
    ALLOCATE ( VMRGRA(NATM,NPSI,NVMR) )
    DO IVMR = 1, NVMR
      DO IPSI = 1, NPSI
        VMRGRA(IDXOLD,IPSI,IVMR) = R3DSAV(:,IPSI,IVMR)
        VMRGRA(IDXATM,IPSI,IVMR) = &
           VAL1DI ( GRDORG, GRDLEV, R3DSAV(:,IPSI,IVMR), .NOT. LINVMR(IVMR) )
      END DO
    END DO
!
    IF ( NVIB .GT. 0 ) THEN
      CALL MOVE_ALLOC ( VIBGRA, R3DSAV ) 
      ALLOCATE ( VIBGRA(NATM,NPSI,NVIB) )
      DO IVIB = 1, NVIB
        DO IPSI = 1, NPSI
          VIBGRA(IDXOLD,IPSI,IVIB) = R3DSAV(:,IPSI,IVIB)
          VIBGRA(IDXATM,IPSI,IVIB) = &
            VAL1DI ( GRDORG, GRDLEV, R3DSAV(:,IPSI,IVMR) )
        END DO
      END DO
   END IF
!
  END IF
!
! Jacobian perturbation levels
  IF ( NJAC .GT. 0 ) THEN  
    WHERE ( JAC%ILO .GE. IDXATM ) JAC%ILO = JAC%ILO + 1
    WHERE ( JAC%IAT .GE. IDXATM ) JAC%IAT = JAC%IAT + 1
    WHERE ( JAC%IUP .GE. IDXATM ) JAC%IUP = JAC%IUP + 1
  END IF
!  
  WHERE ( TAN%IAT .GE. IDXATM ) TAN%IAT = TAN%IAT + 1
!
  IF ( IATOBS .GE. IDXATM ) IATOBS = IATOBS + 1
!
END SUBROUTINE ADDATM
END MODULE ADDATM_SUB
