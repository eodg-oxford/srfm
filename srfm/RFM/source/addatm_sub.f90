MODULE ADDATM_SUB
CONTAINS
SUBROUTINE ADDATM ( LEV, USEHGT, IDXATM ) 
!
! VERSION
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
    REAL(R4), ALLOCATABLE :: GRDORG(:)     ! Original profile altitudes
    REAL(R4), ALLOCATABLE :: R1DSAV(:)     ! Saved value of 1d profile
    REAL(R4), ALLOCATABLE :: R2DSAV(:,:)   ! Saved value of VMRATM, VIBATM
    REAL(R4), ALLOCATABLE :: R3DSAV(:,:,:) ! Saved value of VMRATM, VIBATM
    REAL(R4), POINTER     :: GRDATM(:)     ! Vertical grid of profiles
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( USEHGT ) THEN
    GRDATM => HGTATM
    GRDLEV = LEV
  ELSE
    GRDATM => LNPATM
    GRDLEV = LOG ( LEV ) 
  END IF
  IATM = IBRAKT ( GRDATM, GRDLEV )
  IF ( IATM .LT. 1 .OR. IATM .GT. NATM-1 ) STOP 'F-ADDATM: Logical error'
  IDXATM = IATM + 1
!
  NATM = NATM + 1
  IF ( USEHGT ) THEN
    CALL MOVE_ALLOC ( HGTATM, GRDORG ) 
    ALLOCATE ( HGTATM(NATM) ) 
    GRDATM => HGTATM
  ELSE
    CALL MOVE_ALLOC ( LNPATM, GRDORG ) 
    ALLOCATE ( LNPATM(NATM) ) 
    GRDATM => LNPATM
  END IF
  GRDATM = (/ GRDORG(1:IATM), GRDLEV, GRDORG(IATM+1:) /)
!
  IF ( .NOT. USEHGT ) THEN
    CALL MOVE_ALLOC ( HGTATM, R1DSAV ) 
    ALLOCATE ( HGTATM(NATM) ) 
    HGTATM = (/ R1DSAV(1:IATM), VAL1DI ( GRDORG, GRDLEV, R1DSAV, .TRUE. ), &
                R1DSAV(IATM+1:NATM-1) /)
  END IF
!
  CALL MOVE_ALLOC ( PREATM, R1DSAV ) 
  ALLOCATE ( PREATM(NATM) ) 
  PREATM(1:IATM) = R1DSAV(1:IATM)
  IF ( USEHGT ) THEN
    PREATM(IATM+1) = VAL1DI ( GRDORG, GRDLEV, R1DSAV, .TRUE. ) 
  ELSE
    PREATM(IATM+1) = LEV
  END IF
  PREATM(IATM+2:NATM) = R1DSAV(IATM+1:NATM-1) 
!
  CALL MOVE_ALLOC ( TEMATM, R1DSAV ) 
  ALLOCATE ( TEMATM(NATM) ) 
  TEMATM = (/ R1DSAV(1:IATM), VAL1DI ( GRDORG, GRDLEV, R1DSAV, .FALSE. ), &
              R1DSAV(IATM+1:NATM-1) /)
!
  IF ( ALLOCATED ( EXTATM ) ) THEN
    CALL MOVE_ALLOC ( EXTATM, R1DSAV ) 
    ALLOCATE ( EXTATM(NATM) ) 
    EXTATM = (/ R1DSAV(1:IATM), VAL1DI ( GRDORG, GRDLEV, R1DSAV, .FALSE. ), &
                R1DSAV(IATM+1:NATM-1) /)
  END IF
!
  CALL MOVE_ALLOC ( VMRATM, R2DSAV ) 
  ALLOCATE ( VMRATM(NATM,NVMR) )
  DO IVMR = 1, NVMR
    VMRATM(:,IVMR) = (/ R2DSAV(1:IATM,IVMR), &
          VAL1DI ( GRDORG, GRDLEV, R2DSAV(:,IVMR), .NOT. LINVMR(IVMR) ), &
                        R2DSAV(IATM+1:NATM-1,IVMR) /)
  END DO
!
  IF ( NQFN .GT. 0 ) THEN
    CALL MOVE_ALLOC ( QFNATM, R2DSAV ) 
    ALLOCATE ( QFNATM(NATM,NQFN) )
    DO IQFN = 1, NQFN
      QFNATM(:,IQFN) = (/ R2DSAV(1:IATM,IQFN), &
                          VAL1DI ( GRDORG, GRDLEV, R2DSAV(:,IQFN) ), &
                          R2DSAV(IATM+1:NATM-1,IQFN) /)
    END DO
  END IF   
!
  IF ( NVIB .GT. 0 ) THEN
    CALL MOVE_ALLOC ( VIBATM, R2DSAV ) 
    ALLOCATE ( VIBATM(NATM,NVIB) )
    DO IVIB = 1, NVIB
      VIBATM(:,IVIB) = (/ R2DSAV(1:IATM,IVIB), &
                          VAL1DI ( GRDORG, GRDLEV, R2DSAV(:,IVIB), .FALSE. ), &
                          R2DSAV(IATM+1:NATM-1,IVIB) /)
    END DO
  END IF
!
  IF ( NPSI .GT. 0 ) THEN    ! 2D atmosphere
    CALL MOVE_ALLOC ( PREGRA, R2DSAV ) 
    ALLOCATE ( PREGRA(NATM,NPSI) ) 
    DO IPSI = 1, NPSI
      PREGRA(:,IPSI) = (/ R2DSAV(1:IATM,IPSI), &
                          VAL1DI ( GRDORG, GRDLEV, R2DSAV(:,IPSI), .TRUE. ), &
                          R2DSAV(IATM+1:NATM-1,IPSI) /)
    END DO
!
    CALL MOVE_ALLOC ( TEMGRA, R2DSAV ) 
    ALLOCATE ( TEMGRA(NATM,NPSI) ) 
    DO IPSI = 1, NPSI
      TEMGRA(:,IPSI) = (/ R2DSAV(1:IATM,IPSI), &
                          VAL1DI ( GRDORG, GRDLEV, R2DSAV(:,IPSI), .FALSE. ), &
                          R2DSAV(IATM+1:NATM-1,IPSI) /)
    END DO
!
    IF ( ALLOCATED ( EXTGRA ) ) THEN
      CALL MOVE_ALLOC ( EXTGRA, R2DSAV ) 
      ALLOCATE ( EXTGRA(NATM,NPSI) ) 
      DO IPSI = 1, NPSI
        EXTGRA(:,IPSI) = (/ R2DSAV(1:IATM,IPSI), &
                          VAL1DI ( GRDORG, GRDLEV, R2DSAV(:,IPSI), .FALSE. ), &
                            R2DSAV(IATM+1:NATM-1,IPSI) /)
      END DO
    END IF
!
    CALL MOVE_ALLOC ( VMRGRA, R3DSAV ) 
    ALLOCATE ( VMRGRA(NATM,NPSI,NVMR) )
    DO IVMR = 1, NVMR
      DO IPSI = 1, NPSI
        VMRGRA(:,IPSI,IVMR) = (/ R3DSAV(1:IATM,IPSI,IVMR), &
           VAL1DI ( GRDORG, GRDLEV, R3DSAV(:,IPSI,IVMR), .NOT. LINVMR(IVMR) ), &
                                 R3DSAV(IATM+1:NATM-1,IPSI,IVMR) /)
      END DO
    END DO
!
    IF ( NVIB .GT. 0 ) THEN
      CALL MOVE_ALLOC ( VIBGRA, R3DSAV ) 
      ALLOCATE ( VIBGRA(NATM,NPSI,NVIB) )
      DO IVIB = 1, NVIB
        DO IPSI = 1, NPSI
          VIBGRA(:,IPSI,IVIB) = (/ R3DSAV(1:IATM,IPSI,IVIB), &
                              VAL1DI ( GRDORG, GRDLEV, R3DSAV(:,IPSI,IVMR) ), &
                                   R3DSAV(IATM+1:NATM-1,IPSI,IVIB) /)
        END DO
      END DO
   END IF
!
  END IF
!
! Jacobian perturbation levels
  IF ( NJAC .GT. 0 ) THEN  
    WHERE ( JAC%ILO .GT. IATM ) JAC%ILO = JAC%ILO + 1
    WHERE ( JAC%IAT .GT. IATM ) JAC%IAT = JAC%IAT + 1
    WHERE ( JAC%IUP .GT. IATM ) JAC%IUP = JAC%IUP + 1
  END IF
!  
  WHERE ( TAN%IAT .GT. IATM ) TAN%IAT = TAN%IAT + 1
!
  IF ( IATOBS .GT. IATM ) IATOBS = IATOBS + 1
!
END SUBROUTINE ADDATM
END MODULE ADDATM_SUB
