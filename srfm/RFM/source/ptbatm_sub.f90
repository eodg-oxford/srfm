MODULE PTBATM_SUB
CONTAINS
SUBROUTINE PTBATM ( IJAC )
!
! VERSION
!   20SEP19 AD Bug#23: Call PTBPRF with linear ptb for aerosol. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Perturb/unperturb atmospheric profiles for Jacobian calc
!   Called by JACPTH for each Jacobian element.
!   Routine starts by undoing previous perturbation (except for first call)
!   Last call should be with IJAC = 0, in which case the only action
!   is to undo the previous perturbation.
!   Pertubation sizes (PTB*) are stored in RFMCON_DAT
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE GRACOM_DAT ! Atmospheric 2-D field
    USE JACCOM_DAT ! Jacobian data
    USE PTBCON_DAT ! Jacobian perturbation sizes
!
! SUBROUTINES
    USE ATMAUX_SUB ! Set up auxiliary profiles of atmospheric parameters
    USE PTBPRF_SUB ! Apply perturbation to a profile
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IJAC ! Index of Jacobian element (1:NJAC) or 0
!
! LOCAL VARIABLES
    LOGICAL     :: FIRST  = .TRUE. ! T=first call of this subroutine
    LOGICAL     :: NEWAUX     ! T=need to re-run ATMAUX
    INTEGER(I4) :: ILOW       ! Atmospheric level# just below perturbed levels
    INTEGER(I4) :: IPSI       ! Counter for horiz.profile locations
    INTEGER(I4) :: IPTB       ! Atmospheric level# of Jacobian retrieval element
    INTEGER(I4) :: IUPP       ! Atmospheric level# just above perturbed levels
    INTEGER(I4) :: IVMR       ! Index of VMR profile 
    INTEGER(I4) :: JDX        ! Index of Jacobian species
    INTEGER(I4) :: JDXPRV = 0 ! Previous value of JDX
    REAL(R4), SAVE, ALLOCATABLE :: SAVATM(:)   ! Saved unperturbed .atm profile
    REAL(R4), SAVE, ALLOCATABLE :: SAVGRA(:,:) ! Saved unperturbed 2D atm prfs
!
! EXECUTABLE CODE ------------------------------------------------------------
!
  NEWAUX = .FALSE.
!
  IF ( FIRST ) THEN
    ALLOCATE ( SAVATM(NATM) )
    IF ( NPSI .GT. 0 ) ALLOCATE ( SAVGRA(NATM,NPSI) )
    FIRST = .FALSE.
  ELSE                                    !  Restore last perturbed profile
    SELECT CASE ( JDXPRV ) 
    CASE ( JDXTEM ) 
      TEMATM = SAVATM
      IF ( NPSI .GT. 0 ) TEMGRA = SAVGRA
      NEWAUX = .TRUE.
    CASE ( JDXPRE ) 
      PREATM = SAVATM
      IF ( NPSI .GT. 0 ) PREGRA = SAVGRA
      NEWAUX = .TRUE.
    CASE DEFAULT
      IVMR = JDXPRV
      IF ( IVMR .EQ. IAXVMR ) THEN
        EXTATM = SAVATM
        IF ( NPSI .GT. 0 ) EXTGRA = SAVGRA
        NEWAUX = .TRUE.
      ELSE IF ( IVMR .GE. 1 .AND. IVMR .LE. NVMR ) THEN
        VMRATM(:,IVMR) = SAVATM
        IF ( NPSI .GT. 0 ) VMRGRA(:,:,IVMR) = SAVGRA
      END IF
    END SELECT  ! Ignore other values of JDXPRV
  END IF
!
! If just undoing last perturbation, check if new aux data reqd and exit
  IF ( IJAC .EQ. 0 ) THEN
    IF ( NEWAUX ) CALL ATMAUX
    RETURN
  END IF
!
  IPTB = JAC(IJAC)%IAT
  ILOW = JAC(IJAC)%ILO
  IUPP = JAC(IJAC)%IUP
  JDX  = JAC(IJAC)%JDX
  JDXPRV = JDX
!
  SELECT CASE ( JDX ) 
  CASE ( JDXTEM ) 
    SAVATM = TEMATM
    CALL PTBPRF ( ILOW, IPTB, IUPP, PTBTEM, .FALSE., HGTATM, TEMATM ) 
    IF ( NPSI .GT. 0 ) THEN
      SAVGRA = TEMGRA
      DO IPSI = 1, NPSI
        CALL PTBPRF ( ILOW, IPTB, IUPP, PTBTEM, .FALSE., HGTATM, TEMGRA(:,IPSI)) 
      END DO
    END IF
    NEWAUX = .TRUE.
  CASE ( JDXPRE )
    SAVATM = PREATM
    CALL PTBPRF ( ILOW, IPTB, IUPP, PTBPRE, .TRUE., HGTATM, PREATM ) 
    IF ( NPSI .GT. 0 ) THEN
      SAVGRA = PREGRA
      DO IPSI = 1, NPSI
        CALL PTBPRF ( ILOW, IPTB, IUPP, PTBPRE, .TRUE., HGTATM, PREGRA(:,IPSI) ) 
      END DO
    END IF
    NEWAUX = .TRUE.
  CASE DEFAULT
    IVMR = JDX
    IF ( IVMR .EQ. IAXVMR ) THEN
      SAVATM = EXTATM
      CALL PTBPRF ( ILOW, IPTB, IUPP, PTBEXT, .FALSE., HGTATM, EXTATM ) 
      IF ( NPSI .GT. 0 ) THEN
        SAVGRA = EXTGRA 
        DO IPSI = 1, NPSI
          CALL PTBPRF ( ILOW, IPTB, IUPP, PTBEXT, .FALSE., HGTATM, EXTGRA(:,IPSI) )
        END DO
      END IF
      NEWAUX = .TRUE.
    ELSE IF ( IVMR .GE. 1 .AND. IVMR .LE. NVMR ) THEN
      SAVATM = VMRATM(:,IVMR)
      CALL PTBPRF ( ILOW, IPTB, IUPP, PTBVMR, .TRUE., HGTATM, VMRATM(:,IVMR) )
      IF ( NPSI .GT. 0 ) THEN
        SAVGRA = VMRGRA(:,:,IVMR) 
        DO IPSI = 1, NPSI
          CALL PTBPRF ( ILOW, IPTB, IUPP, PTBVMR, .TRUE., HGTATM, &
                        VMRGRA(:,IPSI,IVMR) )
        END DO
      END IF
    END IF
  END SELECT
!
  IF ( NEWAUX ) CALL ATMAUX
!
! Also need to allow for the fact that if profile is changed at level ILOW+1
! then the calculation in *layer* ILOW will also be affected
!      IF ( ILOW .GE. 1 ) CHGATM(ILOW) = .TRUE.
!
END SUBROUTINE PTBATM
END MODULE PTBATM_SUB

