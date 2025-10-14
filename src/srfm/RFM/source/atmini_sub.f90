MODULE ATMINI_SUB
CONTAINS
SUBROUTINE ATMINI ( LEVGRD, PREGRD )
!
! VERSION
!   02SEP24 AD Checked.
!   24AUG23 AD Bug#41: Allocate EXTATM
!   18APR22 AD Add IAIVMR and initialise for air
!   02DEC19 AD Correction: set LINVMR=T for aerosol. Checked.
!   21JUN17 AD Allow argument to be lnp or alt.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Initialise atmospheric profile data in ATMCOM
!   Called once by ATMGRD or by ATMPAR (with HOMFLG)
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE GASCOM_DAT ! Molecule and isotope data
    USE FLGCOM_DAT, ONLY: LINFLG ! T = Assume VMR varies linearly with altitude
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: LEVGRD(:) ! Profile grid levels (km or ln(p) [mb])
    LOGICAL,  INTENT(IN) :: PREGRD    ! T=pressure grid, F=height grid
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NATM = SIZE ( LEVGRD )
  FIXPRE = PREGRD
!
  ALLOCATE ( HGTATM(NATM) ) 
  ALLOCATE ( TEMATM(NATM) ) 
  ALLOCATE ( PREATM(NATM) ) 
!
  IF ( FIXPRE ) THEN
    PREATM = LEVGRD 
    SETPRE = .TRUE.
    PRESFC = PREATM(1)
    ALLOCATE ( LNPATM(NATM) ) 
    LNPATM = LOG ( PREATM ) 
    LEVATM => LNPATM
  ELSE
    HGTATM = LEVGRD
    SETHGT = .TRUE.
    HGTTOA = HGTATM(NATM)
    HGTSFC = HGTATM(1)
    LEVATM => HGTATM
  END IF
!
  NVMR = NGAS
  ALLOCATE ( SETVMR(NVMR) ) ; SETVMR = .FALSE.
  ALLOCATE ( LINVMR(NVMR) ) ; LINVMR = LINFLG
  IAXVMR = IAXGAS
  IF ( IAXVMR .GT. 0 ) THEN
    LINVMR(IAXVMR) = .TRUE.  ! Aer always linear interp.
    ALLOCATE ( EXTATM(NATM) )
  END IF
  ALLOCATE ( NAMVMR(NVMR) ) ; NAMVMR = GAS%COD
  ALLOCATE ( NTEVMR(NVMR) ) ; NTEVMR = GAS%NTE
!
  ALLOCATE ( VMRATM(NATM,NVMR) )
!
  IAIVMR = IAIGAS
  IF ( IAIVMR .GT. 0 ) THEN
    LINVMR(IAIVMR) = .FALSE.                    ! Air always log interp.
    VMRATM(:,IAIVMR) = 1.0E6                    ! 'air' mixing ratio = 1ppv
    SETVMR(IAIVMR) = .TRUE.
  END IF
!
END SUBROUTINE ATMINI
END MODULE ATMINI_SUB
