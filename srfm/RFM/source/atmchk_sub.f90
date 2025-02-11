MODULE ATMCHK_SUB
CONTAINS
SUBROUTINE ATMCHK ( FAIL, ERRMSG )
!
! VERSION
!   30SEP24 AD Checked.
!   18APR22 AD Bug#33 Add check for IAIVMR. Checked.
!   21JUN17 AD Also check *HGT profile loaded. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Check all required atmospheric profiles loaded
!   Called once by DRVATM.
!   Also sets default surface temperature to match base of atmosphere
!
! VARIABLE KINDS
    USE KIND_DAT  
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE FLGCOM_DAT, ONLY: HYDFLG ! T = use hydros. equlbm for plane par.CG paths
    USE TANCOM_DAT, ONLY: LIMTAN ! T=limb-viewing geometry
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IATM   ! Counter for profile levels
    INTEGER(I4)   :: IVMR   ! Counter over absorbing species
    REAL(R4)      :: SUM    ! Summation for VMR total
    CHARACTER(80) :: WRNMSG ! Warning message sent to log file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Check all required profiles have been found
  FAIL = .TRUE.
  IF ( .NOT. SETTEM ) THEN
    ERRMSG = 'F-ATMCHK: No temperature profile supplied'
    RETURN
  ELSE IF ( .NOT. SETPRE ) THEN
    ERRMSG = 'F-ATMCHK: No pressure profile supplied'
    RETURN
  ELSE IF ( .NOT. SETHGT .AND. ( LIMTAN .OR. HYDFLG ) ) THEN
    ERRMSG = 'F-ATMCHK: No altitude profile supplied'
    RETURN
  ELSE
    DO IVMR = 1, NVMR
      IF ( .NOT. SETVMR(IVMR) ) THEN
        ERRMSG = 'F-ATMCHK: No climatology profile for ' // NAMVMR(IVMR)
        RETURN
      END IF
    END DO
  END IF
!
! Check that total VMR doesn't exceed 1 ppv at any level
  DO IATM = 1, NATM
    SUM = 0.0
    DO IVMR = 1, NVMR
      IF ( IVMR .EQ. IAXVMR ) CYCLE       ! exclude extinction profile
      IF ( IVMR .EQ. IAIVMR ) CYCLE       ! exclude 'air' profile 
      SUM = SUM + VMRATM(IATM,IVMR)  
      IF ( SUM .GT. 1.0E6 ) THEN
        WRNMSG = & 
          'W-ATMCHK: Total VMR > 1ppv adding GAS=' // TRIM(NAMVMR(IVMR)) // &
          ' at HGT=' // C9REAL ( HGTATM(IATM) ) // ' km'
        CALL WRTLOG ( WRNMSG )
      END IF
    END DO
  END DO
!
  FAIL = .FALSE.
!
END SUBROUTINE ATMCHK
END MODULE ATMCHK_SUB
