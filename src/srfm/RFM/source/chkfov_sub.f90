MODULE CHKFOV_SUB
CONTAINS
SUBROUTINE CHKFOV ( FAIL, ERRMSG )
!
! VERSION
!   28SEP24 AD Checked.
!   28MAR22 AD Possible bug fixes. Checked.
!   03DEC21 AD Change log message and slightly rearrange logic.
!   27MAR21 AD Checked.
!   05MAR19 AD Copy ITAN when calling ADDTAN to ensure RAD%SKY etc set.
!   30MAY18 AD Bug#7: Assign new tan paths for all FOV elements
!   01MAY17 AD F90 conversion of fovtan.for. Checked.
!
! DESCRIPTION
!   Check FOV tangent heights
!   Called once by DRVCHK if FOV flag is enabled. 
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FOVCOM_DAT ! Field of View data
    USE TANCOM_DAT ! Tangent path data
    USE FLGCOM_DAT, ONLY: CLCFLG, FVZFLG, PTHFLG ! Option flags
!
! SUBROUTINES
    USE ADDTAN_SUB ! Add tangent ray path
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE CHKLIM_SUB ! Check limb-viewing tangent paths
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE 
!
! ARGUMENTS
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: TOLFAC = 0.01 ! Frac. tol. allowed on FOV spacing
!
! LOCAL VARIABLES
    LOGICAL       :: GEOFOV  ! T=Treat FOV as Geom.Tangent Heights.
    INTEGER(I4)   :: IFOV    ! Counter for FOV points
    INTEGER(I4)   :: ITAN    ! Counter for output tangent heights
    INTEGER(I4)   :: JTAN    ! Counter for all tangent heights
    REAL(R4)      :: FOVDIF  ! Min difference wrt tabulated tan paths
    REAL(R4)      :: TANFOV  ! Value of tangent path + FOV point
    CHARACTER(80) :: MESSGE  ! Text message for Log file
    REAL(R4), ALLOCATABLE :: TOLFOV(:) ! Hgt Tol for matching FOV convol pts
    REAL(R4), ALLOCATABLE :: TANPRV(:) ! Previously added FOV paths
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Set up tangent height tolerances for FOV convolution
! NB NFOV must be at least 3 to accomodate two zeros plus one non-zero value.
  ALLOCATE ( TOLFOV(NFOV) )
  IF ( CLCFLG ) THEN         ! If CLC flag, only use exact matches
    TOLFOV = 0.0
  ELSE                       ! Otherwise use TOLFOV * spacing
    DO IFOV = 2, NFOV-1
      TOLFOV(IFOV) = TOLFAC * &
           MIN ( ABS ( FOV(IFOV+1)%ALT - FOV(IFOV)%ALT ), &
                 ABS ( FOV(IFOV)%ALT - FOV(IFOV-1)%ALT )  )  
    END DO
    TOLFOV(1)    = TOLFOV(2)
    TOLFOV(NFOV) = TOLFOV(NFOV-1)
  END IF
!
! Switch off calculation for set of nominal tangent paths unless required for path
! diagnostic output
  TAN(1:NTAN)%CLC = PTHFLG
!
! Add new tangent paths required for FOV convolution as NTAN+1:MTAN
! After settting up FOV paths for first nom.tan path (ITAN=1), check if 
! subsequent FOV tangent paths can re-use any previously defined paths
!
  GEOFOV = .NOT. ( ELEFOV .OR. FVZFLG )
  ALLOCATE ( ITNFOV(NTAN,NFOV) )
  DO ITAN = 1, NTAN
    DO IFOV = 1, NFOV
      IF ( ELEFOV ) THEN          
        TANFOV = TAN(ITAN)%ELE + FOV(IFOV)%ALT
      ELSE IF ( GEOFOV ) THEN
        TANFOV = TAN(ITAN)%GEO + FOV(IFOV)%ALT
      ELSE
        TANFOV = TAN(ITAN)%HGT + FOV(IFOV)%ALT 
      END IF
      IF ( ITAN .GT. 1 ) THEN
        FOVDIF = MINVAL ( ABS ( TANPRV - TANFOV ) ) 
        IF ( FOVDIF .LE. TOLFOV(IFOV) ) THEN ! re-use existing tangent path
          JTAN   = MINLOC ( ABS ( TANPRV - TANFOV ), 1 ) 
          ITNFOV(ITAN,IFOV) = NTAN + JTAN
          CYCLE 
        END IF
      ENDIF
      CALL ADDTAN ( ITAN, .TRUE. )                     ! add new tangent path  
      TAN(MTAN)%USR = TANFOV
      ITNFOV(ITAN,IFOV) = MTAN
      CALL CHKLIM ( MTAN, ELEFOV, GEOFOV, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    END DO
    IF ( ALLOCATED ( TANPRV ) ) DEALLOCATE ( TANPRV )
    ALLOCATE ( TANPRV(MTAN-NTAN) )
    IF ( ELEFOV ) THEN          ! FOV expressed as elevation angles ...
      TANPRV = TAN(NTAN+1:MTAN)%ELE
    ELSE IF ( GEOFOV ) THEN     ! FOV expressed as geom.tan.hts
      TANPRV = TAN(NTAN+1:MTAN)%GEO
    ELSE                        ! FOV expressed in refracted tan.hts
      TANPRV = TAN(NTAN+1:MTAN)%HGT
    END IF
  END DO
!
  MESSGE = 'I-CHKFOV: Add '// TRIM ( C11INT ( MTAN - NTAN ) ) // &
    ' tan.paths for FOV Convl. New total=' // C11INT ( MTAN )
  CALL WRTLOG ( MESSGE )
!
END SUBROUTINE CHKFOV
END MODULE CHKFOV_SUB

