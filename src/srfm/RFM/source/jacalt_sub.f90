MODULE JACALT_SUB
CONTAINS
SUBROUTINE JACALT ( JDX, TGT, FINISH, ALT, FAIL, ERRMSG, ITAN )
!
! VERSION
!   27SEP24 AD Checked.
!   05AUG19 AD Add USEHGT argument for ATMLEV. Check SETHGT. Checked.
!   14JAN19 AD Add TRIM( ) around TGT in error message
!   01MAY17 AD F90 conversion of part of jacchk.for. Checked.
!
! DESCRIPTION
!   Set Jacobian perturbation altitudes
!   Called by DRVJAC, JACTAN
!   Check details of Jacobians and set up related data.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE JACCOM_DAT ! Jacobian data
!
! SUBROUTINES
    USE ATMLEV_SUB ! Find/insert atmospheric level for given altitude.
    USE C9REAL_GEN ! Write real number as C*9 string
    USE NAMGAS_FNC ! Return molecule name + (iso) associated with GASCOM index
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: JDX    ! Index of Jacobian target species
    CHARACTER(*),  INTENT(IN)  :: TGT    ! Target species for Jacobian
    LOGICAL,       INTENT(IN)  :: FINISH ! T=No more altitudes to add
    REAL(R4),      INTENT(IN)  :: ALT    ! Jacobian altitude [km] 
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
    INTEGER(I4), OPTIONAL, &
                   INTENT(IN)  :: ITAN   ! Index of nominal tangent path
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: ALTTOL = 0.001 ! Altitude matching tolerance [km]
!
! LOCAL VARIABLES
    INTEGER(I4) :: IAT  = -1 ! Atmos profile level at Jacobian ptb maximum
    INTEGER(I4) :: ILO  = -1 ! Atmos profile level at Jacobian ptb lower limit
    INTEGER(I4) :: ITN  = 0  ! Tan.ray index for Jacobian
    INTEGER(I4) :: IUP  = -1 ! Atmos profile level at Jacobian ptb upper limit
    INTEGER(I4) :: NALT = 0  ! Number of stored altitudes for current species
    REAL(R4)    :: ALTTMP    ! Previous ALT altitude value [km] 
    TYPE(JACTYP), ALLOCATABLE :: JACSAV(:) ! Saved JAC during reallocation
!
! EXECUTABLE CODE ------------------------------------------------------------
!
  FAIL = .FALSE.
  IF ( FINISH ) THEN
    SELECT CASE ( NALT ) 
    CASE ( 0 ) 
      IF ( ALLOCATED ( JAC ) ) CALL MOVE_ALLOC ( JAC, JACSAV ) 
      NJAC = NJAC + 1
      ALLOCATE ( JAC(NJAC) ) 
      IF ( ALLOCATED ( JACSAV ) ) JAC(1:NJAC-1) = JACSAV
      JAC(NJAC)%JDX = JDX
      JAC(NJAC)%ILO = 0
      JAC(NJAC)%IAT = 1
      JAC(NJAC)%IUP = NATM + 1
      JAC(NJAC)%ITN = 0 
      JAC(NJAC)%COL = .TRUE.
      JAC(NJAC)%HGT = 0.0
      IF ( JDX .GE. 1 .AND. JDX .LE. NVMR ) THEN
        JAC(NJAC)%COD = NAMGAS(JDX)   ! convert eg 'h2o(4)' to 'h2oi4'
      ELSE
        JAC(NJAC)%COD = TGT
      END IF
    CASE ( 1:2 ) 
      ERRMSG = 'F-JACALT: at least 3 altitudes required for Jacobians, ' &
               // 'target=' // TRIM(TGT)
      FAIL = .TRUE.
    CASE DEFAULT              ! reset for next target
      NALT = 0     
      ILO = -1
      IAT = -1
      IUP = -1
      ITN = 0
    END SELECT
  ELSE                        ! update list of altitudes
    ILO = IAT
    IAT = IUP
    NALT = NALT + 1
    IF ( ALT .LE. HGTSFC - ALTTOL ) THEN
      IUP = 0
    ELSE IF ( ALT .GE. HGTTOA + ALTTOL ) THEN
      IUP = NATM + 1
    ELSE   ! Since increasing alts, adding levels shouldn't change ILO or IAT
      ALTTMP = ALT  ! ATMLEV may adjust ALT slightly, so set local variable
      CALL ATMLEV ( ALTTMP, .TRUE., IUP )  ! T = ALTTMP is altitude 
    END IF
    IF ( IUP .LE. IAT ) THEN
      ERRMSG = 'F-JACALT: ' // TRIM(TGT) // ' ptb altitude=' // C9REAL(ALT) // &
               ' not increasing monotonically'
      FAIL = .TRUE.
    ELSE IF ( JDX .EQ. JDXSFT ) THEN
      ERRMSG = 'F-JACALT: SFCTEM retrieval cannot have altitude levels'
      FAIL = .TRUE.
    ELSE IF ( JDX .EQ. JDXSFE ) THEN
      ERRMSG = 'F-JACALT: SFCEMS retrieval cannot have altitude levels'
      FAIL = .TRUE.
    ELSE IF ( NATM .EQ. 1 ) THEN
      ERRMSG = 'F-JACALT: HOMogeneous path Jacobians cannot have alt.levels'
      FAIL = .TRUE.
    ELSE IF ( .NOT. SETHGT ) THEN
      ERRMSG = 'F-JACALT: Jacobian altitudes require a *HGT profile to be set'
      FAIL = .TRUE.
    ELSE IF ( NALT .GE. 3 ) THEN ! start defining alt-dependent Jacobians
      IF ( ALLOCATED ( JAC ) ) CALL MOVE_ALLOC ( JAC, JACSAV ) 
      NJAC = NJAC + 1
      ALLOCATE ( JAC(NJAC) ) 
      IF ( ALLOCATED ( JACSAV ) ) JAC(1:NJAC-1) = JACSAV
      JAC(NJAC)%JDX = JDX
      JAC(NJAC)%ILO = ILO
      JAC(NJAC)%IAT = IAT
      JAC(NJAC)%IUP = IUP
      JAC(NJAC)%ITN = ITN
      JAC(NJAC)%COL = .FALSE.
      JAC(NJAC)%HGT = HGTATM(IAT) 
      IF ( JDX .GE. 1 .AND. JDX .LE. NVMR ) THEN
        JAC(NJAC)%COD = NAMGAS(JDX)
      ELSE
        JAC(NJAC)%COD = TGT
      END IF
    END IF  
! ITAN argument defines tan.ray index for IUP
    IF ( PRESENT(ITAN) ) THEN
      ITN = ITAN
    ELSE
      ITN = 0
    END IF
!
  END IF  
! 
END SUBROUTINE JACALT
END MODULE JACALT_SUB
