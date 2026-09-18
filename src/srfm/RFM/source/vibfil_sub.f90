MODULE VIBFIL_SUB
  IMPLICIT NONE
  LOGICAL, SAVE :: VIB_FIRST_CALL = .TRUE.
  PUBLIC :: VIBFIL
  PUBLIC :: VIBFIL_RESET
CONTAINS
SUBROUTINE VIBFIL ( NAMVIB, FAIL, ERRMSG )
!
! VERSION
!   24APR26 AD Original.
!
! DESCRIPTION
!   Read Vibrational Level Index file
!   Called by DRVHIT for every Vib file in *HIT section
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE VIBCOM_DAT ! Vibrational level indices
    USE GASCOM_DAT, ONLY: IGSMOL ! IGAS for HITRAN index (or 0). 
    USE HITCOM_DAT, ONLY: USEVIB ! T=Use Vib.Index file
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE IDXARR_GEN ! Index of value within an array, else 0
    USE NAMMOL_FNC ! Return molecule name + i[#iso]
    USE NXTREC_SUB ! Load next record from input file
    USE OPNFIL_SUB ! Open input file
    USE WRTLOG_SUB ! Write text message to log file
!   
  IMPLICIT NONE 
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMVIB ! Name of .vib file 
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL       :: ENDSEC ! T= '*' marker or EOF reached (dummmy)
    LOGICAL       :: USECLS ! T=at least one reqd molecule in this class
    LOGICAL       :: USEFIL ! T=file contains at least one reqd molecule
    INTEGER(I4)   :: ICLS   ! Counter for molec. classes in file
    INTEGER(I4)   :: IDM    ! HITRAN index of molecule
    INTEGER(I4)   :: IDX    ! Vib.level index
    INTEGER(I4)   :: IDXCLS ! Index of molec.class (dummy)
    INTEGER(I4)   :: ILEV   ! Counter for vib levels within class
    INTEGER(I4)   :: IMOL   ! Counter for different molecules 
    INTEGER(I4)   :: IOS    ! IOSTAT saved for error message
    INTEGER(I4)   :: IVIB   ! Index of Vib.Temp profile
    INTEGER(I4)   :: JMOL   ! Secondary molecule counter
    INTEGER(I4)   :: JVIB   ! Index in VIB for inserting this class
    INTEGER(I4)   :: MMOL   ! No.molecules within file class
    INTEGER(I4)   :: NCLS   ! No.different molec.classes in file
    INTEGER(I4)   :: NLEV   ! No.of profile levels 
    CHARACTER(200) :: RECORD ! Record read from file
    TYPE(VIBTYP)  :: CLS    ! Molecular class structure read from file
    TYPE(VIBTYP), ALLOCATABLE :: VIBSAV(:) ! Saved version of VIB
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  USEFIL = .FALSE.
! On first call set up list of molecules in VIB
  IF ( VIB_FIRST_CALL ) THEN 
    NMOL = COUNT ( IGSMOL .NE. 0 ) ! No.different molecules reqd
    ALLOCATE ( IDXMOL(NMOL), IVBMOL(NMOL) ) 
    IVBMOL = 0    ! 0=no index data yet assigned for molecule
    IMOL = 0
    DO IDM = 1, SIZE(IGSMOL)
      IF ( IGSMOL(IDM) .NE. 0 ) THEN
        IMOL = IMOL + 1
        IDXMOL(IMOL) = IDM
      END IF
    END DO
    VIB_FIRST_CALL = .FALSE.
  END IF    

!
  CALL OPNFIL ( LUNTMP, NAMVIB, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
! First record should be number of classes in file
  CALL NXTREC ( LUNTMP, RECORD, ENDSEC, FAIL, ERRMSG ) 
  READ ( RECORD, *, IOSTAT=IOS, ERR=900 ) NCLS
!
  DO ICLS = 1, NCLS
    USECLS = .FALSE.
    CALL NXTREC ( LUNTMP, RECORD, ENDSEC, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    READ ( RECORD, *, IOSTAT=IOS, ERR=900 ) IDXCLS, MMOL, NLEV
    DO JMOL = 1, MMOL
      CALL NXTREC ( LUNTMP, RECORD, ENDSEC, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      READ ( RECORD, *, IOSTAT=IOS, ERR=900 ) IDM
      IMOL = IDXARR ( IDXMOL, IDM ) 
      IF ( IMOL .GT. 0 ) THEN
        USECLS = .TRUE.
        USEFIL = .TRUE.
        IVBMOL(IMOL) = -1            ! flag for being updated
      END IF
    END DO
    ALLOCATE ( CLS%IDX(NLEV), CLS%STR(NLEV) ) 
    DO ILEV = 1, NLEV
      CALL NXTREC ( LUNTMP, RECORD, ENDSEC, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      CLS%STR(ILEV) = RECORD(2:16)
      READ ( RECORD(18:), *, IOSTAT=IOS, ERR=900 ) IDX
      CLS%IDX(ILEV) = IDX
    END DO

    IF ( USECLS ) THEN
! Determine if any previously stored class can be replaced since all
! molecules now assigned to this class instead, or add as new class
      JVIB = NVIB + 1
      DO IVIB = 1, NVIB
        IF ( ANY ( IVBMOL .EQ. IVIB ) ) CYCLE  ! still required
        JVIB = IVIB                            ! no longer required
        DEALLOCATE ( VIB(JVIB)%IDX, VIB(JVIB)%STR )  
      END DO
      IF ( JVIB .GT. NVIB ) THEN          ! Add new class
        IF ( NVIB .GT. 0 ) CALL MOVE_ALLOC ( VIB, VIBSAV ) 
        NVIB = NVIB + 1
        ALLOCATE ( VIB(NVIB) )
        IF ( NVIB .GT. 1 ) VIB(1:NVIB-1) = VIBSAV
      END IF
      ALLOCATE ( VIB(JVIB)%IDX(NLEV), VIB(JVIB)%STR(NLEV) )
      VIB(JVIB)%NLV = NLEV
      VIB(JVIB)%IDX = CLS%IDX
      VIB(JVIB)%STR = CLS%STR
      DO IMOL = 1, NMOL
        IDM = IDXMOL(IMOL)
        IF ( IVBMOL(IMOL) .EQ. -1 ) THEN
          IVBMOL(IMOL) = JVIB
          CALL WRTLOG ( 'I-VIBFIL: Loaded' // TRIM(C11INT(NLEV)) // &
                        ' Vib.Lev indices for ' // TRIM(NAMMOL(IDM)) ) 
        END IF
      END DO
    END IF  
    DEALLOCATE ( CLS%IDX, CLS%STR ) 
  END DO
!
  IF ( USEFIL ) THEN
    USEVIB = .TRUE.  ! Flag in HITCOM for vib indices loaded
  ELSE 
    CALL WRTLOG ( 'W-VIBFIL: file not used - no required molecules' )
  END IF
!
  CLOSE ( LUNTMP )
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-VIBFIL: I/O error on VIB data. IOSTAT=', IOS
!
END SUBROUTINE VIBFIL

SUBROUTINE VIBFIL_RESET()
  VIB_FIRST_CALL = .TRUE.
END SUBROUTINE VIBFIL_RESET

END MODULE VIBFIL_SUB
