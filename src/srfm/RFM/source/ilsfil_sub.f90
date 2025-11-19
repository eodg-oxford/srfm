MODULE ILSFIL_SUB
CONTAINS
SUBROUTINE ILSFIL ( NAMILS, FAIL, ERRMSG )
!
! VERSION
!   12NOV23 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Read ILS data file
!   Called by DRVILS for each file in *ILS section of driver table
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ILSCOM_DAT ! Instrument Lineshape data
    USE FLGCOM_DAT, ONLY: GHZFLG ! T=use GHz spectral axis, F=wavenumber axis
    USE PHYCON_DAT, ONLY: GHZ2CM ! GHz to cm-1 conv.fac
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE NXTREC_SUB ! Load next record from input file
    USE OPNFIL_SUB ! Open input file  
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE 
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMILS ! Name of ILS data file
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!      
! LOCAL VARIABLES      
    LOGICAL       :: ENDSEC ! T=reached end of file, F=continue
    INTEGER(I4)   :: IOS    ! Saved value of IOSTAT for error message
    INTEGER(I4)   :: NPT    ! No. of points in tabulated ILS function
    REAL(R8)      :: PT1    ! Lower limit of ILS fn [cm-1] - relative 
    REAL(R8)      :: PTD    ! Axis spacing [cm-1] of ILS fn
    REAL(R8)      :: WNL    ! Lower Wno for application of ILS function
    REAL(R8)      :: WNU    ! Upper Wno for application of ILS function
    CHARACTER(80) :: RECORD ! Record read from ILS file
    TYPE(ILSTYP), ALLOCATABLE :: ILSSAV(:) ! Copy of ILS during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL OPNFIL ( LUNTMP, NAMILS, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
  DO                   ! loop over different ILS functions in file
    CALL NXTREC ( LUNTMP, RECORD, ENDSEC, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( ENDSEC ) EXIT                 ! reached end-of-file 
!      
    WNL = 0.0D0
    WNU = 0.0D0
    READ ( RECORD, *, IOSTAT=IOS ) NPT, PT1, PTD, WNL, WNU
    IF ( IOS .GT. 0 ) THEN
      FAIL = .TRUE.
      WRITE ( ERRMSG, * ) &
        'F-ILSFIL: Error reading ILS file header. IOSTAT=', IOS
      RETURN 
    END IF 
!
! See if this is the default ILS fn applicable to any spectral range
    IF ( WNL .EQ. 0.0D0 .AND. WNU .EQ. 0.0D0 ) THEN
      IF ( IDFILS .GT. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-ILSFIL: No range specified for ILS Fn. but ' // &
                 'default ILS Fn already loaded.'
        RETURN        
      ELSE 
        IDFILS = NILS + 1
        CALL WRTLOG ( 'I-ILSFIL: Loading default ILS Fn for any Wno.range' )
      END IF
    END IF 
!
    IF ( ALLOCATED ( ILS ) ) CALL MOVE_ALLOC ( ILS, ILSSAV )
    NILS = NILS + 1
    ALLOCATE ( ILS(NILS) )
    IF ( ALLOCATED ( ILSSAV ) ) ILS(1:NILS-1) = ILSSAV
!
! Warn if ILS file specified in different spectral coords to RFM calcs
    IF ( GHZFLG ) THEN 
      IF ( NPT .GT. 0 ) &
        CALL WRTLOG ( 'W-ILSFIL: ILS specified in Wno - converting to GHz' )
    ELSE
      IF ( NPT .LT. 0 ) &
        CALL WRTLOG ( 'W-ILSFIL: ILS specified in GHz - converting to Wno' )
    END IF
!
! If NPts is negative, convert from GHz to Wno
    IF ( NPT .LT. 0 ) THEN
      NPT = ABS ( NPT ) 
      WNL = GHZ2CM * WNL
      WNU = GHZ2CM * WNU
      PT1 = GHZ2CM * PT1
      PTD = GHZ2CM * PTD
    END IF
!
    ILS(NILS)%NPT = NPT
    ILS(NILS)%PT1 = PT1
    ILS(NILS)%PTD = PTD
    ILS(NILS)%PT2 = PT1 + ( NPT - 1 ) * PTD
    ILS(NILS)%WNL = WNL
    ILS(NILS)%WNU = WNU
    ALLOCATE ( ILS(NILS)%FNC(NPT) )
!
    READ ( LUNTMP, *, IOSTAT=IOS ) ILS(NILS)%FNC 
    IF ( IOS .NE. 0 ) THEN
      FAIL = .TRUE.
      WRITE ( ERRMSG, * ) 'F-ILSFIL: Error reading ILS data. IOSTAT=', IOS
      RETURN
    END IF
!
  END DO
!
  CLOSE ( LUNTMP )
  FAIL = .FALSE.
!
END SUBROUTINE ILSFIL
END MODULE ILSFIL_SUB
