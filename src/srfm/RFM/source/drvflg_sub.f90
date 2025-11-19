MODULE DRVFLG_SUB
  USE KIND_DAT
  IMPLICIT NONE
  PRIVATE

  INTEGER(I4), PARAMETER :: MAXFLG = 30 ! Max no.different flags
  INTEGER(I4) :: FLAG_COUNT = 0
  CHARACTER(3) :: FLAG_HISTORY(MAXFLG) = ''

  PUBLIC :: DRVFLG
  PUBLIC :: DRVFLG_RESET
CONTAINS
SUBROUTINE DRVFLG ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   15MAR23 AD Add C32 flag.
!   29APR21 AD Allow SFC+ZEN flags
!   24JUN19 AD Remove CIA flag. Checked.
!   20DEC17 AD F90 conversion of inpflg.for. Checked.
!
! DESCRIPTION
!   Read RFM driver table *FLG section.
!   Called once by RFMDRV. DRVFLG_RESET clears persistent bookkeeping so
!   repeated runs within the same session start with no active flags.
!
! VARIABLE KINDS
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
!
! SUBROUTINES
    USE NXTFLD_SUB ! Load next field from section of driver file
    USE UPCASE_FNC ! Convert text string to upper case
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for Driver File
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL CONSTANTS
! LOCAL VARIABLES
    INTEGER(I4)   :: LENGTH   ! Length of character field read from driver table
    CHARACTER(40) :: FIELD    ! Field read from driver table
    CHARACTER(3)  :: FLG      ! Flag read from driver table (upper case)
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Construct message for log file listing status of all flags
  CALL WRTLOG ( 'I-DRVFLG: Enabled flags:', .TRUE. )
!
! Read first field in *FLG section of driver table
  DO
    CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
!
    IF ( LENGTH .EQ. 0 ) EXIT   ! end of *FLG section reached
!
    IF ( LENGTH .NE. 3 ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-DRVFLG: Flag not a C*3 string, ='//FIELD
      RETURN
    END IF
!
    FLG = UPCASE ( FIELD(1:3) )
    IF ( ANY ( FLAG_HISTORY .EQ. FLG ) ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-DRVFLG: *FLG section contains repeated flag: '//FLG
      RETURN
    END IF
    FLAG_COUNT = FLAG_COUNT + 1
    IF ( FLAG_COUNT .GT. MAXFLG ) STOP 'F-DRVFLG: Logical error'
    FLAG_HISTORY(FLAG_COUNT) = FLG
    CALL WRTLOG ( ' '//FLG, .TRUE. ) 
!
    SELECT CASE ( FLG )
      CASE ( 'ABS' ) ; ABSFLG = .TRUE.
      CASE ( 'AVG' ) ; AVGFLG = .TRUE.
      CASE ( 'BBT' ) ; BBTFLG = .TRUE.
      CASE ( 'BFX' ) ; BFXFLG = .TRUE.
      CASE ( 'BIN' ) ; BINFLG = .TRUE.
      CASE ( 'C32' ) ; C32FLG = .TRUE.
      CASE ( 'CHI' ) ; CHIFLG = .TRUE.
      CASE ( 'CIA' ) ; CALL WRTLOG ( 'W-DRVFLG: CIA flag no longer required' )
      CASE ( 'CLC' ) ; CLCFLG = .TRUE.
      CASE ( 'COO' ) ; COOFLG = .TRUE.
      CASE ( 'CTM' ) ; CTMFLG = .TRUE.
      CASE ( 'DBL' ) ; DBLFLG = .TRUE.
      CASE ( 'FIN' ) ; FINFLG = .TRUE.
      CASE ( 'FLX' ) ; FLXFLG = .TRUE.
      CASE ( 'FOV' ) ; FOVFLG = .TRUE.
      CASE ( 'FVZ' ) ; FVZFLG = .TRUE.
      CASE ( 'GEO' ) ; GEOFLG = .TRUE.
      CASE ( 'GHZ' ) ; GHZFLG = .TRUE.
      CASE ( 'GRA' ) ; GRAFLG = .TRUE.
      CASE ( 'GRD' ) ; GRDFLG = .TRUE.
      CASE ( 'HOM' ) ; HOMFLG = .TRUE.
      CASE ( 'HYD' ) ; HYDFLG = .TRUE.
      CASE ( 'ILS' ) ; ILSFLG = .TRUE.
      CASE ( 'JAC' ) ; JACFLG = .TRUE.
      CASE ( 'JTP' ) ; JTPFLG = .TRUE.
      CASE ( 'LAY' ) ; LAYFLG = .TRUE.
      CASE ( 'LEV' ) ; LEVFLG = .TRUE.
      CASE ( 'LIN' ) ; LINFLG = .TRUE.
      CASE ( 'LOS' ) ; LOSFLG = .TRUE.
      CASE ( 'LUN' ) ; LUNFLG = .TRUE.   ! not used in F90
      CASE ( 'LUT' ) ; LUTFLG = .TRUE.
      CASE ( 'MIX' ) ; MIXFLG = .TRUE.
      CASE ( 'MTX' ) ; MTXFLG = .TRUE.
      CASE ( 'NAD' ) ; NADFLG = .TRUE.
      CASE ( 'NEW' ) ; NEWFLG = .TRUE.
      CASE ( 'NTE' ) ; NTEFLG = .TRUE.
      CASE ( 'OBS' ) ; OBSFLG = .TRUE.
      CASE ( 'OPT' ) ; OPTFLG = .TRUE.
      CASE ( 'PRF' ) ; PRFFLG = .TRUE.
      CASE ( 'PTH' ) ; PTHFLG = .TRUE.
      CASE ( 'QAD' ) ; QADFLG = .TRUE.
      CASE ( 'RAD' ) ; RADFLG = .TRUE.
      CASE ( 'REJ' ) ; REJFLG = .TRUE.
      CASE ( 'REX' ) ; REXFLG = .TRUE.
      CASE ( 'RJT' ) ; RJTFLG = .TRUE.
      CASE ( 'SFC' ) ; SFCFLG = .TRUE.
      CASE ( 'SHH' ) ; SHHFLG = .TRUE.
      CASE ( 'SHP' ) ; SHPFLG = .TRUE.
      CASE ( 'SVD' ) ; SVDFLG = .TRUE.
      CASE ( 'TAB' ) ; TABFLG = .TRUE.
      CASE ( 'TRA' ) ; TRAFLG = .TRUE.
      CASE ( 'VRT' ) ; VRTFLG = .TRUE.
      CASE ( 'VVW' ) ; VVWFLG = .TRUE.
      CASE ( 'WID' ) ; WIDFLG = .TRUE.
      CASE ( 'ZEN' ) ; ZENFLG = .TRUE.     
      CASE DEFAULT
        FAIL = .TRUE.
        ERRMSG = 'F-DRVFLG: *FLG section contains unrecognised flag=' // FLG
        RETURN
    END SELECT
  END DO
!
  CALL WRTLOG ( '', .FALSE. )
!
  IF ( LUNFLG ) CALL WRTLOG ( 'W-DRVFLG: LUN flag redundant in RFM v5' )
! 
! Check for inconsistent flag combinations
  FAIL = .TRUE.
  IF ( ABSFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: ABS and TAB flags are incompatible'
!
  ELSE IF ( AVGFLG .AND. ILSFLG ) THEN
    ERRMSG = 'F-DRVFLG: AVG and ILS flags are incompatible'
  ELSE IF ( AVGFLG .AND. OPTFLG ) THEN
    ERRMSG = 'F-DRVFLG: AVG and OPT flags are incompatible'
  ELSE IF ( AVGFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: AVG and TAB flags are incompatible'
!
  ELSE IF ( BBTFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: BBT and TAB flags are incompatible'
!
  ELSE IF ( BFXFLG .AND. NTEFLG ) THEN
    ERRMSG = 'F-DRVFLG: BFX and NTE flags are incompatible'
  ELSE IF ( BFXFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: BFX and TAB flags are incompatible'
!
  ELSE IF ( CLCFLG .AND. FLXFLG ) THEN
    ERRMSG = 'F-DRVFLG: CLC and FLX flags are incompatible'
  ELSE IF ( CLCFLG .AND. HOMFLG ) THEN
    ERRMSG = 'F-DRVFLG: CLC and HOM flags are incompatible'
  ELSE IF ( CLCFLG .AND. NADFLG ) THEN
    ERRMSG = 'F-DRVFLG: CLC and NAD flags are incompatible'
  ELSE IF ( CLCFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: CLC and TAB flags are incompatible'
  ELSE IF ( CLCFLG .AND. ZENFLG ) THEN
    ERRMSG = 'F-DRVFLG: CLC and ZEN flags are incompatible'
!
  ELSE IF ( COOFLG .AND. VRTFLG ) THEN
    ERRMSG = 'F-DRVFLG: COO and VRT flags are incompatible'
!
  ELSE IF ( DBLFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: DBL and TAB flags are incompatible'
!
  ELSE IF ( FLXFLG .AND. FOVFLG ) THEN
    ERRMSG = 'F-DRVFLG: FLX and FOV flags are incompatible'
  ELSE IF ( FLXFLG .AND. GEOFLG ) THEN
    ERRMSG = 'F-DRVFLG: FLX and GEO flags are incompatible'
  ELSE IF ( FLXFLG .AND. GRAFLG ) THEN
    ERRMSG = 'F-DRVFLG: FLX and GRA flags are incompatible'
  ELSE IF ( FLXFLG .AND. HOMFLG ) THEN
    ERRMSG = 'F-DRVFLG: FLX and HOM flags are incompatible'
  ELSE IF ( FLXFLG .AND. LEVFLG ) THEN
    ERRMSG = 'F-DRVFLG: FLX and LEV flags are incompatible'
  ELSE IF ( FLXFLG .AND. LOSFLG ) THEN
    ERRMSG = 'F-DRVFLG: FLX and LOS flags are incompatible'
  ELSE IF ( FLXFLG .AND. OBSFLG ) THEN
    ERRMSG = 'F-DRVFLG: FLX and OBS flags are incompatible'
  ELSE IF ( FLXFLG .AND. OPTFLG .AND. .NOT. VRTFLG ) THEN
    ERRMSG = 'F-DRVFLG: FLX and OPT flags also require VRT flag'
  ELSE IF ( FLXFLG .AND. PTHFLG ) THEN
    ERRMSG = 'F-DRVFLG: FLX and PTH flags are incompatible'
  ELSE IF ( FLXFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: FLX and TAB flags are incompatible'
!
  ELSE IF ( FOVFLG .AND. HOMFLG ) THEN
    ERRMSG = 'F-DRVFLG: FOV and HOM flags are incompatible'
  ELSE IF ( FOVFLG .AND. NADFLG ) THEN
    ERRMSG = 'F-DRVFLG: FOV and NAD flags are incompatible'
  ELSE IF ( FOVFLG .AND. OPTFLG ) THEN
    ERRMSG = 'F-DRVFLG: FOV and OPT flags are incompatible'
  ELSE IF ( FOVFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: FOV and TAB flags are incompatible'
  ELSE IF ( FOVFLG .AND. ZENFLG ) THEN
    ERRMSG = 'F-DRVFLG: FOV and ZEN flags are incompatible'
!
  ELSE IF ( GEOFLG .AND. HOMFLG ) THEN
    ERRMSG = 'F-DRVFLG: GEO and HOM flags are incompatible'
  ELSE IF ( GEOFLG .AND. NADFLG ) THEN
    ERRMSG = 'F-DRVFLG: GEO and NAD flags are incompatible'
  ELSE IF ( GEOFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: GEO and TAB flags are incompatible'
  ELSE IF ( GEOFLG .AND. ZENFLG ) THEN
    ERRMSG = 'F-DRVFLG: GEO and ZEN flags are incompatible'
!
  ELSE IF ( GRAFLG .AND. HOMFLG ) THEN
    ERRMSG = 'F-DRVFLG: GRA and HOM flags are incompatible'
  ELSE IF ( GRAFLG .AND. NADFLG ) THEN
    ERRMSG = 'F-DRVFLG: GRA and NAD flags are incompatible'
  ELSE IF ( GRAFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: GRA and TAB flags are incompatible'
  ELSE IF ( GRAFLG .AND. ZENFLG ) THEN
    ERRMSG = 'F-DRVFLG: GRA and ZEN flags are incompatible'
!
  ELSE IF ( HOMFLG .AND. JTPFLG ) THEN
    ERRMSG = 'F-DRVFLG: HOM and JTP flags are incompatible'
  ELSE IF ( HOMFLG .AND. LEVFLG ) THEN
    ERRMSG = 'F-DRVFLG: HOM and LEV flags are incompatible'
  ELSE IF ( HOMFLG .AND. LOSFLG ) THEN
    ERRMSG = 'F-DRVFLG: HOM and LOS flags are incompatible'
  ELSE IF ( HOMFLG .AND. NADFLG ) THEN
    ERRMSG = 'F-DRVFLG: HOM and NAD flags are incompatible'
  ELSE IF ( HOMFLG .AND. OBSFLG ) THEN
    ERRMSG = 'F-DRVFLG: HOM and OBS flags are incompatible'
  ELSE IF ( HOMFLG .AND. ZENFLG ) THEN
    ERRMSG = 'F-DRVFLG: HOM and ZEN flags are incompatible'
!
  ELSE IF ( ILSFLG .AND. OPTFLG ) THEN
    ERRMSG = 'F-DRVFLG: ILS and OPT flags are incompatible'
  ELSE IF ( ILSFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: ILS and TAB flags are incompatible'
!
  ELSE IF ( JACFLG .AND. LEVFLG ) THEN 
    ERRMSG = 'F-DRVFLG: JAC and LEV flags are incompatible'
  ELSE IF ( JACFLG .AND. MTXFLG ) THEN 
    ERRMSG = 'F-DRVFLG: JAC and MTX flags are incompatible'
  ELSE IF ( JACFLG .AND. TABFLG ) THEN 
    ERRMSG = 'F-DRVFLG: JAC and TAB flags are incompatible'
!
  ELSE IF ( JTPFLG .AND. NADFLG ) THEN 
    ERRMSG = 'F-DRVFLG: JTP and NAD flags are incompatible'
  ELSE IF ( JTPFLG .AND. ZENFLG ) THEN 
    ERRMSG = 'F-DRVFLG: JTP and ZEN flags are incompatible'
!
  ELSE IF ( LEVFLG .AND. LOSFLG ) THEN
    ERRMSG = 'F-DRVFLG: LEV and LOS flags are incompatible'
  ELSE IF ( LEVFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: LEV and TAB flags are incompatible'
!
  ELSE IF ( LOSFLG .AND. MTXFLG ) THEN
    ERRMSG = 'F-DRVFLG: LOS and MTX flags are incompatible'
  ELSE IF ( LOSFLG .AND. NADFLG ) THEN
    ERRMSG = 'F-DRVFLG: LOS and NAD flags are incompatible'
  ELSE IF ( LOSFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: LOS and TAB flags are incompatible'
  ELSE IF ( LOSFLG .AND. ZENFLG ) THEN
    ERRMSG = 'F-DRVFLG: LOS and ZEN flags are incompatible'
!
  ELSE IF ( MTXFLG .AND. NADFLG ) THEN
    ERRMSG = 'F-DRVFLG: MTX and NAD flags are incompatible'
  ELSE IF ( MTXFLG .AND. ZENFLG ) THEN
    ERRMSG = 'F-DRVFLG: MTX and ZEN flags are incompatible'
!
  ELSE IF ( NADFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: NAD and TAB flags are incompatible'
  ELSE IF ( NADFLG .AND. ZENFLG ) THEN
    ERRMSG = 'F-DRVFLG: NAD and ZEN flags are incompatible'
!
  ELSE IF ( OBSFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: OBS and TAB flags are incompatible'
!
  ELSE IF ( OPTFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: OPT and TAB flags are incompatible'
!
  ELSE IF ( PTHFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: PTH and TAB flags are incompatible'
!
  ELSE IF ( RADFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: RAD and TAB flags are incompatible'
!
  ELSE IF ( REXFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: REX and TAB flags are incompatible'
!
  ELSE IF ( RJTFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: RJT and TAB flags are incompatible'
!
  ELSE IF ( SFCFLG .AND. TABFLG ) THEN
    ERRMSG = 'F-DRVFLG: SFC and TAB flags are incompatible'
!
  ELSE IF ( TABFLG .AND. TRAFLG ) THEN
    ERRMSG = 'F-DRVFLG: TAB and TRA flags are incompatible'
!
  ELSE IF ( TABFLG .AND. ZENFLG ) THEN
    ERRMSG = 'F-DRVFLG: TAB and ZEN flags are incompatible'
!
  ELSE IF ( BFXFLG .AND. .NOT. ( BBTFLG .OR. RADFLG .OR. RJTFLG ) ) THEN
    ERRMSG = 'F-DRVFLG: BFX flag also requires BBT, RAD or RJT flags'
  ELSE IF ( C32FLG .AND. .NOT. CTMFLG ) THEN
    ERRMSG = 'F-DRVFLG: C32 flag also requires CTM flag'
  ELSE IF ( COOFLG .AND. .NOT. FLXFLG ) THEN
    ERRMSG = 'F-DRVFLG: COO flag also requires FLX flag'
  ELSE IF ( FINFLG .AND. .NOT. ( ILSFLG .OR. AVGFLG ) ) THEN
    ERRMSG = 'F-DRVFLG: FIN flag also requires ILS or AVG flag'
  ELSE IF ( FLXFLG .AND. ABSFLG .AND. .NOT. MTXFLG &
            .AND. .NOT. ( NADFLG .OR. ZENFLG ) ) THEN
    ERRMSG = 'F-DRVFLG: ABS+FLX-MTX flags also requires ZEN or NAD'
  ELSE IF ( FLXFLG .AND. TRAFLG .AND. .NOT. MTXFLG &
            .AND. .NOT. ( NADFLG .OR. ZENFLG ) ) THEN
    ERRMSG = 'F-DRVFLG: TRA+FLX-MTX flags also requires ZEN or NAD'
  ELSE IF ( FLXFLG .AND. .NOT. ZENFLG .AND. .NOT. SFCFLG ) THEN
    ERRMSG = 'F-DRVFLG: FLX-ZEN flags also require SFC flag'
  ELSE IF ( FVZFLG .AND. .NOT. FOVFLG ) THEN
    ERRMSG = 'F-DRVFLG: FVZ flag also requires FOV flag'
  ELSE IF ( HYDFLG .AND. .NOT. ( NADFLG .OR. ZENFLG ) ) THEN
    ERRMSG = 'F-DRVFLG: HYD flag also requires NAD or ZEN flag'
  ELSE IF ( JTPFLG .AND. .NOT. JACFLG ) THEN
    ERRMSG = 'F-DRVFLG: JTP flag also requires JAC flag'
  ELSE IF ( MTXFLG .AND. .NOT. FLXFLG ) THEN
    ERRMSG = 'F-DRVFLG: MTX flag also requires FLX flag'
  ELSE IF ( MTXFLG .AND. RADFLG .AND. .NOT. SFCFLG ) THEN
    ERRMSG = 'F-DRVFLG: MTX+RAD flags also requires SFC flag'
  ELSE IF ( NADFLG .AND. .NOT. SFCFLG ) THEN
    ERRMSG = 'F-DRVFLG: NAD flag also requires SFC flag'
  ELSE IF ( VRTFLG .AND. .NOT. FLXFLG ) THEN
    ERRMSG = 'F-DRVFLG: VRT flag also requires FLX flag'
  ELSE
    FAIL = .FALSE.
  END IF
!
END SUBROUTINE DRVFLG

SUBROUTINE DRVFLG_RESET()
    USE FLGCOM_DAT, ONLY: FLGCOM_RESET
    FLAG_COUNT   = 0
    FLAG_HISTORY = ''
    CALL FLGCOM_RESET()
END SUBROUTINE DRVFLG_RESET

END MODULE DRVFLG_SUB
