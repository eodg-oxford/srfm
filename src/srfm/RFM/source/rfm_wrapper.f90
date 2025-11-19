SUBROUTINE rfm_run(status, run_id_in, driver_lines, enable_capture)
!
! Wrapper around the main RFM executable logic so it can be called repeatedly
! from Python via f2py. Returns 0 on success or a non-zero status code on
! failure. Prior to each invocation shared-module state is reset so that
! sequential runs within the same interpreter see a clean Fortran environment.
!
    USE KIND_DAT
    USE FLGCOM_DAT
    USE RUN_ID_DAT
    USE SPCCOM_DAT, ONLY: SPCCOM_RESET, NSPC, SPC
    USE GASCOM_DAT, ONLY: GASCOM_RESET
    USE ATMCOM_DAT, ONLY: ATMCOM_RESET
    USE TANCOM_DAT, ONLY: TANCOM_RESET
    USE LEVCHK_SUB, ONLY: LEVCHK_RESET
    USE PTHCOM_DAT, ONLY: PTHCOM_RESET
    USE HITCOM_DAT, ONLY: HITCOM_RESET
    USE ILSCOM_DAT, ONLY: ILSCOM_RESET
    USE XSCCOM_DAT, ONLY: XSCCOM_RESET
    USE HFLCOM_DAT, ONLY: HFLCOM_RESET
    USE LEVCOM_DAT, ONLY: LEVCOM_RESET
    USE HDRCOM_DAT, ONLY: VIDHDR, HDRCOM_RESET
    USE RFMLUN_DAT, ONLY: LUNLOG, RFMLUN_RESET
    USE LENREC_DAT, ONLY: LENREC
    USE RFMDAL_SUB
    USE RFMDRV_SUB
    USE RFMPRF_SUB
    USE RFMPTH_SUB
    USE RFMSPC_SUB
    USE WRTLOG_SUB
    USE PYOUT_SUB
    USE DRVKEY_SUB, ONLY: DRVKEY_RESET
    USE DRVFLG_SUB, ONLY: DRVFLG_RESET
    USE DRVBUF_DAT, ONLY: DRVBUF_SET, DRVBUF_CLEAR
    USE ATMFIL_SUB, ONLY: ATMFIL_RESET
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(OUT) :: status
    CHARACTER(*), INTENT(IN), OPTIONAL :: run_id_in
    CHARACTER(LENREC), INTENT(IN), OPTIONAL :: driver_lines(:)
    LOGICAL, INTENT(IN), OPTIONAL :: enable_capture
!
! LOCAL CONSTANTS
    LOGICAL, PARAMETER :: PROMPT = .FALSE.
!
! LOCAL VARIABLES
    LOGICAL       :: fail
    LOGICAL       :: capture_enabled
    INTEGER(I4)   :: ios
    INTEGER(I4)   :: ispc
    CHARACTER(80) :: errmsg
    CHARACTER(80) :: logmsg
!
! ---------------------------------------------------------------------------
!
  status = 0
  fail   = .FALSE.
  errmsg = ''
  capture_enabled = .FALSE.
  IF ( PRESENT ( enable_capture ) ) THEN
    capture_enabled = enable_capture
  END IF
!
! Reset any prior inline driver data before each run.
  CALL DRVBUF_CLEAR()
  IF ( PRESENT ( driver_lines ) ) THEN
    IF ( SIZE ( driver_lines ) .GT. 0 ) CALL DRVBUF_SET ( driver_lines )
  END IF
  CALL RFMLUN_RESET()
  CALL DRVKEY_RESET()
  CALL DRVFLG_RESET()
  CALL HDRCOM_RESET()
  CALL SPCCOM_RESET()
  CALL GASCOM_RESET()
  CALL ATMCOM_RESET()
  CALL ATMFIL_RESET()
  CALL TANCOM_RESET()
  CALL LEVCOM_RESET()
  CALL LEVCHK_RESET()
  CALL PTHCOM_RESET()
  CALL HITCOM_RESET()
  CALL HFLCOM_RESET()
  CALL ILSCOM_RESET()
  CALL XSCCOM_RESET()
!
  IF ( PRESENT(run_id_in) ) THEN
    RUN_ID = TRIM(run_id_in)
  ELSE
    RUN_ID = ''
  END IF
!
  VIDHDR = '5.21_23SEP'
  logmsg = 'R-RFM: Running RFM v' // VIDHDR
  WRITE ( *, '(A)' ) logmsg
!
  IF ( PROMPT ) THEN
    WRITE ( *, '(A)' ) 'Optional ID to be appended to filenames (<CR>=none):'
    READ ( *, '(A)' ) RUN_ID
    IF ( RUN_ID .NE. '' ) &
      WRITE ( *, '(A)' ) 'R-RFM:    Filename append string=' // RUN_ID
  END IF
!
! Open log file
  OPEN ( UNIT=LUNLOG, FILE='rfm.log'//RUN_ID, ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=ios )
  IF ( ios .NE. 0 ) THEN
    WRITE ( *, '(A)' ) 'F-RFM: Error opening rfm.log file. IOSTAT=', ios
    status = ios
    RETURN
  END IF
  CALL WRTLOG ( logmsg )
!
! Read driver table contents
  CALL RFMDRV ( fail, errmsg )
  IF ( fail ) GOTO 900
  IF ( capture_enabled ) CALL PYOUT_RESET ( NSPC )
!
! If PRF flag, output atmospheric profile being used
  IF ( PRFFLG ) THEN
    CALL RFMPRF ( fail, errmsg )
    IF ( fail ) GOTO 900
  END IF
!
! Calculate equivalent CG paths
  CALL RFMPTH ( fail, errmsg )
  IF ( fail ) GOTO 900
!
! Loop over required spectral ranges
  DO ispc = 1, NSPC
    IF ( NSPC .GT. 1 ) THEN
      logmsg = 'I-RFM: Calculation for spectral range: ' // SPC(ispc)%LAB
      IF ( .NOT. SHHFLG ) WRITE ( *, '(A)' ) logmsg
      CALL WRTLOG ( logmsg )
    END IF
    CALL RFMSPC ( ispc, fail, errmsg )
    IF ( fail ) GOTO 900
    IF ( capture_enabled ) CALL PYOUT_COLLECT_SPECTRUM ( ispc )
  END DO
!
! Deallocate pointers
  CALL RFMDAL
  CALL LEVCOM_RESET()
  CALL LEVCHK_RESET()
!
900 CONTINUE
  IF ( fail ) THEN
    logmsg = errmsg
    status = 1
  ELSE
    logmsg = 'R-RFM: Successful completion'
    status = 0
  END IF
  CALL WRTLOG ( logmsg )
  WRITE ( *, '(A)' ) logmsg
  CLOSE ( UNIT=LUNLOG )
!
END SUBROUTINE rfm_run

SUBROUTINE rfm_get_optical_grid_size(status, n_levels, n_points, ispc)
!
! Return the dimensions of the captured optical-depth grid.
!
    USE KIND_DAT
    USE PYOUT_SUB, ONLY: PYOUT_GET_SPECTRUM_SIZE
!
    INTEGER(I4), INTENT(OUT) :: status
    INTEGER(I4), INTENT(OUT) :: n_levels
    INTEGER(I4), INTENT(OUT) :: n_points
    INTEGER(I4), INTENT(IN), OPTIONAL :: ispc
!
    INTEGER(I4) :: idx
!
  IF ( PRESENT ( ispc ) ) THEN
    idx = ispc
  ELSE
    idx = 1
  END IF
!
  CALL PYOUT_GET_SPECTRUM_SIZE ( idx, status, n_levels, n_points )
!
END SUBROUTINE rfm_get_optical_grid_size
!
SUBROUTINE rfm_get_optical_grid(status, n_levels, n_points, altitudes, &
                                pressures, temperatures, wavenumbers, &
                                cumulative, ispc)
!
! Fill the provided arrays with the captured optical-depth grid.
!
    USE KIND_DAT
    USE PYOUT_SUB, ONLY: PYOUT_FILL_SPECTRUM
!
    INTEGER(I4), INTENT(OUT) :: status
    INTEGER(I4), INTENT(IN) :: n_levels
    INTEGER(I4), INTENT(IN) :: n_points
    REAL(R8), INTENT(OUT) :: altitudes(n_levels)
    REAL(R8), INTENT(OUT) :: pressures(n_levels)
    REAL(R8), INTENT(OUT) :: temperatures(n_levels)
    REAL(R8), INTENT(OUT) :: wavenumbers(n_points)
    REAL(R8), INTENT(OUT) :: cumulative(n_levels, n_points)
    INTEGER(I4), INTENT(IN), OPTIONAL :: ispc
!
    INTEGER(I4) :: idx
!
  IF ( PRESENT ( ispc ) ) THEN
    idx = ispc
  ELSE
    idx = 1
  END IF
!
  CALL PYOUT_FILL_SPECTRUM ( idx, status, altitudes, pressures, &
                             temperatures, wavenumbers, cumulative )
!
END SUBROUTINE rfm_get_optical_grid
