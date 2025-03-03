MODULE DRVOBS_SUB
CONTAINS
SUBROUTINE DRVOBS ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   04OCT24 AD Checked.
!   27SEP22 AD Bug#36: Add check for unrecognised PARAM
!   05AUG19 AD Rewritten for PARAM=VALUE fields. Allow PREOBS. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Read RFM driver table *OBS section
!   Called by RFMDRV following *OBS marker in Driver table.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE OBSCOM_DAT ! Observer location data
    USE FLGCOM_DAT, ONLY: GRAFLG ! T = use horizontal gradients
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
    USE NXTFLD_SUB ! Load next field from section of driver file
    USE OBSCHK_SUB ! Check and set observer altitude
    USE PARFLD_SUB ! Extract Parameter=Value string from record
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for driver table
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: PSIMAX = 90.0 ! Max abs.val[deg]  of obs horiz.angle 
!                                          relative to reference profile
! LOCAL VARIABLES
    LOGICAL           :: GOTHGT = .FALSE. ! T=HGTOBS has been specified
    LOGICAL           :: GOTPAR           ! T=FIELD is PARAM=VALUE format
    LOGICAL           :: GOTPRE = .FALSE. ! T=HGTPRE has been specified
    LOGICAL           :: GOTPSI = .FALSE. ! T=HGTPSI has been specified
    INTEGER(I4)       :: IFLD   = 0       ! Counter for fields
    INTEGER(I4)       :: IOS    ! Saved value of IOSTAT for error message
    INTEGER(I4)       :: LENGTH ! No.of characters written in FIELD 
    REAL(R4)          :: RVAL   ! Value read from VALUE of FIELD
    CHARACTER(LENREC) :: FIELD  ! Data field from driver file
    CHARACTER(6)      :: PARAM  ! Parameter from PARAM=VALUE pair
    CHARACTER(LENREC) :: VALUE  ! Value from PARAM=VALUE pair
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Read next field in *OBS section
  DO 
    CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
!
! Check for PARAM=VALUE form
    CALL PARFLD ( FIELD, GOTPAR, PARAM, VALUE ) 
!
    IF ( .NOT. GOTPAR ) THEN
      IFLD = IFLD + 1
      SELECT CASE ( IFLD ) 
      CASE ( 1 ) 
        PARAM = 'HGTOBS'
      CASE ( 2 ) 
        PARAM = 'PSIOBS'
      CASE DEFAULT
        ERRMSG = 'F-DRVOBS: Unexpected extra (3rd) field in *OBS section'
        FAIL = .TRUE.
        RETURN
      END SELECT
      VALUE = FIELD
    END IF
!
    READ ( VALUE, *, IOSTAT=IOS ) RVAL
    IF ( IOS .NE. 0 ) THEN
      FAIL = .TRUE.
      WRITE ( ERRMSG, * ) 'F-DRVOBS: Error reading ' // PARAM // '=value.' &
             // ' IOSTAT=', IOS
      RETURN
    END IF 
!
    FAIL = .TRUE.
    SELECT CASE ( PARAM ) 
!
    CASE ( 'HGTOBS' )         ! Observer altitude specified
      IF ( GOTHGT ) THEN
        ERRMSG = 'F-DRVOBS: *OBS section has repeated entry for HGTOBS'
      ELSE IF ( GOTPRE ) THEN
        ERRMSG = 'F-DRVOBS: *OBS section cannot have both HGTOBS and PREOBS'
      ELSE 
        GOTHGT = .TRUE.
        CALL OBSCHK ( RVAL, .TRUE., FAIL, ERRMSG )
        HGTOBS = RVAL
        CALL WRTLOG ( 'I-DRVOBS: Setting Observer Altitude = ' // &
                      TRIM ( C9REAL(HGTOBS) ) // ' [km]' )
      END IF
!
    CASE ( 'PREOBS' )          ! Observer pressure specified
      IF ( GOTPRE ) THEN
        ERRMSG = 'F-DRVOBS: *OBS section has repeated entry for PREOBS'
      ELSE IF ( GOTHGT ) THEN
        ERRMSG = 'F-DRVOBS: *OBS section cannot have both PREOBS and HGTOBS'
      ELSE 
        GOTPRE = .TRUE.
        CALL OBSCHK ( RVAL, .FALSE., FAIL, ERRMSG )
        PREOBS = RVAL
        CALL WRTLOG ( 'I-DRVOBS: Setting Observer Pressure = ' // &
                      TRIM ( C9REAL(PREOBS) ) // ' [hPa]' )
      END IF
!
    CASE ( 'PSIOBS' )          ! Observer Horizontal angle specified
      IF ( GOTPSI ) THEN
        ERRMSG = 'F-DRVOBS: *OBS section has repeated entry for PSIOBS'
      ELSE
        GOTPSI = .TRUE.
        IF ( GRAFLG ) THEN
          PSIOBS = RVAL
          IF ( ABS ( PSIOBS ) .GT. PSIMAX ) THEN
            ERRMSG = 'F-DRVOBS: Abs(Obs.Horiz.angle)=' // &
                     C9REAL ( ABS(PSIOBS) ) // ' > PSIMAX=' // &
                     C9REAL ( PSIMAX ) // ' [deg]'
          ELSE 
            CALL WRTLOG ( 'I-DRVOBS: Setting Observer Horizontal Angle = ' // &
                          TRIM ( C9REAL(PREOBS) ) // ' [deg]' )
            FAIL = .FALSE.
          END IF
        ELSE
          CALL WRTLOG ( 'W-DRVOBS: GRA flag not enabled, so ignoring PSIOBS' )
          FAIL = .FALSE.
        END IF
      END IF
    CASE DEFAULT
      ERRMSG = 'F-DRVOBS: Unrecognised PARAM field: ' // PARAM
    END SELECT
    IF ( FAIL ) RETURN
!
  END DO
!
  FAIL = .TRUE.
  IF ( .NOT. GOTHGT .AND. .NOT. GOTPRE ) THEN
    ERRMSG = 'F-DRVOBS: No observer altitude or pressure in *OBS section'
  ELSE IF ( GRAFLG .AND. .NOT. GOTPSI ) THEN
    ERRMSG = 'F-DRVOBS: No observer horizontal angle supplied ' // &
             '(required with GRA flag)'
  ELSE
    FAIL = .FALSE.
  END IF
!
END SUBROUTINE DRVOBS
END MODULE DRVOBS_SUB

