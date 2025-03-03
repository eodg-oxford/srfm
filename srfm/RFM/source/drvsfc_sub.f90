MODULE DRVSFC_SUB
CONTAINS
SUBROUTINE DRVSFC ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   22MAY24 AD Checked.
!   29APR21 AD Allow for being called with ZEN flag for upward viewing
!   12APR19 AD Ensure changing surface level also adjusts TEMREL. Checked.
!   05MAR19 AD Add RFLSFC parameter. 
!   14DEC17 AD Restructured to allow PARAM=VALUE specifications as well.
!   01MAY17 AD Original. Checked.
!
! DESCRIPTION
!   Read RFM driver table *SFC section
!   Called once by RFMDRV if *SFC section is present and SFC flag enabled
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE SFCCOM_DAT ! Surface parameters
    USE ATMCOM_DAT, ONLY: IATSFC, TEMATM ! Atmospheric temperature profile
    USE FLGCOM_DAT, ONLY: ZENFLG ! T=upward viewing
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
    USE NXTFLD_SUB ! Load next field from section of driver file
    USE PARFLD_SUB ! Extract Parameter=Value string from record
    USE SFCEMS_SUB ! Read surface emissivity data
    USE SFCLEV_SUB ! Read surface height or pressure
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
! LOCAL VARIABLES
    LOGICAL           :: GOTEMS = .FALSE. ! T=surface emissivity specified
    LOGICAL           :: GOTLEV = .FALSE. ! T=surface hgt/press level set
    LOGICAL           :: GOTPAR           ! T=FIELD is PARAM=VALUE format
    LOGICAL           :: GOTREL = .FALSE. ! T=value for TEMREL specified
    LOGICAL           :: GOTRFL = .FALSE. ! T=surface reflec. model specified
    LOGICAL           :: GOTTEM = .FALSE. ! T=surface temperature specified
    LOGICAL           :: LRELTV = .FALSE. ! T=relative temperature specified
    INTEGER(I4)       :: IFLD   = 1       ! Counter for fields within section
    INTEGER(I4)       :: IOS    = 0       ! Saved value of IOSTAT for error message
    INTEGER(I4)       :: LENGTH           ! Length of FIELD
    REAL(R4)          :: TEMREL           ! Relative temperature [K] (wrt Atmos.T)
    CHARACTER(LENREC) :: FIELD            ! Field read from Driver file
    CHARACTER(6)      :: PARAM            ! Parameter from PARAM=VALUE pair
    CHARACTER(LENREC) :: VALUE            ! Value from PARAM=VALUE pair
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Read next field from *SFC section
  DO
    CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
!
! Check for PARAM=VALUE form
    CALL PARFLD ( FIELD, GOTPAR, PARAM, VALUE ) 
!
! If not, then old structure is to identify fields by position in sequence
    IF ( .NOT. GOTPAR ) THEN 
      IF ( UPCASE ( FIELD ) .EQ. 'REL' ) THEN
        IF ( IFLD .NE. 1 ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-DRVSFC: ''REL'' should be the first field in *SFC section'
          RETURN
        END IF
        LRELTV = .TRUE.
        CYCLE
      ELSE IF ( IFLD .EQ. 1 ) THEN
        IF ( LRELTV ) THEN
          PARAM = 'TEMREL'
        ELSE
          PARAM = 'TEMSFC'
        END IF
        VALUE = FIELD
        IFLD = 2
      ELSE IF ( IFLD .EQ. 2 ) THEN
        PARAM = 'EMSSFC'
        VALUE = FIELD
        IFLD = 3
      ELSE IF ( IFLD .EQ. 3 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-DRVSFC: Unexpected extra field in *SFC section: ' // &
                  TRIM ( FIELD(1:30) )
        RETURN
      END IF
    END IF
!
! Check for surface parameters which have no effect if upward-viewing
    IF ( ZENFLG ) THEN
      IF ( PARAM .EQ. 'TEMREL' .OR. PARAM .EQ. 'TEMSFC' .OR. &
           PARAM .EQ. 'EMSSFC' .OR. PARAM .EQ. 'RFLSFC' ) THEN
        CALL WRTLOG ( 'W-DRVSFC: ignoring ' // PARAM //  &  
                      ' - no effect with ZEN flag' )
        PARAM = 'ignore'
      END IF
    END IF
!
    FAIL = .FALSE.
    SELECT CASE ( PARAM ) 
    CASE ( 'TEMREL' )
      IF ( GOTTEM ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-DRVSFC: Surface temperature already specified'
      ELSE
        GOTTEM = .TRUE.
        GOTREL = .TRUE.
        READ ( VALUE, *, IOSTAT=IOS, ERR=900 ) TEMREL
        CALL WRTLOG ( 'I-DRVSFC: Reading Rel. Surface Temperature = ' // &
                       C9REAL(TEMREL) // ' K' ) 
       END IF
    CASE ( 'TEMSFC' ) 
      IF ( GOTTEM ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-DRVSFC: Surface temperature already specified'
        RETURN
      ELSE
        GOTTEM = .TRUE.
        READ ( VALUE, *, IOSTAT=IOS, ERR=900 ) TEMSFC
        CALL WRTLOG ( 'I-DRVSFC: Reading Surface Temperature = ' // &
                      C9REAL(TEMSFC) // ' K' )
      END IF
    CASE ( 'EMSSFC' ) 
      IF ( GOTEMS ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-DRVSFC: Surface emissivity already specified'
        RETURN
      ELSE 
        GOTEMS = .TRUE.
        CALL SFCEMS ( VALUE, FAIL, ERRMSG )
      END IF
    CASE ( 'HGTSFC' ) 
      IF ( GOTLEV ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-DRVSFC: Surface height/pressure already specified'
        RETURN
      ELSE
        GOTLEV = .TRUE.
        CALL SFCLEV ( PARAM, VALUE, FAIL, ERRMSG )
      END IF
    CASE ( 'PRESFC' ) 
      IF ( GOTLEV ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-DRVSFC: Surface height/pressure already specified'
        RETURN
      ELSE
        GOTLEV = .TRUE.
        CALL SFCLEV ( PARAM, VALUE, FAIL, ERRMSG )
      END IF
    CASE ( 'RFLSFC' ) 
      IF ( GOTRFL ) THEN 
        FAIL = .TRUE.
        ERRMSG = 'F-DRVSFC: Surface reflectivity model already specified'
        RETURN
      ELSE 
        GOTRFL = .TRUE.
        SELECT CASE ( UPCASE ( VALUE ) ) 
        CASE ( 'DIFFUSE' ) 
          DIFSFC = .TRUE.
          CALL WRTLOG ( 'I-DRVSFC: Setting surface reflectivity=diffuse' )
        CASE ( 'SPECULAR' ) 
          DIFSFC = .FALSE.
          CALL WRTLOG ( 'I-DRVSFC: Setting surface reflectivity=specular' )
        CASE DEFAULT
          FAIL = .TRUE.
          ERRMSG = 'F-DRVSFC: Unrecognised value for RFLSFC='//FIELD(1:40)
          RETURN
        END SELECT  
      END IF
    CASE ( 'ignore' )
      CONTINUE
    CASE DEFAULT
      FAIL = .TRUE.
      ERRMSG = 'F-DRVSFC: Unrecognised parameter: ' // PARAM
    END SELECT      
    IF ( FAIL ) RETURN
  END DO
!
! Final part deals with setting/checking surface temperature
!
! If relative temperature specified, calculate from (new) atmos.base level
  IF ( GOTREL ) THEN
    TEMSFC = TEMATM(IATSFC) + TEMREL
    CALL WRTLOG ( 'I-DRVSFC: Lower atm.temp=' // &
             C9REAL(TEMATM(IATSFC)) // ' K so TEMREL gives Tsfc=' // &
             C9REAL(TEMSFC) // ' K' ) 
! If surface level adjusted but no surface temperature specified, adjust 
! default to new base of atmosphere value
  ELSE IF ( GOTLEV .AND. .NOT. GOTTEM ) THEN 
    TEMSFC = TEMATM(IATSFC)
    CALL WRTLOG ( 'I-DRVSFC: New surface level gives adjusted Tsfc=' // &
                  C9REAL(TEMSFC) // ' K' )
  END IF
!
! Check if within allowed range (1000K arbitrary upper limit)
  IF ( TEMSFC .LT. 0.0 .OR. TEMSFC .GT. 1000.0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVSFC: Surface Temperature outside allowed range 0:1000 K' 
    RETURN
  END IF
!
! Warn if large surface-atmosphere temperature discontinuity
  TEMREL = TEMSFC - TEMATM(IATSFC) 
  IF ( ABS ( TEMREL ) .GT. 10.0 ) &  ! 10K is arbitrary
    CALL WRTLOG ( 'W-DRVSFC: Large SFC-ATM Temperature discontinuity, =' &
                  // C9REAL(TEMREL) // ' K' )
!
900 CONTINUE
  FAIL = IOS .NE. 0
  IF ( FAIL ) WRITE ( ERRMSG, * ) 'F-DRVSFC: Error reading ' // PARAM // &
    '. IOSTAT=', IOS
!
END SUBROUTINE DRVSFC
END MODULE DRVSFC_SUB
