MODULE DRVTAN_SUB
CONTAINS
SUBROUTINE DRVTAN ( LUNDRV, KEY, FAIL, ERRMSG )
!
! VERSION
!   14JUL26 AD Checked.
!   01AUG25 AD Add KEY argument. Simplified to just deal with tan.paths
!   20MAY24 AD Checked.
!   25MAR19 AD Add PARFLD, TANUNI. Checked.
!   04FEB19 AD Add subroutine TANLOS
!   01MAY17 AD F90 conversion of inptan.for. Tested.
!
! DESCRIPTION
!   Read RFM driver table *TAN/*GEO/*ELE section
!   Called by RFMDRV once.
!   Reads inputs for tangent heights.
!
! VARIABLE KINDS
    USE KIND_DAT
! 
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE TANCOM_DAT ! Tangent path data
    USE ATMCOM_DAT, ONLY: SETHGT ! T=height profile has been specified
    USE FLGCOM_DAT, ONLY: LOSFLG ! T = Calculate elev. pointing Jacobian spectra
!
! SUBROUTINES
    USE NXTFFL_SUB ! Load next field from rfm.drv, expanding filenames
    USE TANFLD_SUB ! Check if string is valid *TAN entry and insert in TANCOM
    USE TANLOS_SUB ! Set up LOS Jacobians
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for Driver File
    CHARACTER(4),  INTENT(IN)  :: KEY    ! '*TAN', '*ELE' or '*GEO'
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)       :: ITAN   ! Counter for tangent paths
    INTEGER(I4)       :: LENGTH ! Length of field read from driver table
    REAL(R4)          :: VALUE  ! numerical value read from FIELD
    CHARACTER(LENREC) :: FIELD  ! Tangent height or Tan.Hgt filename
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  IF ( .NOT. SETHGT ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVTAN: Limb-viewing requires height profile in *ATM section'
    RETURN
  END IF
!
  CALL WRTLOG ( 'I-DRVTAN: ', .TRUE. ) 
!
! Read each field in *TAN section
  DO
    CALL NXTFFL ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
    CALL TANFLD ( FIELD, VALUE, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    CALL WRTLOG ( ' '//FIELD, .TRUE. ) 
  END DO
!
  IF ( LOSFLG ) CALL TANLOS
!
! Check at least one tangent height supplied
  IF ( NTAN .EQ. 0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVTAN: No entries in ' // KEY // ' section'
  ELSE IF ( LOSFLG .AND. NTAN .EQ. 1 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVTAN: LOS flag requires at least 2 tangent heights'
  ELSE
    CALL WRTLOG ( '', .FALSE. )
    IF ( LOSFLG .AND. NTAN .EQ. 2 ) CALL WRTLOG ( 'W-INPTAN: Only 2 tan.hts, ' &
                              // 'so LOS Jacobians from linear interpolation' )
  END IF
!
  LIMTAN = .TRUE.
  USRELE = KEY .EQ. '*ELE'
  USRGEO = KEY .EQ. '*GEO'
!
  TAN%CLC = .TRUE.
  TAN%SEC = 1.0   ! Not used for scaling entire path
  DO ITAN = 1, NTAN
    TAN(ITAN)%ITN = ITAN
  END DO
!
END SUBROUTINE DRVTAN
END MODULE DRVTAN_SUB

