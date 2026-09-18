MODULE DRVHGT_SUB
CONTAINS
SUBROUTINE DRVHGT ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   09JUL26 AD Checked.
!   01AUG25 AD Original. Extracted from DRVTAN
!
! DESCRIPTION
!   Read RFM driver table *HGT section
!   Called by RFMDRV once.
!   Reads inputs for tangent heights, sec zenith angles, matrix levels
!   according to flags in *FLG section. A separate module is used to read
!   tabulation axes for the TAB flag.
!
! VARIABLE KINDS
    USE KIND_DAT
! 
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE TANCOM_DAT ! Tangent path data
    USE ATMCOM_DAT, ONLY: HGTSFC, HGTTOA ! Hgts of surface and top of atmos.
    USE FLGCOM_DAT, ONLY: MTXFLG ! T = inter-level matrix of flux calculations
!
! SUBROUTINES
    USE ATMLEV_SUB ! Find/insert atmospheric level for given altitude
    USE C9REAL_GEN ! Write real number as C*9 string
    USE LOCASE_FNC ! Convert text string to lower case
    USE NXTFFL_SUB ! Load next field from rfm.drv, expanding filenames
    USE TANFLD_SUB ! Check if string is valid *TAN entry and insert in TANCOM
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for Driver File
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)       :: IATM   ! Atmos profile level
    INTEGER(I4)       :: ITAN   ! Counter for height levels
    INTEGER(I4)       :: JTAN   ! Index of secondary ray paths
    INTEGER(I4)       :: LENGTH ! Length of field read from driver table
    REAL(R4)          :: VALUE  ! numerical value read from FIELD
    CHARACTER(LENREC) :: FIELD  ! Tangent height or Tan.Hgt filename
    CHARACTER(3)      :: STRBOA ! BOA height as string
    CHARACTER(3)      :: STRTOA ! TOA height as string
    TYPE(TANTYP), ALLOCATABLE :: TANSAV(:) ! Saved version of TAN
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  CALL WRTLOG ( 'I-DRVHGT: ', .TRUE. ) 
!
! Read each field in *HGT section
  DO
    CALL NXTFFL ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
    IF ( LOCASE ( FIELD ) .EQ. 'boa' ) THEN
      FIELD = TRIM ( C9REAL(HGTSFC) )
      STRBOA = FIELD(1:3)
    ELSE IF ( LOCASE ( FIELD ) .EQ. 'toa' ) THEN
      FIELD = TRIM ( C9REAL(HGTTOA) )
      STRTOA = FIELD(1:3)
    END IF          
    CALL TANFLD ( FIELD, VALUE, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    IF ( VALUE .LT. HGTSFC ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-DRVHGT: Level = ' // TRIM ( FIELD(1:20) ) // &
         ' km below base of atmosphere (' // TRIM ( C9REAL(HGTSFC) ) // ' km)'
      RETURN
    ELSE IF ( VALUE .GT. HGTTOA ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-DRVHGT: Level = ' // TRIM ( FIELD(1:20) ) // &
         ' km above top of atmosphere (' // TRIM ( C9REAL(HGTTOA) ) //  ' km)'
      RETURN
    END IF
    CALL WRTLOG ( ' '//FIELD, .TRUE. ) 
  END DO
  CALL WRTLOG ( '', .FALSE. )
  IF ( STRBOA .NE. '' ) &
    CALL WRTLOG ( 'I-DRVTAN: Interpreted BOA as ' // STRBOA // ' km' )
  IF ( STRTOA .NE. '' ) &
    CALL WRTLOG ( 'I-DRVTAN: Interpreted TOA as ' // STRTOA // ' km' )
!
  IF ( NTAN .EQ. 0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVHGT: No entries in *HGT section'
    RETURN
  ELSE IF ( MTXFLG .AND. NTAN .EQ. 1 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVHGT: MTX flag requires at least 2 height levels'
    RETURN
  END IF
!
  TAN%HGT = TAN%USR
! Ensure each output level corresponds to an ATM profile level
! Since this is bottom upwards, %IAT levels won't be changed if further
! levels are interpolated. But assign %IAT via IATM to avoid accidental
! reassignment within atmlev.
  DO ITAN = 1, NTAN
    CALL ATMLEV ( TAN(ITAN)%HGT, .TRUE., IATM ) 
    TAN(ITAN)%IAT = IATM
  END DO
!
  IF ( MTXFLG ) THEN  ! extend TAN and replicate 1:NTAN
    MTAN = NTAN + NTAN**2
    CALL MOVE_ALLOC ( TAN, TANSAV ) 
    ALLOCATE ( TAN(MTAN) )
    TAN(1:NTAN) = TANSAV
    JTAN = 0
    DO ITAN = 1, NTAN
      JTAN = JTAN + NTAN
      TAN(JTAN+1:JTAN+NTAN) = TANSAV
      TAN(JTAN+1:JTAN+NTAN)%STR = TAN(ITAN)%STR
    END DO
  END IF
!
END SUBROUTINE DRVHGT
END MODULE DRVHGT_SUB

