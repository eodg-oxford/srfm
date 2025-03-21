MODULE RECHDB_SUB
CONTAINS
SUBROUTINE RECHDB ( LUNHIT, HIT, USEIDM, WNOMAX, HDB, EOF, FAIL, ERRMSG ) 
!
! VERSION
!   24AUG24 AD Checked.
!   11AUG23 AD Pass data via arguments
!   31MAY23 AD Original. 
!
! DESCRIPTION
!   Read record from HITRAN database file
!   Called by HITREC and INIHIT.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE HDBCOM_DAT ! HDB data structure
    USE SETCOM_DAT ! HITRAN line parameter sets
    USE HITCOM_DAT, ONLY: HITTYP ! HITRAN line data structure
    USE IDXCON_DAT, ONLY: IDXH2O ! RFM/HITRAN index for H2O
    USE PHYCON_DAT, ONLY: AVOG   ! Avogradro's number [kmol/cm2]
!
  IMPLICIT NONE
!
! ARGUMENTS      
    INTEGER(I4),   INTENT(IN)  :: LUNHIT    ! Logical unit number of HITRAN file
    TYPE(HITTYP),  INTENT(OUT) :: HIT       ! HITRAN Line parameters
    LOGICAL,       INTENT(IN)  :: USEIDM(:) ! T=use molecule#
    REAL(R8),      INTENT(IN)  :: WNOMAX    ! Highest wno [cm-1] reqd
    TYPE(HDBTYP),  INTENT(IN)  :: HDB       ! HITRAN Database structure info
    LOGICAL,       INTENT(OUT) :: EOF       ! T = reached end of file/spc.range
    LOGICAL,       INTENT(OUT) :: FAIL      ! T = fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG    ! Error message written if FAIL=T 
!
! LOCAL VARIABLES
    INTEGER(I4)        :: IDM     ! HITRAN Index of molecule
    INTEGER(I4)        :: IFLD    ! Counter for usable fields in data record
    INTEGER(I4)        :: IOS     ! Saved value of IOSTAT for error messages
    INTEGER(I4)        :: IPAR    ! Index of HITRAN parameter type
    INTEGER(I4)        :: IPT     ! Position of start of field in record
    INTEGER(I4)        :: JPT     ! Position of end of field in record
    REAL(R8)           :: DSTR    ! Line strength allowing for < 1.0E-38
    REAL(R8)           :: WNO     ! Wavenumber [cm-1] 
    CHARACTER(HDB%LEN) :: RECORD  ! Input record read as text
    CHARACTER(MAXLEN)  :: SUBSTR  ! Data field from record
    LOGICAL     :: GOTPAR(MAXPAR) ! T=valid param found in record
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
  EOF = .FALSE.
!
! Cycle until a record with all required parameters is found
  DO 
    READ ( LUNHIT, '(A)', IOSTAT=IOS, ERR=900, END=800 ) RECORD
    READ ( RECORD, HDB%WNOFMT, IOSTAT=IOS, ERR=900 ) WNO
    IF ( WNO .GT. WNOMAX ) GOTO 800
    READ ( RECORD, HDB%IDMFMT, IOSTAT=IOS, ERR=900 ) IDM
    IF ( IDM .LE. 0 .OR. IDM .GT. SIZE(USEIDM) ) CYCLE
    IF ( .NOT. USEIDM(IDM) ) CYCLE
!
    GOTPAR = .FALSE.
    DO IFLD = 1, HDB%NFD
      IPT = HDB%IPT(IFLD)
      JPT = HDB%JPT(IFLD)
      SUBSTR = RECORD(IPT:JPT)
      IF ( SUBSTR(1:1) .NE. '#' ) THEN
        IPAR = HDB%IPR(IFLD)
        GOTPAR(IPAR) = .TRUE.
        SELECT CASE ( IPAR )
          CASE ( IPAR_WNO ) ; HIT%WNO = WNO
          CASE ( IPAR_IDM ) ; HIT%IDM = IDM
          CASE ( IPAR_IDI ) 
            READ ( SUBSTR, *, IOSTAT=IOS, ERR=900 ) HIT%IDI
! CO2 Isotope#10 returned as 0, although Iso#11 and 12 OK
            IF ( HIT%IDI .EQ. 0 ) HIT%IDI = 10
          CASE ( IPAR_STR ) 
            READ ( SUBSTR, *, IOSTAT=IOS, ERR=900 ) DSTR
            HIT%STR = SNGL ( DSTR * AVOG )
          CASE ( IPAR_ELS ) ; READ ( SUBSTR, *, IOSTAT=IOS, ERR=900 ) HIT%ELS
          CASE ( IPAR_HWA ) ; READ ( SUBSTR, *, IOSTAT=IOS, ERR=900 ) HIT%HWA
          CASE ( IPAR_HWS ) ; READ ( SUBSTR, *, IOSTAT=IOS, ERR=900 ) HIT%HWS
          CASE ( IPAR_TCA ) ; READ ( SUBSTR, *, IOSTAT=IOS, ERR=900 ) HIT%TCA
          CASE ( IPAR_TCS ) ; READ ( SUBSTR, *, IOSTAT=IOS, ERR=900 ) HIT%TCS
          CASE ( IPAR_PSA ) ; READ ( SUBSTR, *, IOSTAT=IOS, ERR=900 ) HIT%PSA
          CASE ( IPAR_PSS ) ; READ ( SUBSTR, *, IOSTAT=IOS, ERR=900 ) HIT%PSS
          CASE ( IPAR_LMA ) ; READ ( SUBSTR, *, IOSTAT=IOS, ERR=900 ) HIT%LMA
          CASE ( IPAR_LMS ) ; READ ( SUBSTR, *, IOSTAT=IOS, ERR=900 ) HIT%LMS
          CASE DEFAULT ; STOP 'F-RECHDB: Logical error'
        END SELECT
      END IF
    END DO
! If all required parameters have been found, set optional fields & exit
    IF ( ALL ( GOTPAR(1:NRQPAR) ) ) THEN
! Set optional fields whether specified for this file or not
      IF ( .NOT. GOTPAR(IPAR_TCS) ) HIT%TCS = HIT%TCA
      IF ( .NOT. GOTPAR(IPAR_PSA) ) HIT%PSA = 0.0
      IF ( .NOT. GOTPAR(IPAR_PSS) ) HIT%PSS = HIT%PSA
      IF ( .NOT. GOTPAR(IPAR_LMA) ) HIT%LMA = 0.0
      IF ( .NOT. GOTPAR(IPAR_LMS) ) HIT%LMS = HIT%LMA
      RETURN    ! Exit with HITCOM variables set
    END IF
!
  END DO
!
! Exit with EOF or past WNOMAX
800 CONTINUE
  EOF = .TRUE. 
  RETURN
!
900 CONTINUE
  FAIL = .TRUE.
  WRITE ( ERRMSG, * ) &
    'F-RECHDB: Failed to read HITRAN database file. IOSTAT=', IOS
!
END SUBROUTINE RECHDB
END MODULE RECHDB_SUB

