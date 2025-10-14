MODULE HDBCHK_SUB
CONTAINS
SUBROUTINE HDBCHK ( LUNHIT, HDB, FAIL, ERRMSG )
!
! VERSION
!   07AUG24 AD Checked.
!   11AUG23 AD Pass HDB structure as argument.
!   31MAY23 AD Original.
!
! DESCRIPTION
!   Check HITRAN database file header record
!   Called by HDBFIL
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE HDBCOM_DAT ! HDB data structure
    USE SETCOM_DAT ! HITRAN line parameter sets
!
! SUBROUTINES
    USE IDXPAR_FNC ! Return IPAR value associated with HITRAN line parameter
    USE IDXSET_FNC ! Return index of line parameter set
    USE LENREC_FNC ! Determine length of next record in file
    USE LOCASE_FNC ! Convert text string to lower case
    USE TXTFLD_SUB ! Identify start and end points of text field
    USE TXTPAR_FNC ! Return text string describing HITRAN database line parameter
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNHIT ! Name of HITRAN file
    TYPE(HDBTYP),  INTENT(OUT) :: HDB    ! HDB file aux.data 
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: MAXFLD = 50 ! Max.no fields from database file
!
! LOCAL VARIABLES
    INTEGER(I4)    :: IEND   ! Location of end of field in header record
    INTEGER(I4)    :: IFLD   ! Counter for database fields
    INTEGER(I4)    :: IOS    ! Saved value of IOSTAT for error messages
    INTEGER(I4)    :: IPAR   ! Index of HITRAN parameter
    INTEGER(I4)    :: ISTA   ! Location of start of field in header record
    INTEGER(I4)    :: JEND   ! Location of end of field in data record
    INTEGER(I4)    :: JFLD   ! Counter for usable database fields 
    INTEGER(I4)    :: JSTA   ! Location of start of field in data record
    INTEGER(I4)    :: LENHDR ! Length of header record 
    INTEGER(I4)    :: NFLD   ! Counter for number of usable fields
    INTEGER(I4)    :: IPRFLD(MAXFLD)    ! Index of field typ
    INTEGER(I4)    :: IPTFLD(2,MAXFLD)  ! Start,End pos of field in record
    LOGICAL        :: GETPAR(MAXPAR)    ! T=get this parameter from file
    CHARACTER(MAXLEN) :: FIELD          ! Name of field extracted from header rec
    CHARACTER(:), ALLOCATABLE :: DATREC ! Data record from file
    CHARACTER(:), ALLOCATABLE :: HDRREC ! Header record from file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Read header record
  LENHDR = LENREC ( LUNHIT ) 
  ALLOCATE ( CHARACTER(LENHDR) :: HDRREC )
  READ ( LUNHIT, '(A)', IOSTAT=IOS, ERR=900 ) HDRREC
!
! Read first data record
  HDB%LEN = LENREC ( LUNHIT )
  ALLOCATE ( CHARACTER(HDB%LEN) :: DATREC ) 
  READ ( LUNHIT, '(A)', IOSTAT=IOS, ERR=900 ) DATREC
! 
! Go through fields
  GETPAR = .FALSE.
  IFLD = 0
  NFLD = NRQPAR
! Certain data fields might have spaces at start (or in middle - tbd)
! so assume field starts at 2 places beyond previous end
  JEND = -1
  DO 
    IFLD = IFLD + 1
    CALL TXTFLD ( HDRREC, IFLD, ISTA, IEND )
    IF ( IEND .GT. ISTA + MAXLEN-1 ) STOP 'F-HDBCHK: Logical error'
    IF ( ISTA .EQ. 0 ) EXIT
    FIELD = LOCASE ( HDRREC(ISTA:IEND) )
    JSTA = JEND + 2
! A couple of fields containg spaces, so set JEND according to actual size
    IF ( FIELD(1:6) .EQ. 'statep' ) THEN   ! 256 characters
      JEND = JSTA + 255
    ELSE IF ( FIELD .EQ. 'par_line' ) THEN ! 160 characters
      JEND = JSTA + 159
    ELSE
      CALL TXTFLD ( DATREC, IFLD, ISTA, JEND )  ! ISTA is dummy here
    END IF
!
! Test for a recognised/usable parameter
    IPAR = IDXPAR ( FIELD ) 
    IF ( IPAR .EQ. 0 ) CYCLE
    GETPAR(IPAR) = .TRUE.
    IF ( IPAR .LE. NRQPAR ) THEN ! Match Field# 1:NRQPAR with IPAR_ index
      JFLD = IPAR 
    ELSE                         ! Optional fields, arbitrary order
      NFLD = NFLD + 1
      JFLD = NFLD
    END IF
    IPRFLD(JFLD) = IPAR
    IPTFLD(:,JFLD) = (/ JSTA, JEND /)

! Create format string for reading wavenumber
    IF ( IPAR .EQ. IPAR_WNO ) THEN
      IF ( JSTA .EQ. 1 ) THEN 
        HDB%WNOFMT = '(F12.6)'
      ELSE
        WRITE ( HDB%WNOFMT, '(A,I4,A)' ) '(', JSTA-1, 'X,F12.6)'
      END IF
! Create format string for reading molecule ID
    ELSE IF ( IPAR .EQ. IPAR_IDM ) THEN
      IF ( JSTA .EQ. 1 ) THEN
        HDB%IDMFMT = '(I2)' 
      ELSE
        WRITE ( HDB%IDMFMT, '(A,I4,A)' ) '(', JSTA-1, 'X,I2)'
      END IF
    END IF
  END DO
  HDB%NFD = NFLD
  ALLOCATE ( HDB%IPR(NFLD), HDB%IPT(NFLD), HDB%JPT(NFLD) )
  HDB%IPR = IPRFLD(1:NFLD)
  HDB%IPT = IPTFLD(1,1:NFLD)
  HDB%JPT = IPTFLD(2,1:NFLD)
!
  CALL WRTLOG ( 'I-HDBCHK: File contains additional line parameters: ', &
                .TRUE. )
!
! Check all required fields are present
  DO IPAR = 1, MAXPAR
    IF ( GETPAR(IPAR) ) THEN 
      IF ( IPAR .GT. NRQPAR ) CALL WRTLOG ( TXTPAR(IPAR), .TRUE. ) 
    ELSE
      IF ( IPAR .LE. NRQPAR ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-HDBCHK: HITRAN data missing required field: ' // &
                 TRIM ( TXTPAR(IPAR) )
        RETURN
      END IF
    END IF  
  END DO
!
! Set flags in HITCOM showing which optional parameters might be set
  HDB%IST = IDXSET ( GETPAR ) 
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-HDBCHK: Read failure on HITRAN file. IOSTAT=', IOS
!
END SUBROUTINE HDBCHK
END MODULE HDBCHK_SUB
