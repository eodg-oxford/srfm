MODULE HITFIL_SUB
CONTAINS
SUBROUTINE HITFIL ( NAMHIT, TYPHIT, WNOREQ, IDMREQ, FAIL, ERRMSG )
!
! VERSION
!   13AUG24 AD Checked.
!   11AUG23 AD Adpapted for multiple files
!   30MAY23 AD Renamed from OPNHIT and simplified structure.
!   29JAN20 AD Open all forms of HITRAN data file. Checked.
!   12JUN17 AD Allow for ASCII files. Checked.
!   O1MAY17 AD F90 conversion. 
!
! DESCRIPTION
!   Open HITRAN line data file and check contents
!   Called by DRVHIT once for each file
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE HFLCOM_DAT ! HITRAN file data
    USE HDBCOM_DAT ! HDB data structure
    USE RFMLUN_DAT, ONLY: LUNNXT ! Next free LUN
!
! SUBROUTINES
    USE BINFIL_SUB ! Check HITRAN binary data file
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE C9REAL_GEN ! Write real number as C*9 string
    USE HDBFIL_SUB ! Check HITRAN line database file
    USE HITQAL_SUB ! Extract qualifiers from .hit filename
    USE IDXSET_FNC ! Return index of line parameter set
    USE LEXIST_FNC ! Check if file exists
    USE PARFIL_SUB ! Check HITRAN line data .par file
    USE TWOQAL_SUB ! Extract pair of numbers from '( : )' string
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMHIT ! Name of HITRAN file
    CHARACTER(6),  INTENT(IN)  :: TYPHIT ! Type of file, eg 'BINFIL'
    REAL(R8),      INTENT(IN)  :: WNOREQ(2) ! Required wavenumber range [cm-1]
    LOGICAL,       INTENT(IN)  :: IDMREQ(:) ! T=use molec# from file
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if fatal error detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL      :: FIXWNO    ! T=User-specified wavenumber range
    LOGICAL      :: LEXCLD    ! T=list of excluded molecules for file
    LOGICAL      :: USEFIL    ! T=file contains usable contents
    INTEGER(I4)  :: IOS       ! Saved value of IOSTAT for error messages
    TYPE(HDBTYP) :: HDB       ! HITRAN database file auxiliary data
    REAL(R8)     :: WNOFIL(2) ! Specified Wavenumber range [cm-1] for file
    LOGICAL      :: IDMFIL(SIZE(IDMREQ)) ! T=use molec# from file
    LOGICAL      :: IDMSUB(SIZE(IDMREQ)) ! Subset of reqd mols for this file
    LOGICAL      :: QALIDM(SIZE(IDMREQ)) ! T=molec# specified in file qualifier
    CHARACTER(LEN(NAMHIT))    :: NAMFIL  ! Local version of NAMHIT
    TYPE(HFLTYP), ALLOCATABLE :: HFLSAV(:) ! Saved HFL during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NAMFIL = NAMHIT
!
! Check for wavenumber range specified as filename qualifier
  WNOFIL = WNOREQ
  CALL TWOQAL ( NAMFIL, .TRUE., FIXWNO, WNOFIL, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN 
  IF ( FIXWNO ) THEN
    CALL WRTLOG ( 'I-HITFIL: Specified wavenumber range for file: ' // & 
      TRIM(C9REAL(WNOFIL(1))) // '-' // TRIM(C9REAL(WNOFIL(2))) // ' [cm-1]' )
    IF ( WNOFIL(1) .GT. WNOREQ(1) .OR. WNOFIL(2) .LT. WNOREQ(2) ) THEN
      CALL WRTLOG ( 'W-HITFIL: Ignore file, outside required wavenumber range' )
      RETURN
    END IF
    WNOFIL = 0.0
  END IF
!
  CALL HITQAL ( NAMFIL, QALIDM, LEXCLD, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN
  IDMSUB = IDMREQ
  IF ( ANY ( QALIDM ) ) THEN 
    IF ( LEXCLD ) THEN 
      WHERE ( QALIDM ) IDMSUB = .FALSE.
    ELSE
      IDMFIL = QALIDM
      IDMSUB = .FALSE.
    END IF
  END IF        
!
  CALL WRTLOG ( 'I-HITFIL: Opening file#' // TRIM(C11INT(NHFL+1)) // ' ' &
                // TRIM(NAMFIL) ) 
  IF ( .NOT. LEXIST ( NAMFIL ) ) THEN 
    FAIL = .TRUE.
    ERRMSG = 'F-HITTYP: file not found'
    RETURN
  END IF
!
  SELECT CASE ( TYPHIT )
    CASE ( 'BINFIL' ) 
      CALL BINFIL ( LUNNXT, NAMFIL, WNOFIL, IDMSUB, &
                    USEFIL, IDMFIL, FAIL, ERRMSG )
    CASE ( 'HDBFIL' ) 
      CALL HDBFIL ( LUNNXT, NAMFIL, WNOFIL, IDMSUB, &
                    USEFIL, IDMFIL, HDB, FAIL, ERRMSG )
    CASE ( 'PARFIL' ) 
      CALL PARFIL ( LUNNXT, NAMFIL, WNOFIL, IDMSUB, &
                    USEFIL, IDMFIL, FAIL, ERRMSG )
    CASE DEFAULT
      STOP 'F-HITFIL: Logical error'
  END SELECT
  IF ( FAIL ) RETURN

! Add new HITRAN file
  IF ( USEFIL ) THEN 
    IF ( NHFL .GT. 0 ) CALL MOVE_ALLOC ( HFL, HFLSAV )
    NHFL = NHFL + 1
    ALLOCATE ( HFL(NHFL) ) 
    IF ( NHFL .GT. 1 ) HFL(1:NHFL-1) = HFLSAV
    HFL(NHFL)%LUN = LUNNXT
    HFL(NHFL)%TYP = TYPHIT(1:3)
    HFL(NHFL)%IDM = IDMFIL
    HFL(NHFL)%HIT%IUV = 0  ! Define as zero, only changed if non-LTE
    HFL(NHFL)%HIT%ILV = 0
    IF ( TYPHIT .EQ. 'HDBFIL' ) THEN
      HFL(NHFL)%HDB = HDB
      HFL(NHFL)%HIT%IST = HDB%IST  ! Assign set of possible line parameters
    ELSE
      HFL(NHFL)%HIT%IST = IDXSET() ! Assign default set of reqd.params only
    END IF
    LUNNXT = LUNNXT + 1
  ELSE
    CLOSE ( LUNNXT, IOSTAT=IOS, ERR=900 ) 
  END IF
!
  RETURN
!
900 CONTINUE
  FAIL = .TRUE.
  WRITE ( ERRMSG, * ) &
    'F-HITFIL: Error closing HITRAN file. IOSTAT=', IOS
!
END SUBROUTINE HITFIL
END MODULE HITFIL_SUB
