MODULE CIAFIL_SUB
CONTAINS
SUBROUTINE CIAFIL ( NAMCIA, FAIL, ERRMSG )
!
! VERSION
!   29JUL24 AD Checked.
!   24JUN19 AD Add ADDAIR if 'Air' specified as broad.gas. Set GAS%CIA. Checked
!   01MAY17 AD F90 conversion. Tested
!
! DESCRIPTION
!   Read Collision-Induced Abs. data file
!   Called by DRVCIA or CIADEF for each .cia file 
!   Treats .cia files as arbitrary concatenations of datasets.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE SPCCOM_DAT ! Spectral range data
    USE IDXCON_DAT, ONLY: IDXAIR ! RFM 'molecular' index for air
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE ADDAIR_SUB ! Add air as absorber and atmospheric profile
    USE CIACHK_SUB ! Check CIA data consistent with other RFM options
    USE CIAMOL_SUB ! Identify molecules from .cia file header record
    USE OPNFIL_SUB ! Open input file
    USE REACIA_SUB ! Read data from .cia file
    USE USEMOL_FNC ! T = molecule is required
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMCIA ! Name of CIA data file
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    REAL(R8), PARAMETER :: WNOTOL = 1.0D0 ! Tol. [cm-1] for matching spc ranges
!
! LOCAL VARIABLES
    LOGICAL       :: NEWSPC    ! T=new spectral range within file data
    INTEGER(I4)   :: ID1,ID2   ! HITRAN/RFM indices of collision molecules
    INTEGER(I4)   :: IGAS      ! Index of molecule in GASCOM, or 0 if not reqd
    INTEGER(I4)   :: IOS       ! Saved value of IOSTAT for error message
    INTEGER(I4)   :: IPT       ! Counter for (wno,data) points
    INTEGER(I4)   :: NPT       ! No. data points in dataset
    INTEGER(I4)   :: NTEM      ! No. temperature tabulations for each spec.range
    REAL(R4)      :: TEM       ! Temperature [K] associated with CIA dataset
    REAL(R4)      :: TMAX,TMIN ! Max/Min Temp [K] for each spectral range
    REAL(R8)      :: WNL,WNU   ! Min/Max Wno [cm-1] read for dataset
    REAL(R8)      :: WNLSAV, WNUSAV ! Saved values of WNL, WNU
    CHARACTER(80) :: MESSGE    ! Message for Log file
    CHARACTER(20) :: MOLID     ! Part of .cia file hdr containing molec. names
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  WNLSAV = 0.0D0
  WNUSAV = 0.0D0
  NTEM = 0
!
! Open file and print first record to log file
  CALL OPNFIL ( LUNTMP, NAMCIA, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN
!
  DO 
    READ ( LUNTMP, '(A20,2F10.4,I7,F7.1,F10.4)', IOSTAT=IOS, ERR=900 ) &
      MOLID, WNL, WNU, NPT, TEM
    IF ( IOS .LT. 0 ) EXIT                ! end-of-file
    CALL CIAMOL ( ADJUSTL(MOLID), ID1, ID2, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
!
! If either molecule not specified for calc, then skip file.
! But allow for second molecule being 'air', which won't be explicitly specified
    IF ( .NOT. USEMOL(ID1) .OR. &
         .NOT. ( USEMOL(ID2) .OR. ID2 .EQ. IDXAIR ) ) EXIT
!
! Set GAS%CIA to indicate that some .cia data has been specified for this
! molecule, even if it does not match any required spectral range
    IGAS = IGSMOL(ID1)
    IF ( IGAS .GT. 0 ) GAS(IGAS)%CIA = .TRUE.
    IGAS = IGSMOL(ID2)
    IF ( IGAS .GT. 0 ) GAS(IGAS)%CIA = .TRUE.
!
! Broadener ID2 may be 'Air', which won't initially be listed as a required
! absorber, so USEMOL will return FALSE.
    IF ( ANY ( WNL .LE. SPC%WXU .AND. WNU .GE. SPC%WXL ) .AND. &
         USEMOL(ID1) .AND. ( USEMOL(ID2) .OR. ID2 .EQ. IDXAIR ) ) THEN
      CALL CIACHK ( ID1, ID2, TEM, WNL, WNU, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
! If Air-broadening is required, now ensure Air is listed as an absorber 
      IF ( ID2 .EQ. IDXAIR ) CALL ADDAIR
!
      CALL REACIA ( LUNTMP, ID1, ID2, NPT, TEM, WNL, WNU, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      NEWSPC = ABS ( WNL - WNLSAV ) .GT. WNOTOL .OR. &
               ABS ( WNU - WNUSAV ) .GT. WNOTOL
      IF ( NEWSPC ) THEN 
        IF ( NTEM .GT. 0 ) CALL WRTLOG ( MESSGE ) ! output previous message
        WNLSAV = WNL
        WNUSAV = WNU
        TMIN = TEM
        TMAX = TEM
        NTEM = 1
        WRITE ( MESSGE, '(A,F11.4,A,F11.4,A,I4,A,F6.1,A)' ) &
          'I-CIAFIL:', WNLSAV, ' -', WNUSAV, ' cm-1 read', NTEM, &
          ' Temp at', TMIN, ' K'
      ELSE 
        TMIN = MIN ( TMIN, TEM ) 
        TMAX = MAX ( TMAX, TEM ) 
        NTEM = NTEM + 1
        WRITE ( MESSGE, '(A,F11.4,A,F11.4,A,I4,A,F6.1,A,F6.1,A)' ) &
          'I-CIAFIL:', WNLSAV, ' -', WNUSAV, ' cm-1 read', NTEM, &
          ' Temps from', TMIN, ' -', TMAX, ' K'
      END IF
    ELSE             ! skip dataset
      DO IPT = 1, NPT
        READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) ! WDUM, RDUM
      END DO
    END IF
  END DO
!
  IF ( NTEM .GT. 0 ) THEN
    CALL WRTLOG ( MESSGE ) ! output last message
  ELSE 
    CALL WRTLOG ( 'W-CIAFIL: File ignored - contains no required data' )
  END IF
!
  CLOSE ( LUNTMP, IOSTAT=IOS, ERR=900 )
!
900 CONTINUE
  FAIL = IOS .GT. 0
  IF ( FAIL ) WRITE ( ERRMSG, * ) & 
    'F-CIAFIL: Error reading CIA file. IOSTAT=', IOS
!
END SUBROUTINE CIAFIL
END MODULE CIAFIL_SUB
