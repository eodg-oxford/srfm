MODULE REAXSC_SUB
CONTAINS
SUBROUTINE REAXSC ( LUNXSC, RFMFMT, MOLEC1, WNLSPC, WNUSPC, &
                    NODATA, FAIL, ERRMSG )
!
! VERSION
!   19FEB20 AD Rewritten to load data into XFLCOM. Checked.
!   08FEB19 AD Bug#16: Extract no.triangles NTRI from TRIANG.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Read data from .xsc file
!   Called by XSCFIL for each file.
!   Loads data into XFLCOM, but only datasets which span WNOMIN:WNOMAX
!   The file is opened and closed outside this subroutine.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE XFLCOM_DAT ! Contents of .xsc file
!
! SUBROUTINES
    USE MOLIDX_SUB ! Give molecule name for HITRAN/RFM index, or vice-versa
!
  IMPLICIT NONE 
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)    :: LUNXSC    ! LUNXSC for accessing file
    LOGICAL,       INTENT(IN)    :: RFMFMT    ! T=RFM format, F=HITRAN format
    CHARACTER(20), INTENT(INOUT) :: MOLEC1    ! Molecule name from 1st dataset 
    REAL(R8),      INTENT(IN)    :: WNLSPC(:) ! Min.Wno required for RFM 
    REAL(R8),      INTENT(IN)    :: WNUSPC(:) ! Max.Wno required for RFM 
    LOGICAL,       INTENT(OUT)   :: NODATA    ! T=No data within reqd spec.range
    LOGICAL,       INTENT(OUT)   :: FAIL      ! TRUE if fatal error is detected
    CHARACTER(80), INTENT(OUT)   :: ERRMSG    ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IDX1   ! Molec. index corresponding to MOLEC1
    INTEGER(I4)   :: IDX2   ! Molec. index corresponding to MOLEC2
    INTEGER(I4)   :: IOF    ! Offset in ABCXFL for loading abs.coeff data
    INTEGER(I4)   :: IOS    ! Saved value of IOSTAT for error message
    INTEGER(I4)   :: IPT    ! Counter for reading tabulated data points
    INTEGER(I4)   :: IXFL   ! Counter for stored tabulations from file
    INTEGER(I4)   :: IXSP   ! Index of spectral band
    INTEGER(I4)   :: NPT    ! No. abs.coeff values within each tabulation
    INTEGER(I4)   :: NXSP   ! No. distinct spectral ranges in file
    INTEGER(I4)   :: NXST   ! No. tabulations for current spec.rng (RFM format)
    REAL(R4)      :: PRE    ! Pressure [Torr] for tabulation
    REAL(R4)      :: RDUM   ! Dummy variable for reading abs.coeff values
    REAL(R4)      :: TEM    ! Temperature [K] for tabulation
    REAL(R8)      :: WN1    ! Lower Wno [cm-1] for particular tabulation
    REAL(R8)      :: WN2    ! Upper Wno [cm-1] for particular tabulation
    CHARACTER(20) :: MOLEC2 ! Molecule name for each tabulation
    REAL(R4),     ALLOCATABLE :: ABCSAV(:) ! Saved ABCXFL during reallocation
    TYPE(XFLTYP), ALLOCATABLE :: XFLSAV(:) ! Saved XFL during reallocation
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  NXFL = 0
  IOF = 0
  NXSP = 0 
  NXST = 0
  IDX1 = 0
  CALL MOLIDX ( IDX1, MOLEC1 ) ! Get IDX1 corresponding to MOLEC1
!
! Loop over all tabulated datasets within .xsc file
  DO
    IF ( RFMFMT ) THEN  ! Skip RFM-specific header for NXST tabulations
      IF ( NXST .EQ. 0 ) &
        READ ( LUNXSC, '(3X,I7,I10)', IOSTAT=IOS, END=900, ERR=900 ) IDX2, NXST
      NXST = NXST - 1
    END IF
    READ ( LUNXSC, '(A20,2F10.4,I7,F7.2,F6.1)', IOSTAT=IOS, END=900, ERR=900 ) &
      MOLEC2, WN1, WN2, NPT, TEM, PRE
    MOLEC2 = ADJUSTL ( MOLEC2 ) 
!
! Compare Molec.index since same molecule may have different forms of name
! Assume this has already been done for RFM-format and molecule ID= IDX2
    IF ( .NOT. RFMFMT ) THEN
      IDX2 = 0
      CALL MOLIDX ( IDX2, MOLEC2 )
    END IF
!
    IF ( IDX1 .NE. IDX2 ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-REAXSC: Inconsistent Molecule names in file: ' // &
               TRIM ( MOLEC1 ) // '&' // TRIM ( MOLEC2 ) 
      RETURN
    END IF 
!
! Check if some overlap with min/max spectral range required for RFM 
! To save space, only store data if it overlaps required RFM range
    IF ( ALL ( WN1 .GE. WNUSPC .OR. WN2 .LE. WNLSPC ) ) THEN  ! skip 
      READ ( LUNXSC, '(10E10.3)', IOSTAT=IOS, ERR=900 ) ( RDUM, IPT = 1, NPT ) 
      CYCLE
    END IF
!
! Identify if overlapping any spectral bands already identified within data
    IXSP = 0
    DO IXFL = 1, NXFL
      IF ( WN1 .LT. XFL(IXFL)%WN2 .AND. WN2 .GT. XFL(IXFL)%WN1 ) THEN
        IF ( IXSP .EQ. 0 ) THEN          ! found first match
          IXSP = XFL(IXFL)%ISP
        ELSE IF ( IXSP .NE. XFL(IXFL)%ISP ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-REAXSC: Ambiguous assigment of spectral bands ' & 
                    // 'in .xsc file'
          RETURN
        END IF
      END IF
    END DO         
    IF ( IXSP .EQ. 0 ) THEN    ! Identified as new spectral band
      NXSP = NXSP + 1
      IXSP = NXSP
    END IF
!
! Save tabulation 
    IF ( NXFL .GT. 0 ) CALL MOVE_ALLOC ( XFL, XFLSAV )
    NXFL = NXFL + 1
    ALLOCATE ( XFL(NXFL) )
    IF ( NXFL .GT. 1 ) XFL(1:NXFL-1) = XFLSAV
    XFL(NXFL)%IOF = IOF
    XFL(NXFL)%ISP = IXSP
    XFL(NXFL)%NPT = NPT      
    XFL(NXFL)%PRE = PRE
    XFL(NXFL)%TEM = TEM
    XFL(NXFL)%WN1 = WN1
    XFL(NXFL)%WN2 = WN2
! Save absorption coefficient data
    IF ( NXFL .GT. 1 ) CALL MOVE_ALLOC ( ABCXFL, ABCSAV )
    ALLOCATE ( ABCXFL(IOF+NPT) ) 
    IF ( NXFL .GT. 1 ) ABCXFL(1:IOF) = ABCSAV
    READ ( LUNXSC, '(10E10.3)', IOSTAT=IOS, ERR=900 ) &
      ( ABCXFL(IPT), IPT = IOF+1, IOF+NPT ) 
    IOF = IOF + NPT
!
  END DO
!
! Jump here when reached end-of-file (with IOS < 0 )
900 CONTINUE
  NODATA = NXFL .EQ. 0
  FAIL = IOS .GT. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-REAXSC: Error reading X/S file. IOSTAT=', IOS
!
END SUBROUTINE REAXSC
END MODULE REAXSC_SUB
