MODULE RFMPRF_SUB
CONTAINS
SUBROUTINE RFMPRF ( FAIL, ERRMSG )
!
! VERSION
!   01APR23 AD Checked.
!   08NOV17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Write out RFM internal atmospheric profile
!   Called once by RFM if PRF flag selected.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE FLGCOM_DAT ! Option flags
    USE HDRCOM_DAT ! Output header data
    USE NAMCOM_DAT ! RFM output filenames
    USE NTECOM_DAT ! Non-LTE data
    USE GRACOM_DAT, ONLY: NPSI, NOMGRA   ! No. horiz prof locs, index of psi=0
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE MAKNAM_SUB ! Construct filename for RFM output files
    USE MOVGRA_SUB ! Move profiles within GRACOM arrays
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL CONSTANTS
    CHARACTER(10), PARAMETER :: FMT = '(1P5G13.5)' ! Format for output values
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IOS    ! Value of IOSTAT saved for error messages.
    INTEGER(I4)   :: IPSI   ! Counter for horizontal profile locations
    INTEGER(I4)   :: IVIB   ! Counter for vibrational temperature profiles
    INTEGER(I4)   :: IVMR   ! Counter for VMR profiles
    CHARACTER(7)  :: STASTR = 'UNKNOWN' ! Status for OPEN statements
    CHARACTER(LENNAM) :: FILNAM ! Name of file actually opened (incl. RUNID)
!
! EXECUTABLE CODE -------------------------------------------------------------
! 
  FAIL = .FALSE.
  IF ( .NOT. PRFFLG ) RETURN   ! Only output if PRF flag enabled
!
  IF ( NEWFLG ) STASTR = 'NEW'
!
  DO IPSI = 1, MAX ( NPSI, 1 )
!
! Construct filename and open file
    IF ( GRAFLG ) THEN       ! Insert angle into output filename
      CALL MAKNAM ( PRFNAM, FILNAM, IPSI=IPSI )
      CALL MOVGRA ( IPSI, 0 )
    ELSE                     ! No angle required
      CALL MAKNAM ( PRFNAM, FILNAM ) 
    END IF
    CALL WRTLOG ( 'I-RFMPRF: Opening output file: ' // FILNAM )
!
    OPEN ( UNIT=LUNTMP, FILE=FILNAM, STATUS=STASTR, ACTION='WRITE', &
           IOSTAT=IOS, ERR=900 )
!
! Write File Header 
    WRITE ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) &
      '! RFM internal profile written by RFM v.' // VIDHDR
    WRITE ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) TXTHDR
!        
    WRITE ( LUNTMP, '(I6,A)', IOSTAT=IOS, ERR=900 ) &
      NATM, ' ! No. profile levels'
    WRITE ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) '*HGT [km]'
    WRITE ( LUNTMP, FMT, IOSTAT=IOS, ERR=900 ) HGTATM
!
    WRITE ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) '*PRE [mb]'
    WRITE ( LUNTMP, FMT, IOSTAT=IOS, ERR=900 ) PREATM
!
    WRITE ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) '*TEM [K]'
    WRITE ( LUNTMP, FMT, IOSTAT=IOS, ERR=900 ) TEMATM
!
    DO IVMR = 1, NVMR
      IF ( IVMR .EQ. IAXVMR ) THEN
        WRITE ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) &
          '*' // NAMVMR(IVMR) // ' [km-1]'
        WRITE ( LUNTMP, FMT, IOSTAT=IOS, ERR=900 ) EXTATM
      ELSE
        WRITE ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) &
          '*' // NAMVMR(IVMR) // ' [ppmv]'
        WRITE ( LUNTMP, FMT, IOSTAT=IOS, ERR=900 ) VMRATM(:,IVMR)
      END IF
    END DO
!
    DO IVIB = 1, NVIB
      WRITE ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) & 
        '*' // NTE(IVIB)%COD // ' [K] (Vib-Kin Temp)' 
      WRITE ( LUNTMP, FMT, IOSTAT=IOS, ERR=900 ) VIBATM(:,IVIB)
    END DO
!
    WRITE ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) '*END'
    CLOSE ( LUNTMP, IOSTAT=IOS, ERR=900 )
!
  END DO
!
  IF ( GRAFLG ) CALL MOVGRA ( NOMGRA, 0 )          
!
  900 CONTINUE
  FAIL = IOS .NE. 0
  IF ( FAIL ) WRITE ( ERRMSG, * ) & 
     'F-RFMPRF: I/O failure on output file. IOSTAT=', IOS
!
END SUBROUTINE RFMPRF
END MODULE RFMPRF_SUB

