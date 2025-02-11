PROGRAM RFM
!
! VERSION (update VIDHDR)
!   23SEP24 AD v5.21_23SEP Bug#48, modified TANCNV, CHKLIM
!   04JUL24 AD v5.21_24JUL Bug#46, modified TABPTH, TABWRT, ICLPTH
!   21JUN24 AD v5.21_21JUN Change CIACHK to warning messages
!   15APR24 AD v5.21_15APR Bug#45, modified DRVGAS
!   16FEB24 AD v5.21_16FEB Bug#44, modified NXTPRF
!   18DEC23 AD v5.21_18DEC Bug#42, modified JACPTH; Bug#43 RFMSPC, SPCFIN
!   24AUG23 AD v5.21_24AUG Bug#41, modified ATMINI
!   11AUG23 AD v5.21_11AUG Multiple HITRAN files. Remove DRVNAM.
!   30JUN23 AD v5.20_30JUN Bug#40, modified CHKABS.
!   22MAY23 AD v5.20_22MAY Bug#39, modified INIHFL.
!   20APR23 AD v5.20_20APR New HITRAN 20 molecules and TIPS data
!   15MAR23 AD v5.13_15MAR Change from MT_CKD v3.2 to v4.1 continuum
!   13MAR23 AD v5.13_13MAR Bug#38
!   19JAN23 AD v5.13_19JAN Bug#37
!   27SEP22 AD v5.12_27SEP Bugs#32-36
!   07MAR22 AD v5.12_07MAR Bug#31
!   03DEC21 AD v5.12_03DEC Speed up Jacobians
!   29APR21 AD v5.11_29APR Bug#29, 30.
!   02APR21 AD v5.11_02APR Bug#27 - new ATMGRD
!   19MAR21 AD v5.10_29MAR Bug#28 - new SFCLEV
!   25MAR21 AD v5.10_25MAR Extend handling of incl/excl profiles in *ATM files
!   19MAR21 AD v5.10_19MAR Add messages to log file from *ATM section
!   23MAY20 AD v5.10_23MAY Avoid NaN for low wno. in Planck fn.
!   25APR20 AD v5.10_25APR Fix bug in HITTYP
!   09MAR20 AD v5.10_09MAR Bug#26, revised MOLIDX
!   03MAR20 AD v5.10_03MAR Modified SPCRAD for NADir viewing OPT,TRA calc.
!   19FEB20 AD v5.10_19FEB Bug#24,#25, revised .xsc file handling
!   29JAN20 AD v5.10_29JAN Revised HITRAN file handling
!   20SEP19 AD v5.10_20SEP Bug#23
!   16AUG19 AD v5.10_16AUG Allow for .atm file qualifiers
!   05AUG19 AD v5.10_05AUG Allow for *OBS to specify pressure level 
!   01JUL19 AD v5.10_01JUL Bug#22
!   24JUN19 AD v5.10_24JUN Revised CIA data handling
!   07JUN19 AD v5.10_07JUN Bug#21
!   03JUN19 AD v5.10_03JUN Extended isolst
!   04APR19 AD v5.10_04APR Bug#20
!   28MAR19 AD v5.10_28MAR Bug#19
!   25MAR19 AD v5.10_25MAR - allows 'UNITS= ...' in *TAN/*LEN section
!   12MAR19 AD v5.10_12MAR - allow wavenumber for refractivity. Bug#18
!   05MAR19 AD v5.10_05MAR 
!   01FEB19 AD v5.02
!   01JUN18 AD v5.01
!   29JAN18 AD v5.00 F90 version. Tested.
!
! DESCRIPTION
!   Reference Forward Model
!   This is the F90 version.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE RUN_ID_DAT ! RFM Run ID string
    USE SPCCOM_DAT ! Spectral range data
    USE HDRCOM_DAT, ONLY: VIDHDR ! RFM version identifier
    USE RFMLUN_DAT, ONLY: LUNLOG ! LUN for log file
!
! SUBROUTINES
    USE RFMDAL_SUB ! Deallocate program level pointers
    USE RFMDRV_SUB ! Read RFM driver table
    USE RFMPRF_SUB ! Write out RFM internal atmospheric profile
    USE RFMPTH_SUB ! Set up RFM path calculations
    USE RFMSPC_SUB ! Calculations for one spectral range
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! LOCAL CONSTANTS
    LOGICAL, PARAMETER :: PROMPT = .FALSE. ! T=prompt user for appended ID
!
! LOCAL VARIABLES
    LOGICAL       :: FAIL   ! Set TRUE if a fatal error is detected
    INTEGER(I4)   :: IOS    ! Value of IOSTAT on OPEN  
    INTEGER(I4)   :: ISPC   ! Counter for spectral ranges
    CHARACTER(80) :: ERRMSG ! Error message if FAIL is TRUE
    CHARACTER(80) :: LOGMSG ! Text message sent to log file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  VIDHDR = '5.21_23SEP'
  LOGMSG = 'R-RFM: Running RFM v' // VIDHDR    
  WRITE ( *, '(A)' ) LOGMSG
!
  IF ( PROMPT ) THEN
    WRITE ( *, '(A)' ) 'Optional ID to be appended to filenames (<CR>=none):'
    READ ( *, '(A)' ) RUN_ID
    IF ( RUN_ID .NE. '' ) &
      WRITE ( *, '(A)' ) 'R-RFM:    Filename append string=' // RUN_ID
  END IF
!
! Open log file
  OPEN ( UNIT=LUNLOG, FILE='rfm.log'//RUN_ID, ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=IOS )
  IF ( IOS .NE. 0 ) THEN
    WRITE ( *, '(A)' ) 'F-RFM: Error opening rfm.log file. IOSTAT=', IOS
    STOP
  END IF
  CALL WRTLOG ( LOGMSG ) 
!
! Read driver table contents
  CALL RFMDRV ( FAIL, ERRMSG ) 
  IF ( FAIL ) GOTO 900
!
! If PRF flag, output atmospheric profile being used
  IF ( PRFFLG ) THEN 
    CALL RFMPRF ( FAIL, ERRMSG )
    IF ( FAIL ) GOTO 900
  END IF
!
! Calculate equivalent CG paths
  CALL RFMPTH ( FAIL, ERRMSG )
  IF ( FAIL ) GOTO 900
!
! Loop over required spectral ranges
  DO ISPC = 1, NSPC
    IF ( NSPC .GT. 1 ) THEN 
      LOGMSG = 'I-RFM: Calculation for spectral range: ' // SPC(ISPC)%LAB
      IF ( .NOT. SHHFLG ) WRITE ( *, '(A)' ) LOGMSG
      CALL WRTLOG ( LOGMSG ) 
    END IF
    CALL RFMSPC ( ISPC, FAIL, ERRMSG ) 
    IF ( FAIL ) GOTO 900
  END DO
!
! Deallocate pointers
  CALL RFMDAL 
!
900 CONTINUE
  IF ( FAIL ) THEN
    LOGMSG = ERRMSG
  ELSE
    LOGMSG = 'R-RFM: Successful completion'
  END IF
  CALL WRTLOG ( LOGMSG )
  WRITE ( *, '(A)' ) LOGMSG
!
END PROGRAM RFM
