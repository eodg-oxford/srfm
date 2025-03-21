MODULE LUTTAB_SUB
CONTAINS
SUBROUTINE LUTTAB ( LUNLUT, NAMLUT, BINFIL, ILUT, NDP, NDT, FAIL, ERRMSG )
!
! VERSION
!   30MAY24 AD Checked.
!   28MAR19 AD Bug#19: Also set LUT%WNL, LUT%WNU. Add GHZ2CM. Checked.
!   01JUL17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Initialise LUT-TAB file
!   Called by INILUT for each TAB LUT file.
!   This loads file header data into TABCOM and keeps the LUN for the file open.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LUTCOM_DAT ! TAB LUT data
    USE TABCOM_DAT ! Axes for tabulated absorption coefficients
    USE PHYCON_DAT, ONLY: GHZ2CM ! GHz to cm-1 conv.fac
!
! SUBROUTINES
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
    USE LUTPTH_SUB ! Set LUT path interpolation
    USE OPNFIL_SUB ! Open input file
    USE SUBAXS_SUB ! Subsample (p,T) axes from LUT TAB input file
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNLUT ! LUN for LUT file
    CHARACTER(*),  INTENT(IN)  :: NAMLUT ! Name of LUT file
    LOGICAL,       INTENT(IN)  :: BINFIL ! T = binary file, F =ASCII
    INTEGER(I4),   INTENT(IN)  :: ILUT   ! Index of LUT
    INTEGER(I4),   INTENT(IN)  :: NDP    ! p-axis subsample fac (1=no subsamp.)
    INTEGER(I4),   INTENT(IN)  :: NDT    ! T-axis subsample fac (1=no subsamp.)
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IDXISO ! Isotopic ID of gas, read from LUT file
    INTEGER(I4)   :: IDXMOL ! HITRAN ID of gas, read from LUT file 
    INTEGER(I4)   :: IGAS   ! Index of gas in LUT file
    INTEGER(I4)   :: IOS    ! Saved value of IOSTAT for error messages
    INTEGER(I4)   :: NLKP   ! No. p,T,q interpolation points required
    REAL(R4)      :: RGAS   ! Molec#ID plus Iso#ID read as a real number
    CHARACTER(5)  :: GASSTR ! Molec/Iso ID part of header record
    CHARACTER(80) :: RECORD ! Record read from LUT file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  LUT(ILUT)%LUN = LUNLUT
  LUT(ILUT)%BIN = BINFIL
!
  IF ( BINFIL ) THEN
    CALL WRTLOG ( 'I-LUTTAB: Opening file: ' // NAMLUT )
    OPEN ( UNIT=LUNLUT, FILE=NAMLUT, STATUS='OLD', ACTION='READ', &
           FORM='UNFORMATTED', IOSTAT=IOS, ERR=900 ) 
    READ ( LUNLUT, IOSTAT=IOS, ERR=900 ) RECORD
    CALL WRTLOG ( RECORD )
    DO WHILE ( RECORD(1:1) .EQ. '!' )
      READ ( LUNLUT, IOSTAT=IOS, ERR=900 ) RECORD
    END DO
  ELSE
    CALL OPNFIL ( LUNLUT, NAMLUT, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    READ ( LUNLUT, '(A)', IOSTAT=IOS, ERR=900 ) RECORD
  END IF     
!
! RECORD contains format identifier variable - checked in TABINF
!
! Next record contains Gas ID and array information
  IF ( BINFIL ) THEN
    READ ( LUNLUT, IOSTAT=IOS, ERR=900 ) RECORD
  ELSE
    READ ( LUNLUT, '(A)', IOSTAT=IOS, ERR=900 ) RECORD
  END IF 
!      
  GASSTR = RECORD(1:5) 
  IF ( INDEX ( GASSTR, '.' ) .EQ. 0 ) THEN
    READ ( GASSTR, * ) IDXMOL
    IDXISO = -1
  ELSE
    READ ( GASSTR, * ) RGAS
    IDXMOL = INT ( RGAS ) 
    IDXISO = NINT ( 10 * ( RGAS - IDXMOL ) )
  END IF            
!
  IGAS = IDXGAS ( IDXMOL, IDXISO ) 
  LUT(ILUT)%IGS = IGAS
!
  READ ( RECORD(6:), * ) NVTAB, V1TAB, V2TAB, DVTAB, NATAB, NPTAB, NTTAB, NQTAB
!
  OFFTAB = NTTAB .LT. 0
  NTTAB = ABS ( NTTAB )
!
  LUT(ILUT)%GHZ = NVTAB .LT. 0
  NVTAB = ABS ( NVTAB ) 
  IF ( LUT(ILUT)%GHZ ) THEN
    V1TAB = V1TAB * GHZ2CM
    V2TAB = V2TAB * GHZ2CM
  END IF
!
  LUT(ILUT)%NX = NATAB
  ALLOCATE ( LUT(ILUT)%TAB(NATAB,2) )
! 
  LUT(ILUT)%WNL = V1TAB
  LUT(ILUT)%WNU = V2TAB
! Set -ve wnos for buffers so REALUT reads in at least 2 records on 1st call
  LUT(ILUT)%WNO = (/ -2.0D0, -1.0D0 /)
  LUT(ILUT)%IBU = 2
!
  ALLOCATE ( PAXTAB(NPTAB), TEMTAB(NPTAB), VMRTAB(NPTAB,1) )
  ALLOCATE ( TAXTAB(NTTAB) )
  ALLOCATE ( QAXTAB(NQTAB) ) 
!
  IF ( BINFIL ) THEN
    READ ( LUNLUT, IOSTAT=IOS, ERR=900 ) PAXTAB
    READ ( LUNLUT, IOSTAT=IOS, ERR=900 ) TEMTAB
    READ ( LUNLUT, IOSTAT=IOS, ERR=900 ) VMRTAB
    READ ( LUNLUT, IOSTAT=IOS, ERR=900 ) TAXTAB
    READ ( LUNLUT, IOSTAT=IOS, ERR=900 ) QAXTAB
  ELSE
    READ ( LUNLUT, *, IOSTAT=IOS, ERR=900 ) PAXTAB
    READ ( LUNLUT, *, IOSTAT=IOS, ERR=900 ) TEMTAB
    READ ( LUNLUT, *, IOSTAT=IOS, ERR=900 ) VMRTAB
    READ ( LUNLUT, *, IOSTAT=IOS, ERR=900 ) TAXTAB
    READ ( LUNLUT, *, IOSTAT=IOS, ERR=900 ) QAXTAB
  END IF
!
  IF ( NDP .NE. 1 .OR. NDT .NE. 1 ) CALL SUBAXS ( NDP, NDT )
!
  NLKP = 1
  IF ( NPTAB .GT. 1 ) NLKP = NLKP * 2  ! interpolate in pressure
  IF ( NTTAB .GT. 1 ) NLKP = NLKP * 2  ! interpolate in tempreature
  IF ( NQTAB .GT. 1 ) NLKP = NLKP * 2  ! interpolate in VMR
  CALL LUTPTH ( ILUT, IGAS, NLKP ) 
!
  DEALLOCATE ( PAXTAB, TEMTAB, VMRTAB, TAXTAB, QAXTAB ) 
  IF ( ALLOCATED ( IXORG ) ) DEALLOCATE ( IXORG ) 
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL )  WRITE ( ERRMSG, * ) &
    'F-LUTTAB: I/O failure on LUT file. IOSTAT=', IOS
!
END SUBROUTINE LUTTAB
END MODULE LUTTAB_SUB
