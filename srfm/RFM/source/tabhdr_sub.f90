MODULE TABHDR_SUB
CONTAINS
SUBROUTINE TABHDR ( LUN, IGAS, FAIL, ERRMSG )
!
! VERSION
!   30OCT23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Write TAB file header records
!   Called once by SPCTAB for each output .tab file
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FULCOM_DAT ! Full grid data
    USE GASCOM_DAT ! Molecule and isotope data
    USE HDRCOM_DAT ! Output header data
    USE TABCOM_DAT ! Axes for tabulated absorption coefficients
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUN    ! LUN for .tab file
    INTEGER(I4),   INTENT(IN)  :: IGAS   ! Absorber index
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error occurs
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: RFMT = 1.0    ! .tab file format identifier
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IDM    ! HITRAN/RFM index of molecule
    INTEGER(I4)   :: IOS    ! Value of IOSTAT saved for error messages.
    INTEGER(I4)   :: ISIGNT ! Sign for NTTAB in output file
    INTEGER(I4)   :: ISIGNV ! Sign for NVTAB in output file
    CHARACTER(5)  :: GASSTR ! Molecule ID, incl. isotope
    CHARACTER(80) :: HEADR1 ! Header records for output files
    CHARACTER(80) :: HEADR2 !
    CHARACTER(80) :: HEADR3 ! 
    CHARACTER(80) :: HEADR4 ! Format identifier
    CHARACTER(80) :: HEADR5 ! Molec/Iso ID, spec range and axis sizes 
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( GHZTAB .EQ. 1.0D0 ) THEN
    ISIGNV = 1              
  ELSE
    ISIGNV = -1   ! indicate output spectral axis in GHz
  END IF
!
  IF ( OFFTAB ) THEN
    ISIGNT = -1 
  ELSE
    ISIGNT = 1
  END IF
!
  HEADR1 = '! ' // GAS(IGAS)%COD // ' Tabulated Absorp.Coeff. ' &
           // 'created by RFM v.' // VIDHDR // '                '
  HEADR2 = TXTHDR
  HEADR3 = '! GAS        NV      V1          V2        DV' // &
           '        NA      NP      NT      NQ'
  WRITE ( HEADR4, '(F10.2)' ) RFMT
!
! Find out if this gas has been split into isotopes
  IDM = GAS(IGAS)%IDM
  IF ( ISOMOL(IDM) ) THEN
    WRITE ( GASSTR, '(I3,A,I1)' ) IDM, '.', GAS(IGAS)%IDI
  ELSE
    WRITE ( GASSTR, '(I5)' ) IDM
  END IF
  WRITE ( HEADR5, '(A5,I10,2F12.4,F8.4,4I8)' ) &
    GASSTR, ISIGNV*NFUL, GHZTAB*WNLFUL,  GHZTAB*WNUFUL, GHZTAB*WNRFUL, &
    NATAB, NPTAB, ISIGNT*NTTAB, NQTAB
!
  IF ( BINTAB ) THEN
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR1
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR2
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR3
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR4
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR5
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) PAXTAB
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) TEMTAB
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) VMRTAB(:,IGAS) 
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) TAXTAB
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) QAXTAB
  ELSE
    WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR1
    WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR2
    WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR3
    WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR4
    WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR5
    WRITE ( LUN, *, IOSTAT=IOS, ERR=900 ) PAXTAB
    WRITE ( LUN, *, IOSTAT=IOS, ERR=900 ) TEMTAB
    WRITE ( LUN, *, IOSTAT=IOS, ERR=900 ) VMRTAB(:,IGAS)
    WRITE ( LUN, *, IOSTAT=IOS, ERR=900 ) TAXTAB
    WRITE ( LUN, *, IOSTAT=IOS, ERR=900 ) QAXTAB
  END IF
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-TABHDR: I/O failure on output file. IOSTAT=', IOS
!
END SUBROUTINE TABHDR
END MODULE TABHDR_SUB
