MODULE SVDINF_SUB
CONTAINS
SUBROUTINE SVDINF ( NAMSVD, TAB, IDXMOL, IDXISO, V1, V2, FAIL, ERRMSG )
!
! VERSION
!   20JUL24 AD Checked.
!   26JAN18 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Extract info from SVD-LUT header 
!   Called by SVDFIL for each SVD-LUT file
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE SVDHDR_SUB ! Read SVD-LUT header record
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMSVD ! Name of SVD-LUT file
    CHARACTER(3),  INTENT(OUT) :: TAB    ! Tabulation function
    INTEGER(I4),   INTENT(OUT) :: IDXMOL ! RFM/HITRAN ID of gas
    INTEGER(I4),   INTENT(OUT) :: IDXISO ! Isotopic ID of gas, or -1
    REAL(R8),      INTENT(OUT) :: V1     ! Lower Wno limit [cm-1] of data
    REAL(R8),      INTENT(OUT) :: V2     ! Upper Wno limit [cm-1] of data
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IDUM   ! Dummy integer (NL, no singular vectors)
    INTEGER(I4)   :: IOS    ! Saved value of IOSTAT for error messages
    INTEGER(I4)   :: NV     ! No.spectral points in SVD-LUT file
    REAL(R8)      :: DV     ! Wavenumber [cm-1] interval for SVD-LUT file
    CHARACTER(80) :: RECORD ! Record read from SVD-LUT file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  OPEN ( UNIT=LUNTMP, FILE=NAMSVD, STATUS='OLD', ACTION='READ', &
         IOSTAT=IOS, ERR=900 )
  READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) RECORD
  READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) RECORD
  CALL WRTLOG ( RECORD )
  IF ( RECORD(1:1) .NE. '#' ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-SVDINF: 2nd record does not start with ''#'''
    RETURN
  END IF
!
  DO WHILE ( RECORD(1:1) .EQ. '#' )
    READ ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) RECORD
  END DO
!
  CALL SVDHDR ( RECORD, TAB, IDXMOL, IDXISO, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN
!
  READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) IDUM, NV, V1, DV
!
  IF ( NV .LT. 1 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-SVDINF: Irreg spectral grid LUT-SVD files not yet implemented'
    RETURN
  END IF
!
  V2 = V1 + ( NV - 1 ) * DV
  CLOSE ( LUNTMP, IOSTAT=IOS, ERR=900 )
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-SVDINF: Error reading LUT file header. IOSTAT=', IOS
!
END SUBROUTINE SVDINF
END MODULE SVDINF_SUB
