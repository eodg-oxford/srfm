MODULE REASVD_SUB
CONTAINS
SUBROUTINE REASVD ( LUNSVD, NAMSVD, ISVD, NLREQ, FAIL, ERRMSG )
!
! VERSION
!   14JUN24 AD Checked.
!   19DEC17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Load SVD-compressed LUT data from file
!   Called by INISVD for each SVD-LUT file.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE SVDCOM_DAT ! SVD-compressed LUT data
!
! SUBROUTINES
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
    USE SVDGRD_SUB ! Thin out SVD LUT matrix to match (irregular) fine grid
    USE SVDHDR_SUB ! Read SVD-LUT header record
    USE SVDLIM_SUB ! Limit No. Singular Vectors used in SVD-LUT
    USE SVDPTH_SUB ! Set SVD path interpolation
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNSVD ! LUN for accessing SVD-LUT data
    CHARACTER(*),  INTENT(IN)  :: NAMSVD ! Name of SVD-LUT file
    INTEGER(I4),   INTENT(IN)  :: ISVD   ! Index in /SVDCOM/ for storing data
    INTEGER(I4),   INTENT(IN)  :: NLREQ  ! No SVs reqd (0=use all)
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IDXISO ! Isotopic ID of gas, read from SVD-LUT file
    INTEGER(I4)   :: IDXMOL ! HITRAN ID of gas, read from SVD-LUT file 
    INTEGER(I4)   :: IOS    ! Saved value of IOSTAT for error messages
    INTEGER(I4)   :: NG     ! No. spectral points after matching to RFM grid
    INTEGER(I4)   :: NL     ! No. of S.Vs.
    INTEGER(I4)   :: NP     ! No. -lnp axis points
    INTEGER(I4)   :: NT     ! No. Tem axis points
    INTEGER(I4)   :: NV     ! No. of spectral points
    INTEGER(I4)   :: NX     ! No. p x No. T points
    REAL(R4)      :: DP     ! -lnp axis increment
    REAL(R4)      :: DT     ! Tem axis increment
    REAL(R4)      :: P1     ! 1st -lnp axis point
    REAL(R4)      :: T1     ! 1st Tem axis point
    REAL(R8)      :: DV     ! Wavenumber axis increment [cm-1]
    REAL(R8)      :: V1     ! Lower wavenumber axis value [cm-1]
    CHARACTER(80) :: RECORD ! Record read from SVD-LUT file
    REAL(R4), ALLOCATABLE :: UT(:,:) ! Transpose of U-matrix
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  OPEN ( UNIT=LUNSVD, FILE=NAMSVD, STATUS='OLD', ACTION='READ', IOSTAT=IOS, &
         ERR=900 ) 
  READ ( LUNSVD, '(A)', IOSTAT=IOS, ERR=900 ) RECORD
  READ ( LUNSVD, '(A)', IOSTAT=IOS, ERR=900 ) RECORD
  IF ( RECORD(1:1) .EQ. '#' ) THEN
    CALL WRTLOG ( RECORD ) 
  ELSE 
    FAIL = .TRUE. 
    ERRMSG = 'F-REASVD: 2nd record of file does not start'// &
             ' with expected ''#'' character'
    RETURN
  END IF
!
  DO WHILE ( RECORD(1:1) .EQ. '#' ) 
    READ ( LUNSVD, '(A)', IOSTAT=IOS, ERR=900, END=900 ) RECORD
  END DO
!     
  CALL SVDHDR ( RECORD, SVD(ISVD)%TAB, IDXMOL, IDXISO, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN
!
  SVD(ISVD)%IGS = IDXGAS ( IDXMOL, IDXISO ) 
!
  READ ( LUNSVD, *, IOSTAT=IOS, ERR=900 ) NL, NV, V1, DV, NP, P1, DP, NT, T1, DT
  IF ( NP .EQ. 1 ) DP = 1.0
  IF ( NT .EQ. 1 ) DT = 1.0
!
  SVD(ISVD)%NL = NL
  ALLOCATE ( UT(NL,NV) ) 
  READ ( LUNSVD, *, IOSTAT=IOS, ERR=900 ) UT
! Sample UT in spectral domain to match full (irregular) grid points
  CALL SVDGRD ( UT, NV, V1, DV, SVD(ISVD)%IGO, NG ) 
!
  ALLOCATE ( SVD(ISVD)%U(NG,NL) )
  SVD(ISVD)%U = TRANSPOSE ( UT(:,1:NG) ) 
!
  NX = NP * NT
  SVD(ISVD)%NX = NX 
  ALLOCATE ( SVD(ISVD)%K(NL,NX) )
  READ ( LUNSVD, *, IOSTAT=IOS, ERR=900 ) SVD(ISVD)%K
!
  CLOSE ( LUNSVD, IOSTAT=IOS, ERR=900 )
!
  IF ( NLREQ .NE. 0 ) THEN  ! Limit no. SV.s if required
    CALL SVDLIM ( NLREQ, SVD(ISVD)%NL ) 
    IF ( SVD(ISVD)%NL .EQ. 0 ) SVD(ISVD)%TAB = 'lin'  ! set 'lin' for zero abs.
  END IF
!
  CALL SVDPTH ( ISVD, NP, P1, DP, NT, T1, DT ) 
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-REASVD: I/O failure on SVD LUT file. IOSTAT=', IOS
!
END SUBROUTINE REASVD
END MODULE REASVD_SUB
