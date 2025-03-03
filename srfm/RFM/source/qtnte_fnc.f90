MODULE QTNTE_FNC
CONTAINS
REAL(R4) FUNCTION QTNTE ( IDXMOL, IDXISO, TEM, QVNTE )
!
! VERSION
!   10NOV23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Calculate total internal partition sums for non-LTE
!   Called by NTECLC.
!   The ratio of total internal partition sums, Q(296K)/Q(Path Temp) is 
!   calculated for a given molecule and isotpic species. The line strength 
!   can then be adjusted by multiplying by this factor.
!   For non-LTE an externally supplied Vib.Part.sum is used with a locally
!   generated rotation component:
!
!        QTNTE = QR296 * QV296 / ( QRTEM * QVNTE )
!
!   Based on subroutines supplied by R.GAMACHE (email).
!
! REFERENCES
!   Fischer, J. and R. R. Gamache
!     Partition sums for non-local thermodynamic equilibrium applications
!     J.Quant.Spect.Rad.Transfer 74, 273-284, 2002.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE PHYCON_DAT, ONLY: TEMREF ! Ref temperature (296K)
!
! SUBROUTINES
    USE QTFCT_FNC  ! Calculate total internal partition sum
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN)  :: IDXMOL ! HITRAN gas ID
    INTEGER(I4), INTENT(IN)  :: IDXISO ! HITRAN isotope ID 
    REAL(R4),    INTENT(IN)  :: TEM    ! Path temperature [K]
    REAL(R4),    INTENT(IN)  :: QVNTE  ! Non-LTE Vibrational Partition Fn.
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: MAXSPE = 40 ! Tot. No.isotope species (all gases)
!
! LOCAL VARIABLES
    INTEGER(I4) :: ISPE             ! Counter for isotope species (all gases)
    INTEGER(I4) :: J                ! Counter for DATA statements
    REAL(R8)    :: QRSTD            ! Rot.Part.Fn. at Ref.Temp 296 K
    REAL(R8)    :: QRTEM            ! Rot.Part.Fn at temperature TEM         
    REAL(R8)    :: QVSTD            ! Vib.Part.Fn. at Ref.Temp 296 K
    REAL(R8)    :: QVTEM            ! Local Double Precision QVNTE 
    REAL(R8)    :: QRCOEF(MAXSPE,4) ! Rotational Partitional functions
    REAL(R8)    :: QVCOEF(MAXSPE,4) ! Vibrational Partitional functions
    REAL(R8)    :: TEMPTH           ! Local version of TEM
!
! DATA STATEMENTS
!
! Rotational components
!...H2O     --   161
    DATA (QRCOEF( 1,J),J=1,4)/-.49589E+01,0.28147E+00, 0.12848E-02,-.58343E-06/
!...H2O     --   181
    DATA (QRCOEF( 2,J),J=1,4)/-.49589E+01,0.28147E+00, 0.12848E-02,-.58343E-06/
!...H2O     --   171
    DATA (QRCOEF( 3,J),J=1,4)/-.29625E+02,0.16822E+01, 0.76781E-02,-.34865E-05/
!...H2O     --   162
    DATA (QRCOEF( 4,J),J=1,4)/-.25676E+02,0.13745E+01, 0.63618E-02,-.28991E-05/
!...H2O     --   182
    DATA (QRCOEF( 5,J),J=1,4)/-.26018E+02,0.13915E+01, 0.64422E-02,-.29359E-05/
!...H2O     --   172 
    DATA (QRCOEF( 6,J),J=1,4)/-.15514E+03,0.83009E+01, 0.38425E-01,-.17511E-04/
!
! M.Lopez-Puertas & C.Mertens (1998) have suggested that assuming a simple
! rigid rotator gives better results than using these coefficients: 
! Qr(T)/Qr(T0) = (T/T0)**m where m=1 for a linear molecule such as CO2 
! (m=3/2 for other molecules with 3 rotational degrees of freedom)
! This can be implemented by setting each QRCOEF to / 0.,1.,0.,0./
!...CO2     --   626
    DATA (QRCOEF( 7,J),J=1,4)/0.16708E+00,0.89066E+00, 0.10000E-04,0.10901E-11/
!...CO2     --   636 
    DATA (QRCOEF( 8,J),J=1,4)/0.33418E+00,0.17812E+01, 0.19997E-04,0.11721E-11/
!...CO2     --   628 
    DATA (QRCOEF( 9,J),J=1,4)/0.33411E+00,0.18879E+01, 0.21192E-04,0.71375E-11/
!...CO2     --   627 
    DATA (QRCOEF(10,J),J=1,4)/0.20048E+01,0.11015E+02, 0.12444E-03,0.35048E-10/
!...CO2     --   638 
    DATA (QRCOEF(11,J),J=1,4)/0.66820E+00,0.37759E+01, 0.42377E-04,0.12864E-10/
!...CO2     --   637  
    DATA (QRCOEF(12,J),J=1,4)/0.40095E+01,0.22031E+02, 0.24705E-03,0.68948E-10/
!...CO2     --   828  
    DATA (QRCOEF(13,J),J=1,4)/0.16702E+00,0.10021E+01, 0.11248E-04,0.50684E-11/
!...CO2     --   728 
    DATA (QRCOEF(14,J),J=1,4)/0.20043E+01,0.11684E+02, 0.13096E-03,0.51947E-10/
! Gamache coefficients for unused isotope 
!...CO2     --   727 
!      DATA (QRCOEF( 0,J),J=1,4)/0.60138E+01,0.34074E+02,
!     &               0.38252E-03,0.11807E-09/
!
!...O3      --   666 
    DATA (QRCOEF(15,J),J=1,4)/-.10882E+03,0.52320E+01, 0.24982E-01,-.11471E-04/
!...O3      --   668 
    DATA (QRCOEF(16,J),J=1,4)/-.23242E+03,0.11171E+02, 0.53349E-01,-.24498E-04/
!...O3      --   686 
    DATA (QRCOEF(17,J),J=1,4)/-.11357E+03,0.54601E+01, 0.26072E-01,-.11972E-04/
!...O3      --   667 
    DATA (QRCOEF(18,J),J=1,4)/-.13512E+04,0.64955E+02, 0.31017E+00,-.14243E-03/
!...O3      --   676 
    DATA (QRCOEF(19,J),J=1,4)/-.66755E+03,0.32096E+02, 0.15326E+00,-.70373E-04/
!...O3      --   886 
!      DATA (QRCOEF( 0,J),J=1,4)/-.24276E+03,0.11668E+02,
!     &               0.55719E-01,-.25586E-04/
!...O3      --   868 
!      DATA (QRCOEF( 0,J),J=1,4)/-.12427E+03,0.59713E+01,
!     &               0.28519E-01,-.13096E-04/
!...O3      --   678 
!      DATA (QRCOEF( 0,J),J=1,4)/-.14264E+04,0.68559E+02,
!     &               0.32741E+00,-.15034E-03/
!...O3      --   687 
!      DATA (QRCOEF( 0,J),J=1,4)/-.14107E+04,0.67815E+02,
!     &               0.32384E+00,-.14870E-03/
!...O3      --   768 
!      DATA (QRCOEF( 0,J),J=1,4)/-.14440E+04,0.69394E+02,
!     &               0.33141E+00,-.15218E-03/
!...O3      --   776 
!      DATA (QRCOEF( 0,J),J=1,4)/-.82908E+04,0.39855E+03,
!     &               0.19032E+01,-.87394E-03/
!...O3      --   767 
!      DATA (QRCOEF( 0,J),J=1,4)/-.41960E+04,0.20168E+03,
!     &               0.96312E+00,-.44226E-03/
!...O3      --   888 
!      DATA (QRCOEF( 0,J),J=1,4)/-.12988E+03,0.62405E+01,
!     &               0.29805E-01,-.13687E-04/
!...O3      --   887 
!      DATA (QRCOEF( 0,J),J=1,4)/-.15077E+04,0.72452E+02,
!     &               0.34602E+00,-.15889E-03/
!...O3      --   878 
!      DATA (QRCOEF( 0,J),J=1,4)/-.76294E+03,0.36658E+02,
!     &               0.17508E+00,-.80399E-04/
!...O3      --   778 
!      DATA (QRCOEF( 0,J),J=1,4)/-.88635E+04,0.42594E+03,
!     &               0.20342E+01,-.93412E-03/
!...O3      --   787 
!      DATA (QRCOEF( 0,J),J=1,4)/-.43824E+04,0.21063E+03,
!     &               0.10059E+01,-.46190E-03/
!...O3      --   777 
!      DATA (QRCOEF( 0,J),J=1,4)/-.25751E+05,0.12377E+04,
!     &               0.59105E+01,-.27141E-02/
!
!...N2O     --    446 
    DATA (QRCOEF(20,J),J=1,4)/0.12952E+02,0.15381E+02, -.10496E-01,0.52270E-04/
!...N2O     --    456 
    DATA (QRCOEF(21,J),J=1,4)/0.20050E+01,0.99541E+01, 0.11354E-03,0.48260E-10/
!...N2O     --    546 
    DATA (QRCOEF(22,J),J=1,4)/0.20048E+01,0.10301E+02, 0.11746E-03,0.50698E-10/
!...N2O     --    448 
    DATA (QRCOEF(23,J),J=1,4)/0.30070E+01,0.15815E+02, 0.18051E-03,0.12112E-09/
!...N2O     --    447 
    DATA (QRCOEF(24,J),J=1,4)/0.18043E+02,0.92299E+02, 0.10529E-02,0.51969E-09/
!
!...CO      --    26 
    DATA (QRCOEF(25,J),J=1,4)/0.33737E+00,0.36156E+00, 0.91570E-06,-.75455E-10/
!...CO      --    36 
    DATA (QRCOEF(26,J),J=1,4)/0.67434E+00,0.75639E+00, 0.18989E-05,-.14272E-09/
!...CO      --    28
    DATA (QRCOEF(27,J),J=1,4)/0.33716E+00,0.37964E+00, 0.95167E-06,-.68683E-10/
!...CO      --    27 
    DATA (QRCOEF(28,J),J=1,4)/0.20234E+01,0.22256E+01, 0.56019E-05,-.43383E-09/
!...CO      --    38  
    DATA (QRCOEF(29,J),J=1,4)/0.67393E+00,0.79605E+00, 0.19825E-05,-.13639E-09/
!...CO      --    37  
    DATA (QRCOEF(30,J),J=1,4)/0.40454E+01,0.46408E+01, 0.36191E-04,-.10571E-08/
!
!...CH4     --    211 
    DATA (QRCOEF(31,J),J=1,4)/-.18312E+02,0.93325E+00, 0.43490E-02,-.20147E-05/
!...CH4     --    311 
    DATA (QRCOEF(32,J),J=1,4)/-.36622E+02,0.18664E+01, 0.86974E-02,-.40291E-05/
!...CH4     --    212 
    DATA (QRCOEF(33,J),J=1,4)/-.14977E+03,0.75332E+01, 0.35251E-01,-.16341E-04/
! 
!...NO      --    46 
    DATA (QRCOEF(34,J),J=1,4)/-.35215E+02,0.27383E+01, 0.56720E-02,-.53671E-05/
!...NO      --    56 
    DATA (QRCOEF(35,J),J=1,4)/-.21715E+02,0.18926E+01, 0.39182E-02,-.37069E-05/
!...NO      --    48 
    DATA (QRCOEF(36,J),J=1,4)/-.37147E+02,0.28869E+01, 0.59844E-02,-.56624E-05/
! 
!...NO2     --   646 
    DATA (QRCOEF(37,J),J=1,4)/-.43332E+03,0.20816E+02, 0.98876E-01,-.46097E-04/
!
!...OH      --    61 
    DATA (QRCOEF(38,J),J=1,4)/0.77363E+01,0.17152E+00, 0.34886E-03,-.33504E-06/
!...OH      --    81 
    DATA (QRCOEF(39,J),J=1,4)/0.76867E+01,0.17321E+00, 0.34997E-03,-.33611E-06/
!...OH      --    62 
    DATA (QRCOEF(40,J),J=1,4)/0.83072E+01,0.45959E+00, 0.99174E-03,-.92504E-06/
! 
! Vibrational components
!...H2O     --   161 
    DATA (QVCOEF( 1,J),J=1,4)/0.99842E+00,0.29713E-04, -.17345E-06,0.32366E-09/
!...H2O     --   181 
    DATA (QVCOEF( 2,J),J=1,4)/0.99841E+00,0.29843E-04, -.17490E-06,0.32764E-09/
!...H2O     --   171 
    DATA (QVCOEF( 3,J),J=1,4)/0.99841E+00,0.29836E-04, -.17447E-06,0.32606E-09/
!...H2O     --   162 
    DATA (QVCOEF( 4,J),J=1,4)/0.99821E+00,0.36179E-04, -.23230E-06,0.48406E-09/
!...H2O     --   182 
    DATA (QVCOEF( 5,J),J=1,4)/0.99821E+00,0.36264E-04, -.23407E-06,0.49050E-09/
!...H2O     --   172 
    DATA (QVCOEF( 6,J),J=1,4)/0.99821E+00,0.36224E-04, -.23319E-06,0.48735E-09/
! 
!...CO2     --   626 
    DATA (QVCOEF( 7,J),J=1,4)/0.10384E+01,-.64414E-03, 0.26438E-05,0.21323E-09/
!...CO2     --   636 
    DATA (QVCOEF( 8,J),J=1,4)/0.10408E+01,-.69739E-03, 0.29539E-05,0.10076E-10/
!...CO2     --   628 
    DATA (QVCOEF( 9,J),J=1,4)/0.10389E+01,-.65708E-03, 0.27190E-05,0.17163E-09/
!...CO2     --   627 
    DATA (QVCOEF(10,J),J=1,4)/0.10385E+01,-.64989E-03, 0.26847E-05,0.18280E-09/
!...CO2     --   638 
    DATA (QVCOEF(11,J),J=1,4)/0.10414E+01,-.70990E-03, 0.30216E-05,0.34863E-11/
!...CO2     --   637 
    DATA (QVCOEF(12,J),J=1,4)/0.10410E+01,-.70345E-03, 0.29864E-05,0.11948E-10/
!...CO2     --   828  
    DATA (QVCOEF(13,J),J=1,4)/0.10397E+01,-.67365E-03, 0.28109E-05,0.11336E-09/
!...CO2     --   728 
    DATA (QVCOEF(14,J),J=1,4)/0.10393E+01,-.66478E-03, 0.27615E-05,0.14840E-09/
!...CO2     --   727 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10390E+01,-.65785E-03,
!     &               0.27232E-05,0.16976E-09/
!
!...O3      --   666  
    DATA (QVCOEF(15,J),J=1,4)/0.10198E+01,-.28986E-03, 0.83514E-06,0.15130E-08/
!...O3      --   668 
    DATA (QVCOEF(16,J),J=1,4)/0.10213E+01,-.31931E-03, 0.98160E-06,0.14407E-08/
!...O3      --   686 
    DATA (QVCOEF(17,J),J=1,4)/0.10215E+01,-.31701E-03, 0.93195E-06,0.15468E-08/
!...O3      --   667 
    DATA (QVCOEF(18,J),J=1,4)/0.10208E+01,-.30728E-03, 0.91502E-06,0.14740E-08/
!...O3      --   676 
    DATA (QVCOEF(19,J),J=1,4)/0.10209E+01,-.30652E-03, 0.89097E-06,0.15275E-08/
!...O3      --   886 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10233E+01,-.35115E-03,
!     &               0.11010E-05,0.14519E-08/
!...O3      --   868 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10231E+01,-.35270E-03,
!     &               0.11405E-05,0.13660E-08/
!...O3      --   678 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10224E+01,-.33660E-03,
!     &               0.10462E-05,0.14394E-08/
!...O3      --   687 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10226E+01,-.33677E-03,
!     &               0.10246E-05,0.14950E-08/
!...O3      --   768 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10227E+01,-.34052E-03,
!     &               0.10617E-05,0.14198E-08/
!...O3      --   776 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10215E+01,-.31895E-03,
!     &               0.95379E-06,0.15073E-08/
!...O3      --   767 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10213E+01,-.32016E-03,
!     &               0.98582E-06,0.14387E-08/
!...O3      --   888 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10247E+01,-.38169E-03,
!     &               0.12752E-05,0.13361E-08/
!...O3      --   887 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10240E+01,-.36730E-03,
!     &               0.11967E-05,0.13866E-08/
!...O3      --   878 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10240E+01,-.36855E-03,
!     &               0.12140E-05,0.13467E-08/
!...O3      --   778 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10234E+01,-.35416E-03,
!     &               0.11240E-05,0.14176E-08/
!...O3      --   787 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10235E+01,-.35268E-03,
!     &               0.11053E-05,0.14557E-08/
!...O3      --   777 
!      DATA (QVCOEF( 0,J),J=1,4)/0.10226E+01,-.33759E-03,
!     &               0.10450E-05,0.14567E-08/
!
!...N2O     --    446 
    DATA (QVCOEF(20,J),J=1,4)/0.10493E+01,-.88790E-03, 0.40828E-05,-.66527E-09/
!...N2O     --    456 
    DATA (QVCOEF(21,J),J=1,4)/0.10513E+01,-.93420E-03, 0.43832E-05,-.86356E-09/
!...N2O     --    546 
    DATA (QVCOEF(22,J),J=1,4)/0.10494E+01,-.89574E-03, 0.41556E-05,-.69976E-09/
!...N2O     --    448 
    DATA (QVCOEF(23,J),J=1,4)/0.10487E+01,-.88676E-03, 0.41275E-05,-.64301E-09/
!...N2O     --    447 
    DATA (QVCOEF(24,J),J=1,4)/0.10496E+01,-.89456E-03, 0.41231E-05,-.66366E-09/
!
!...CO      --    26 
    DATA (QVCOEF(25,J),J=1,4)/0.99944E+00,0.96257E-05, -.50262E-07,0.81728E-10/
!...CO      --    36 
    DATA (QVCOEF(26,J),J=1,4)/0.99936E+00,0.11045E-04, -.57272E-07,0.93156E-10/
!...CO      --    28 
    DATA (QVCOEF(27,J),J=1,4)/0.99939E+00,0.10836E-04, -.57613E-07,0.94711E-10/
!...CO      --    27 
    DATA (QVCOEF(28,J),J=1,4)/0.99940E+00,0.10383E-04, -.53727E-07,0.87169E-10/
!...CO      --    38 
    DATA (QVCOEF(29,J),J=1,4)/0.99930E+00,0.12213E-04, -.64899E-07,0.10734E-09/
!...CO      --    37 
    DATA (QVCOEF(30,J),J=1,4)/0.99944E+00,0.96257E-05, -.50262E-07,0.81728E-10/
!
!...CH4     --    211 
    DATA (QVCOEF(31,J),J=1,4)/0.99109E+00,0.18230E-03, -.11900E-05,0.25319E-08/
!...CH4     --    311 
    DATA (QVCOEF(32,J),J=1,4)/0.99109E+00,0.18230E-03, -.11900E-05,0.25319E-08/
!...CH4     --    212 
    DATA (QVCOEF(33,J),J=1,4)/0.99109E+00,0.18230E-03, -.11900E-05,0.25319E-08/
!
!...NO      --    46 
    DATA (QVCOEF(34,J),J=1,4)/0.99899E+00,0.17970E-04, -.97249E-07,0.16549E-09/
!...NO      --    56 
    DATA (QVCOEF(35,J),J=1,4)/0.99893E+00,0.19170E-04, -.10469E-06,0.17976E-09/
!...NO      --    48 
    DATA (QVCOEF(36,J),J=1,4)/0.99891E+00,0.19684E-04, -.10817E-06,0.18708E-09/
!
!...NO2     --   646 
    DATA (QVCOEF(37,J),J=1,4)/0.10129E+01,-.19035E-03, 0.55618E-06,0.91873E-09/
!
!...OH      --    61 
    DATA (QVCOEF(38,J),J=1,4)/0.99999E+00,0.20075E-06, -.93665E-09,0.13295E-11/
!...OH      --    81 
    DATA (QVCOEF(39,J),J=1,4)/0.99984E+00,0.26827E-05, -.13312E-07,0.20502E-10/
!...OH      --    62 
    DATA (QVCOEF(40,J),J=1,4)/0.99999E+00,0.20075E-06, -.93665E-09,0.13295E-11/
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  SELECT CASE ( IDXMOL ) 
  CASE ( 1 ) ; ISPE = IDXISO       ! H2O
  CASE ( 2 ) 
    ISPE = 6 + IDXISO              ! CO2
    IF ( IDXISO .EQ. 9 ) ISPE = 13 ! No data for 838, so duplicate 828
  CASE ( 3 ) ; ISPE = 14 + IDXISO  ! O3
  CASE ( 4 ) ; ISPE = 19 + IDXISO  ! N2O
  CASE ( 5 ) ; ISPE = 24 + IDXISO  ! CO
  CASE ( 6 ) ; ISPE = 30 + IDXISO  ! CH4
  CASE ( 7 ) ; ISPE = 33 + IDXISO  ! NO
  CASE ( 8 ) ; ISPE = 36 + IDXISO  ! NO2
  CASE ( 9 ) ; ISPE = 37 + IDXISO  ! OH
  CASE DEFAULT
    QTNTE = QTFCT ( IDXMOL, IDXISO, TEM )
  END SELECT
!
  TEMPTH = DBLE ( TEM )
  QVTEM = DBLE ( QVNTE )
!
  QRTEM = QRCOEF(ISPE,1) + TEMPTH * ( QRCOEF(ISPE,2) + &
          TEMPTH * ( QRCOEF(ISPE,3) + TEMPTH * QRCOEF(ISPE,4) ))
  QRSTD = QRCOEF(ISPE,1) + &
          DBLE ( TEMREF ) * ( QRCOEF(ISPE,2) + &
          DBLE ( TEMREF ) * ( QRCOEF(ISPE,3) + &
          DBLE ( TEMREF ) *   QRCOEF(ISPE,4)   ))
  QVSTD = QVCOEF(ISPE,1) + &
          DBLE ( TEMREF ) * ( QVCOEF(ISPE,2) + &
          DBLE ( TEMREF ) * ( QVCOEF(ISPE,3) + &
          DBLE ( TEMREF ) *   QVCOEF(ISPE,4)   ))
!
  QTNTE = SNGL ( QRSTD * QVSTD / QRTEM / QVTEM )
!
END FUNCTION QTNTE
END MODULE QTNTE_FNC
