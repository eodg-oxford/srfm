MODULE WRTSTT_SUB
CONTAINS
SUBROUTINE WRTSTT ( ISPC, FAIL, ERRMSG )
!
! VERSION
!   02APR23 AD Checked.
!   08NOV17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Write widemesh statistics
!   Called by SPCOUT for each spectral interval if WID flag enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE HDRCOM_DAT ! Output header data
    USE IDXCON_DAT ! RFM/HITRAN indices for specific molecules
    USE STTCOM_DAT ! Widemesh statistics
    USE NAMCOM_DAT, ONLY: WIDNAM ! Name of widemesh stats output file
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE OPNOUT_SUB ! Open spectral output files for current spectral range
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: ISPC   ! Spectral range#
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error occurs
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: ICT  ! 1=used ctm data for this interval, 0=no data
    INTEGER(I4) :: IGAS ! Counter for absorbers
    INTEGER(I4) :: IOS  ! Saved value of IOSTAT for error message
    INTEGER(I4) :: IWC1 ! Lower WM index of continuum data
    INTEGER(I4) :: IWC2 ! Upper WM index of continuum data
    INTEGER(I4) :: IWID ! Counter for WM intervals
    INTEGER(I4) :: IXS  ! 1=used x/s data for this interval, 0=no data
    INTEGER(I4) :: MGAS ! Starting absorber index for output files
    INTEGER(I4) :: NEX  ! No lines outside spc.range
    INTEGER(I4) :: NIN  ! No lines inside spc.range
    INTEGER(I4) :: NLF  ! No. lines contributing to finemesh interval
    INTEGER(I4) :: NLI  ! No. lines in widemesh interval
    INTEGER(I4) :: NLW  ! No. lines contributing to widemesh interval
    INTEGER(I4) :: NQD  ! No. paths using quadratic interpolation
    REAL(R8)    :: WNCL ! Lower limit of tabulated continuuum data
    REAL(R8)    :: WNCU ! Upper limit of tabulated continuuum data
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! One file per absorbing species, plus one combined, unless only one absorber
  IF ( NSTG .EQ. 1 ) THEN
    MGAS = 1
  ELSE
    MGAS = 0
  END IF
!
! Continuum data flags rather crudely handled given known ranges of continua
! for various molecules
  DO IGAS = 1, NGAS
    IF ( GAS(IGAS)%CTM ) THEN
      SELECT CASE ( GAS(IGAS)%IDM ) 
      CASE ( IDXH2O ) 
        WNCL = 0.0D0
        WNCU = 20000.0D0
      CASE ( IDXCO2 )
        WNCL = 0.0D0
        WNCU = 4000.0D0
      CASE ( IDXO2 ) 
        WNCL = 1360.0D0
        WNCU = 1805.0D0
      CASE ( IDXN2 ) 
        WNCL = 2120.0d0
        WNCU = 2605.0D0
      CASE DEFAULT
        STOP 'F-WRTSTT: Logical error'  ! Unexpected continuum molecule
      END SELECT
      IWC1 = INT ( ( WNCL - WN1STT ) / WNDSTT ) + 1 
      IWC1 = MAX ( 1, IWC1 ) 
      IWC2 = INT ( ( WNCU - WN1STT ) / WNDSTT ) + 1 
      IWC2 = MIN ( IWC2, NSTW ) 
      STT(IWC1:IWC2,IGAS)%ICT = 1
    END IF
  END DO
!
  DO IGAS = MGAS, NGAS
    CALL OPNOUT ( LUNTMP, WIDNAM, FAIL, ERRMSG, ISPC=ISPC, IGAS=IGAS )
    IF ( IGAS .EQ. 0 ) THEN
      WRITE ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) &
        '! Total Widemesh Line count created by RFM v.'// VIDHDR
    ELSE
      WRITE ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) &
        '! ' // GAS(IGAS)%COD // ' Widemesh Line count created by RFM v.' &
        // VIDHDR
    END IF
    WRITE ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) TXTHDR
    WRITE ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) NSTW, WN1STT, WNDSTT, &
      ' =NWid,Wno1,DWno'
    WRITE ( LUNTMP, '(A)', IOSTAT=IOS, ERR=900 ) 'Itvl#  Wavenumber  ' &
      // 'Tot.L =(WidL + FinL) Loc.L X/S CTM NQd'
!
    DO IWID = 1, NSTW
      IF ( IGAS .EQ. 0 ) THEN
        NLW = SUM ( STT(IWID,:)%NLW )
        NLF = SUM ( STT(IWID,:)%NLF )
        NLI = SUM ( STT(IWID,:)%NLI )
        IXS = SUM ( STT(IWID,:)%IXS )
        ICT = SUM ( STT(IWID,:)%ICT ) 
        NQD = SUM ( STT(IWID,:)%NQD ) 
        NIN = SUM ( NINSTT ) 
        NEX = SUM ( NEXSTT ) 
      ELSE
        NLW = STT(IWID,IGAS)%NLW
        NLF = STT(IWID,IGAS)%NLF
        NLI = STT(IWID,IGAS)%NLI
        IXS = STT(IWID,IGAS)%IXS
        ICT = STT(IWID,IGAS)%ICT
        NQD = STT(IWID,IGAS)%NQD
        NIN = NINSTT(IGAS)
        NEX = NEXSTT(IGAS)
      END IF
      WRITE ( LUNTMP, '(I5,F12.4,4I7,I3,I4,I5)', IOSTAT=IOS, ERR=900 ) &
        IWID, WN1STT + (IWID-1) * WNDSTT, (NLW+NLF), NLW, NLF, NLI, &
        IXS, ICT, NQD
    END DO
!
    WRITE ( LUNTMP, * ) &
      NIN, NEX, ' = Internal, external lines used for whole range'
!
    CLOSE ( LUNTMP, IOSTAT=IOS, ERR=900 )
  END DO
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-WRTSTT: I/O failure on output file. IOSTAT=', IOS
!
END SUBROUTINE WRTSTT
END MODULE WRTSTT_SUB

