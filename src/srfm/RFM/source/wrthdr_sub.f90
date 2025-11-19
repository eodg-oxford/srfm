MODULE WRTHDR_SUB
CONTAINS
SUBROUTINE WRTHDR ( LUN, TYP, IJAC, ILEV, ISPC, ITAN, FAIL, ERRMSG )
!
! VERSION
!   27MAY24 AD Checked.
!   25MAR19 AD Use USRUNI for length units if HOM flag. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Write header records of output files
!   Called by SPCWRT for each spectral output file
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE HDRCOM_DAT ! Output header data
    USE JACCOM_DAT ! Jacobian data
    USE LEVCOM_DAT ! Intermediate output levels
    USE SPCCOM_DAT ! Spectral range data
    USE TANCOM_DAT ! Tangent path data
    USE PHYCON_DAT, ONLY: GHZ2CM ! GHz to cm-1 conv.fac
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUN    ! Next available LUN
    CHARACTER(3),  INTENT(IN)  :: TYP    ! Type of spectrum 'ABS','COO',etc
    INTEGER(I4),   INTENT(IN)  :: IJAC   ! Index of Jacobian, or 0    
    INTEGER(I4),   INTENT(IN)  :: ILEV   ! Index of output level, or 0    
    INTEGER(I4),   INTENT(IN)  :: ISPC   ! Spectral range index 
    INTEGER(I4),   INTENT(IN)  :: ITAN   ! Index of output tan.path
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error occurs
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IGHFAC ! -1 if output in GHz, else 1
    INTEGER(I4)   :: IOS    ! Value of IOSTAT from WRITE statements
    INTEGER(I4)   :: IPT    ! Pointer to start of substring
    REAL(R8)      :: GHZFAC ! Wno to GHz conversion factor
    CHARACTER(80) :: HEADR1 ! 1st Header record for output files
    CHARACTER(80) :: HEADR2 ! Header from driver table *HDR section
    CHARACTER(80) :: HEADR3 ! titles, or Jacobian information (JACFLG)
    CHARACTER(15) :: LABEL  ! Label written to file header
    CHARACTER(32) :: PTHTXT ! Text describing viewing geometry
    CHARACTER(13) :: TYPTXT ! Text describing type of spectrum
    CHARACTER(10) :: WNOGHZ ! 'Wavenumber' or 'Frequency ' in HEADR3
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  SELECT CASE ( TYP )
  CASE ( 'ABS' ) ; TYPTXT = 'Absorption'
  CASE ( 'BBT' ) ; TYPTXT = 'Bright. Temp.'
  CASE ( 'COO' ) ; TYPTXT = 'Cooling Rate'
  CASE ( 'OPT' ) ; TYPTXT = 'Optical Depth'
  CASE ( 'RAD' ) ; TYPTXT = 'Radiance'
  CASE ( 'RJT' ) ; TYPTXT = 'Rayl-J.Temp.'
  CASE ( 'TRA' ) ; TYPTXT = 'Transmission'
  CASE DEFAULT
    STOP 'F-WRTHDR: logical error'
  END SELECT
!
! Set according to whether output is in wavenumber or GHz
  IF ( GHZFLG ) THEN
    WNOGHZ = 'Frequency '
  ELSE
    WNOGHZ = 'Wavenumber'
  END IF
!
! Set header records for output files ('123456789' replaced by tangent height)
  IF ( HOMFLG ) THEN
    PTHTXT = 'Homog. Path Length =' // TAN(ITAN)%STR // ' ' // USRUNI
  ELSE IF ( FLXFLG ) THEN
    PTHTXT = 'for Altitude Level =' // TAN(ITAN)%STR // ' km' 
  ELSE IF ( NADFLG ) THEN
    PTHTXT = 'Nadir Path Airmass =' // TAN(ITAN)%STR 
  ELSE IF ( ZENFLG ) THEN
    PTHTXT = 'Zenith Path Airmass=' // TAN(ITAN)%STR 
  ELSE IF ( USRELE ) THEN
    PTHTXT = 'View Elevation Ang.=' // TAN(ITAN)%STR // ' dg'
  ELSE IF ( USRGEO ) THEN
    PTHTXT = 'Limb Geom.Tang.Hgt =' // TAN(ITAN)%STR // ' km'
  ELSE                                    ! Normal limb-viewing geometry
    PTHTXT = 'Limb Path Tang.Hgt =' // TAN(ITAN)%STR // ' km'
  END IF
  IF ( FLXFLG ) THEN
    HEADR1 = '! Flux '//TYPTXT//' calc.'//PTHTXT//' by RFM v.'//VIDHDR
  ELSE
    HEADR1 = '! '//TYPTXT//' calc. for '//PTHTXT//' by RFM v.'//VIDHDR
  END IF
  IPT = INDEX ( HEADR1, '123456789' )
!
  HEADR2 = TXTHDR
! Note that HEADR3 is modified by JACHDR if JAC,MTX,LEV flags enabled 
  IF ( IJAC .EQ. 0 .AND. ILEV .EQ. 0 ) THEN
    HEADR3 = '!No.Pts  Lower_'//WNOGHZ//'   Delta_'//WNOGHZ//'  '// &
             'Upper_'//WNOGHZ//'  Label'
  ELSE IF ( ILEV .GT. 0 ) THEN            ! Intermediate output levels
    HEADR3 = '! Intermediate spectrum at level ' // &
             TRIM ( C9REAL ( LEV(ILEV)%HGT ) ) // ' km'
    IF ( LEV(ILEV)%IDR .EQ. 1 ) THEN
      HEADR3 = TRIM ( HEADR3 ) // ' upward path'
    ELSE IF ( LEV(ILEV)%IDR .EQ. -1 ) THEN
      HEADR3 = TRIM ( HEADR3 ) // ' downward path'
    END IF
  ELSE IF ( JAC(IJAC)%COL ) THEN     ! Column or surface Jacobian
    HEADR3 = '! Jacobian spectrum for ' // TRIM ( JAC(IJAC)%COD )
  ELSE                               ! Profile level perturbation
    HEADR3 = '! Jacobian spectrum for ' // TRIM ( JAC(IJAC)%COD ) // &
             ' perturbed at altitude =' // TRIM ( C9REAL ( JAC(IJAC)%HGT ) ) &
             // ' [km]'
  END IF
!
! Conversion factor from [cm-1] to [GHz] if required
  IF ( GHZFLG ) THEN
    GHZFAC = 1.0D0 / GHZ2CM 
    IGHFAC = -1
  ELSE
    GHZFAC = 1.0D0
    IGHFAC = 1
  END IF
!
  IF ( SPC(ISPC)%LAB .EQ. '' ) THEN
    LABEL = '''' // TYPTXT // ''''
  ELSE
    LABEL = '''' // SPC(ISPC)%LAB // ''''
  END IF
!
  IF ( BINFLG ) THEN
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR1
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR2
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR3
    WRITE ( LUN, IOSTAT=IOS, ERR=900 ) &
      SPC(ISPC)%NPT*IGHFAC, SPC(ISPC)%WNL*GHZFAC, SPC(ISPC)%WNR*GHZFAC, &
      SPC(ISPC)%WNU*GHZFAC, LABEL
  ELSE
    WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR1
    WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR2
    WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR3
    WRITE ( LUN, '(I7, 3F18.10, 1X, A)', IOSTAT=IOS, ERR=900) &
      SPC(ISPC)%NPT*IGHFAC, SPC(ISPC)%WNL*GHZFAC, SPC(ISPC)%WNR*GHZFAC, &
      SPC(ISPC)%WNU*GHZFAC, LABEL
  END IF
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-WRTHDR: Error writing header records. IOSTAT=', IOS
!
END SUBROUTINE WRTHDR
END MODULE WRTHDR_SUB

