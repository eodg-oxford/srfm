MODULE FLGCOM_DAT
!
! VERSION
!   15MAR23 AD Add C32 flag.
!   24JUN19 AD Remove CIA flag. Checked.
!   19DEC17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   RFM option flags
!   Loaded by DRVFLG or presence of optional sections in driver table.
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL VARIABLES
    LOGICAL :: ABSFLG = .FALSE. ! T = output Absorption spectra
    LOGICAL :: AVGFLG = .FALSE. ! T = spectrally average output
    LOGICAL :: BBTFLG = .FALSE. ! T = output Brightness Temperature spectra
    LOGICAL :: BFXFLG = .FALSE. ! T = Source Fn varies with layer optical depth
    LOGICAL :: BINFLG = .FALSE. ! T = write output spectra as binary files
    LOGICAL :: C32FLG = .FALSE. ! T = use (old) MT_CKD v3.2 H2O continuum model
    LOGICAL :: CHIFLG = .FALSE. ! T = use Chi-Factor for any CO2 lines
    LOGICAL :: CLCFLG = .FALSE. ! T = explicit LBL calc for each path
    LOGICAL :: COOFLG = .FALSE. ! T = Calculate cooling rates
    LOGICAL :: CTMFLG = .FALSE. ! T = use continuum data
    LOGICAL :: DBLFLG = .FALSE. ! T = write output spectra as Double Precision
    LOGICAL :: FINFLG = .FALSE. ! T = user-defined finemesh resolution
    LOGICAL :: FLXFLG = .FALSE. ! T = solid angle integration for flux calc
    LOGICAL :: FOVFLG = .FALSE. ! T = apply FOV convolution
    LOGICAL :: FVZFLG = .FALSE. ! T = FOV fn defined as refracted tan. heights
    LOGICAL :: GEOFLG = .FALSE. ! T = use geometric (ie non-refrac) ray paths
    LOGICAL :: GHZFLG = .FALSE. ! T = use GHz for spectral axis, F=wavenumber
    LOGICAL :: GRAFLG = .FALSE. ! T = use horiontal gradients
    LOGICAL :: GRDFLG = .FALSE. ! T = use irregular grid
    LOGICAL :: HOMFLG = .FALSE. ! T = use homogeneous path
    LOGICAL :: HYDFLG = .FALSE. ! T = use hydros. equlbm for plane par.CG paths
    LOGICAL :: ILSFLG = .FALSE. ! T = use instrument lineshape convolution
    LOGICAL :: JACFLG = .FALSE. ! T = calculate Jacobian spectra
    LOGICAL :: JTPFLG = .FALSE. ! T = Jacobians for tan.pt. pertubrations only
    LOGICAL :: LAYFLG = .FALSE. ! T = use sub-layering
    LOGICAL :: LEVFLG = .FALSE. ! T = spectra at intermediate output levels
    LOGICAL :: LINFLG = .FALSE. ! T = Assume VMR varies linearly with altitude
    LOGICAL :: LOSFLG = .FALSE. ! T = Calculate elev. pointing Jacobian spectra
    LOGICAL :: LUNFLG = .FALSE. ! T = Reuse same LUN for all output spectra
    LOGICAL :: LUTFLG = .FALSE. ! T = Use Look-Up Tables for absorption coeffs.
    LOGICAL :: MIXFLG = .FALSE. ! T = use line-mixing
    LOGICAL :: MTXFLG = .FALSE. ! T = inter-level matrix of flux calculations
    LOGICAL :: NADFLG = .FALSE. ! T = nadir-viewing
    LOGICAL :: NEWFLG = .FALSE. ! T = open new files with STATUS='NEW'
    LOGICAL :: NTEFLG = .FALSE. ! T = non-LTE calculations
    LOGICAL :: OBSFLG = .FALSE. ! T = observer located within atmosphere
    LOGICAL :: OPTFLG = .FALSE. ! T = output Optical Depth spectra
    LOGICAL :: PRFFLG = .FALSE. ! T = output internal atmospheric profile
    LOGICAL :: PTHFLG = .FALSE. ! T = output Path diagnostics
    LOGICAL :: QADFLG = .FALSE. ! T = use simple quadratic fit to line wings
    LOGICAL :: RADFLG = .FALSE. ! T = output Radiance spectra
    LOGICAL :: REJFLG = .FALSE. ! T = reject weak lines
    LOGICAL :: REXFLG = .FALSE. ! T = implement Rayleigh extinction
    LOGICAL :: RJTFLG = .FALSE. ! T = output Rayleigh-Jeans Bright.Temp spectra
    LOGICAL :: SFCFLG = .FALSE. ! T = Allow for opaque surface under atmosphere
    LOGICAL :: SHHFLG = .FALSE. ! T = suppress progress messages
    LOGICAL :: SHPFLG = .FALSE. ! T = use non-Voigt line shapes
    LOGICAL :: SVDFLG = .FALSE. ! T = SVD-compressed LUTs
    LOGICAL :: TABFLG = .FALSE. ! T = output tabulated absorption coefficients
    LOGICAL :: TRAFLG = .FALSE. ! T = output Transmittance spectra
    LOGICAL :: VRTFLG = .FALSE. ! T = flux-type calcs without hemisph.integ.
    LOGICAL :: VVWFLG = .FALSE. ! T = apply Van Vleck-Weisskopf correction
    LOGICAL :: WIDFLG = .FALSE. ! T = output Wide Mesh Calc diagnostics
    LOGICAL :: ZENFLG = .FALSE. ! T = zenith-viewing
!
END MODULE FLGCOM_DAT
