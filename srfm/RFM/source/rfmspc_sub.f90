MODULE RFMSPC_SUB
CONTAINS
SUBROUTINE RFMSPC ( ISPC, FAIL, ERRMSG ) 
!
! VERSION
!   18DEC23 AD Bug#43 - Pass NEWSPC rather than IWID to SPCFIN
!   24JUN19 AD Remove CIAFLG - always call SPCCIA. Checked.
!   03MAR19 AD Bug#17 - Add NOGRID argument to SPCGRD
!   04FEB19 AD Remove subroutine SPCLOS - now part of SPCOUT
!   29JAN18 AD F90 original.
!
! DESCRIPTION
!   Calculations for one spectral range
!   Called by RFM for each spectral range
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
!
! SUBROUTINES
    USE SPCCIA_SUB ! Calculate the Collision-induced absorption
    USE SPCCTM_SUB ! Calculate continuum absorption 
    USE SPCDAL_SUB ! Deallocate spectral range level pointers
    USE SPCFIN_SUB ! Perform a fine pass over the wavenumber grid
    USE SPCFLX_SUB ! Spectral flux calculation
    USE SPCFOV_SUB ! Convolve spectra with FOV
    USE SPCFUL_SUB ! Initialise full grid
    USE SPCGRD_SUB ! Sets up fine-resolution grid for each widemesh interval
    USE SPCILS_SUB ! Convolve spectra with ILS
    USE SPCINI_SUB ! Initialise widemesh grid
    USE SPCINT_SUB ! Interpolate spectra to regular grid
    USE SPCJAC_SUB ! Calculate Jacobians
    USE SPCLUT_SUB ! Calculate the absorption using external Abs.Coeff. Tables.
    USE SPCOUT_SUB ! Write spectral output data
    USE SPCRAD_SUB ! Radiative transfer calculation
    USE SPCREX_SUB ! Calculate the contribution of Rayleigh extinction
    USE SPCSVD_SUB ! Calculate absorption using SVD-compressed Look-Up Tables
    USE SPCTAB_SUB ! Write TABulated absorption coefficients
    USE SPCWID_SUB ! Perform wide mesh absorption calculations
    USE SPCWNG_SUB ! Interpolate wide mesh absorption across fine mesh
    USE SPCXSC_SUB ! Absorption of x/s molecules
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: ISPC   ! Index of spectral range
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error occurs
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL     :: NEWSPC ! T=start of new spectral range
    LOGICAL     :: NOGRID ! T=no grid pts in current widemesh interval
    INTEGER(I4) :: IWID   ! Counter for widemesh intervals
    INTEGER(I4) :: NWID   ! Total number of widemesh intervals
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  CALL SPCFUL ( ISPC, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
  CALL SPCINI ( ISPC, NWID, FAIL, ERRMSG )  
  IF ( FAIL ) RETURN
! 
  IF ( .NOT. SHHFLG ) WRITE (*,*) 'I-RFMSPC: Widemesh calculation ...'
  CALL SPCWID ( FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
  CALL SPCCTM  
!
  IF ( .NOT. SHHFLG ) WRITE (*,*) 'I-RFMSPC: Finemesh calculation ...'
  DO IWID = 1, NWID
    IF ( .NOT. SHHFLG ) WRITE (*,*) IWID, NWID
    NEWSPC = .TRUE.
    CALL SPCGRD ( IWID, NOGRID )
    IF ( NOGRID ) CYCLE
    CALL SPCWNG ( IWID ) 
!
    CALL SPCCIA
!
    IF ( REXFLG ) CALL SPCREX
!
    CALL SPCXSC ( IWID ) 
!
    CALL SPCFIN ( NEWSPC, FAIL, ERRMSG )  ! Resets NEWSPC to FALSE
    IF ( FAIL ) RETURN
!
    IF ( SVDFLG ) CALL SPCSVD
!
    IF ( LUTFLG ) THEN
       CALL SPCLUT ( FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    END IF   
!
    IF ( TABFLG ) THEN
      CALL SPCTAB ( ISPC, IWID, NWID, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    ELSE IF ( FLXFLG ) THEN
      CALL SPCFLX
    ELSE
      CALL SPCRAD 
    END IF
  END DO
!
  IF ( JACFLG ) CALL SPCJAC    ! SPCJAC has to be before SPCFOV
!
  IF ( FOVFLG ) CALL SPCFOV
!
  IF ( AVGFLG .OR. ILSFLG ) THEN
    CALL SPCILS ( ISPC )
  ELSE
    CALL SPCINT 
  END IF
!
  CALL SPCOUT ( ISPC, FAIL, ERRMSG )
!
  CALL SPCDAL 
!
END SUBROUTINE RFMSPC
END MODULE RFMSPC_SUB
