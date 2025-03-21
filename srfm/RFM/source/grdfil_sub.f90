MODULE GRDFIL_SUB
CONTAINS
SUBROUTINE GRDFIL ( NAMGRD, FAIL, ERRMSG )
!
! VERSION
!   14NOV23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Open .grd file and check contents
!   Called by DRVGRD and GRDDEF for each file specified in *GRD section
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE SPCCOM_DAT ! Spectral range data
!
! SUBROUTINES
    USE ADDGFL_SUB ! Add irreg. grid file to GFLCOM
    USE REAGRD_SUB ! Read irregular grid data from .grd file into GRDCOM
    USE REASPC_SUB ! Read irregular grid data from .spc file into GRDCOM
    USE SPCRNG_SUB ! Extract spectral range/resln from file
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMGRD ! Name of .grd file
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)  :: IGFL ! Index of .grd file in GFLCOM
    INTEGER(I4)  :: ISPC ! Spectral range counter (1:NSPC)
    INTEGER(I4)  :: NPT  ! No.points in grid file (dummy)
    CHARACTER(3) :: TYP  ! Type of grid file ('GRD', 'SPC' or 'BIN')
    REAL(R8)     :: WN1  ! Wavenumber [cm-1] of first grid pt (lowest Wno).
    REAL(R8)     :: WN2  ! Wavenumber [cm-1] of last grid pt (highest Wno).
    REAL(R8)     :: WND  ! Nominal Wavenumber spacing of grid (or 0)
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL WRTLOG ( 'I-GRDFIL: Opening .grd file: ' // NAMGRD )
!
  CALL SPCRNG ( NAMGRD, .FALSE., TYP, NPT, WN1, WN2, WND, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN
!
! Check if grid completely spans any selected output spectral range 
! Only check spec ranges with regular output grids - those with irreg output
! grids will already have grid file assigned from *SPC section
! NB: same grid may be usable for more than one spectral range
  IGFL = 0
!         
  DO ISPC = 1, NSPC
    IF ( SPC(ISPC)%WNR .EQ. 0.0D0 .OR. &      ! irregular output grid or
         SPC(ISPC)%WNL .LT. WN1   .OR. &      ! outside range
         SPC(ISPC)%WNU .GT. WN2         ) CYCLE
!
! Check no irreg grid already assigned to this spectral range
    IF ( SPC(ISPC)%IGD .NE. 0 ) THEN
      FAIL = .TRUE.                         
      ERRMSG = 'F-GRDFIL: Another Grid already assigned ' // &
               'for spectral range: ' // SPC(ISPC)%LAB
      RETURN
    END IF
!
! If first time this file is required, add to list
    IF ( IGFL .EQ. 0 ) CALL ADDGFL ( NAMGRD, TYP, WN1, WN2, IGFL ) 
!
! Assign this grid file to spectral range
    SPC(ISPC)%IGD = IGFL
    CALL WRTLOG ( 'I-GRDFIL: GRD assigned to spectral range: '//SPC(ISPC)%LAB )
!
! Read in all grid points required for this spectral range to check that file
! can be read OK and that there is sufficient array space to store grid
    IF ( TYP .EQ. 'GRD' ) THEN           ! .grd file
      CALL REAGRD ( NAMGRD, SPC(ISPC)%WXL, SPC(ISPC)%WXU, FAIL, ERRMSG )
    ELSE                               ! spectral file
      CALL REASPC ( NAMGRD, (TYP .EQ. 'BIN'), &
                    SPC(ISPC)%WXL, SPC(ISPC)%WXU, FAIL, ERRMSG )  
    END IF
    IF ( FAIL ) RETURN
  END DO
!
  IF ( IGFL .EQ. 0 ) &
    CALL WRTLOG ( 'W-GRDFIL: File ignored - grid not applicable' )
!
END SUBROUTINE GRDFIL
END MODULE GRDFIL_SUB
