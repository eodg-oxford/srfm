MODULE SPCFUL_SUB
CONTAINS
SUBROUTINE SPCFUL ( ISPC, FAIL, ERRMSG ) 
!
! VERSION
!   12JUL24 AD Checked.
!   07JUN19 AD Bug#21 Adjust WNRFUL to match irreg.grid. Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Initialise full grid 
!   Called by RFMSPC at the start of each spectral range.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE FULCOM_DAT ! Full grid data
    USE GFLCOM_DAT ! Irregular grid files
    USE GRDCOM_DAT ! Irregular grid
    USE SPCCOM_DAT ! Spectral range data
    USE FINCOM_DAT, ONLY: NOMFIN ! Finemesh data
    USE TANCOM_DAT, ONLY: MTAN   ! No.calculated tangent paths
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
    USE REAGRD_SUB ! Read irregular grid data from .grd file into GRDCOM
    USE REASPC_SUB ! Read irregular grid data from .spc file into GRDCOM
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: ISPC   ! Spectral range index
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IFUL ! Counter for full grid
    INTEGER(I4) :: IGFL ! Index of irregular grid file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( SPC(ISPC)%WNR .EQ. 0.0D0 ) THEN
    WNRFUL = 1.0D0 / NOMFIN
  ELSE IF ( ILSFLG .OR. AVGFLG ) THEN  
    WNRFUL = SPC(ISPC)%WNR / NINT ( NOMFIN * SPC(ISPC)%WNR ) 
  ELSE
    WNRFUL = SPC(ISPC)%WNR 
  END IF
!
! WNOFUL etc should be allocated for ISPC > 1
  IF ( ALLOCATED ( WNOFUL ) ) DEALLOCATE ( OPTFUL, RADFUL, TRAFUL, WNOFUL )
  IF ( ALLOCATED ( COOFUL ) ) DEALLOCATE ( COOFUL ) ! only if COOFLG=T
!
  IRRFUL = SPC(ISPC)%IGD .GT. 0
  IF ( IRRFUL ) THEN               ! irregular grid for this range
    IGFL = SPC(ISPC)%IGD
    IF ( GFL(IGFL)%GRD ) THEN
      CALL REAGRD ( GFL(IGFL)%NAM, SPC(ISPC)%WXL, SPC(ISPC)%WXU, FAIL, ERRMSG )
    ELSE
      CALL REASPC ( GFL(IGFL)%NAM, GFL(IGFL)%BIN, SPC(ISPC)%WXL, SPC(ISPC)%WXU,&
                    FAIL, ERRMSG )
    END IF
    IF ( FAIL ) RETURN 
! Bug#21 .grd files have their own value of WNRGRD which may be different 
! from WNRFUL (established with NOMFIN), so switch if necessary
    IF ( ABS ( WNRFUL - WNRGRD ) .GT. 1.0E-6 * WNRFUL ) THEN
      CALL WRTLOG ( 'W-SPCFUL: Changing fine grid from ' // &
                    TRIM ( C9REAL(WNRFUL) ) // ' to ' // &
                    TRIM ( C9REAL(WNRGRD) ) // ' cm-1 to match irreg.grd' )
      WNRFUL = WNRGRD
    END IF
    NFUL = NGRD
    ALLOCATE ( WNOFUL(NFUL) ) 
    WNOFUL = GRD(1:NGRD)%WNO
  ELSE 
    NFUL = NINT ( ( SPC(ISPC)%WXU - SPC(ISPC)%WXL ) / WNRFUL ) + 1
    ALLOCATE ( WNOFUL(NFUL) )
    DO IFUL = 1, NFUL
      WNOFUL(IFUL) = SPC(ISPC)%WXL + (IFUL-1)*WNRFUL
    END DO
  END IF
!
  IFUL1 = 0
  IFUL2 = 0
  WNLFUL = WNOFUL(1)
  WNUFUL = WNOFUL(NFUL)
!
  ALLOCATE ( OPTFUL(NFUL,MTAN), RADFUL(NFUL,MTAN), TRAFUL(NFUL,MTAN) )
  OPTFUL = 0.0D0
  RADFUL = 0.0D0
  TRAFUL = 1.0D0
  IF ( COOFLG ) THEN
    ALLOCATE ( COOFUL(NFUL,MTAN) ) 
    COOFUL = 0.0D0
  END IF
!
END SUBROUTINE SPCFUL
END MODULE SPCFUL_SUB
