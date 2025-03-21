MODULE CHKABS_SUB
CONTAINS
SUBROUTINE CHKABS ( FAIL, ERRMSG )
!
! VERSION
!   25JUN24 AD Checked.
!   30JUN23 AD Bug#40 Check if GAS%CIA is set TRUE as well.
!   20DEC17 AD F90 conversion of abschk.for. Checked.
!
! DESCRIPTION
!   Check information available for all absorbers
!   Called once by DRVCHK.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE LFLCOM_DAT ! Look-Up Table files
    USE SFLCOM_DAT ! SVD-compressed LUT files
    USE SPCCOM_DAT ! Spectral range data
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL     :: LSOURC ! T=spec.data source defined by LUT or SVD
    INTEGER(I4) :: IGAS   ! Absorber counter
    INTEGER(I4) :: ISPC   ! Spectran range counter
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO IGAS = 1, NGAS
    IF ( GAS(IGAS)%COD .EQ. 'air' ) CYCLE     ! No absorber data reqd for air
    DO ISPC = 1, NSPC
      LSOURC = .FALSE.
      IF ( NLFL .GT. 0 ) THEN                      ! Using (some) LUTs
        LSOURC = IDXLFL(IGAS,ISPC) .GT. 0 
      END IF
      IF ( NSFL .GT. 0 ) THEN
        IF ( IDXSFL(IGAS,ISPC) .GT. 0 ) THEN
          IF ( LSOURC ) THEN
            FAIL = .TRUE.
            ERRMSG = 'F-CHKABS: Both LUT and SVD defined for ' // &
                     TRIM ( GAS(IGAS)%COD ) // ' for spectral range ' // &
                     TRIM ( SPC(ISPC)%LAB ) 
            RETURN
          END IF
          LSOURC = .TRUE.
        END IF
      END IF
      IF ( .NOT. ( LSOURC .OR. GAS(IGAS)%HIT .OR. GAS(IGAS)%XSC &
                   .OR. GAS(IGAS)%CIA ) ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-CHKABS: No spectroscopic data for ' // &
                 TRIM ( GAS(IGAS)%COD ) // ' for spectral range ' // &
                 TRIM ( SPC(ISPC)%LAB )
        RETURN
      END IF
    END DO
  END DO 
!
! Normal exit
  FAIL = .FALSE.
!
END SUBROUTINE CHKABS
END MODULE CHKABS_SUB
