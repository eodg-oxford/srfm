MODULE ILSSPC_SUB
CONTAINS
SUBROUTINE ILSSPC ( FAIL, ERRMSG )
!
! VERSION
!   15NOV23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Assign ILS Functions to each spectral range
!   Called once by DRVILS
! 
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ILSCOM_DAT ! Instrument Lineshape data
    USE SPCCOM_DAT ! Spectral range data
!
  IMPLICIT NONE 
!
! ARGUMENTS 
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IILS   ! Counter for ILS Functions
    INTEGER(I4) :: ISPC   ! Counter for spectral ranges for RFM calcs.
    REAL(R8)    :: WNOTOL ! Tolerance [cm-1] for matching wavenumbers
!
! EXECUTABLE CODE ------------------------------------------------------------
!
  DO ISPC = 1, NSPC
!
! Check that one, and only one, ILS function can be assigned to every spc.range
    SPC(ISPC)%ILS = 0                ! no ILS fn assigned yet
    DO IILS = 1, NILS 
      IF ( IILS .EQ. IDFILS ) CYCLE  ! Skip default ILS for now
      WNOTOL = 0.1D0 * ILS(IILS)%PTD
      IF ( SPC(ISPC)%WNL .GE. ILS(IILS)%WNL - WNOTOL .AND. &
           SPC(ISPC)%WNU .LE. ILS(IILS)%WNU + WNOTOL      )  THEN ! match found
        IF ( SPC(ISPC)%ILS .NE. 0 ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-ILSSPC: Ambiguous: at least 2 ILS Fns.' // &
                   ' for spectral range: ' // SPC(ISPC)%LAB
          RETURN
        END IF
        SPC(ISPC)%ILS = IILS
        EXIT
      END IF
    END DO
!
! If still not assigned, use any default ILS fn
    IF ( SPC(ISPC)%ILS .EQ. 0 ) SPC(ISPC)%ILS = IDFILS  ! IDFILS=0 if no default
!
! If still unassigned, fatal error
    IF ( SPC(ISPC)%ILS .EQ. 0 ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-ILSSPC: No ILS Fn supplied to cover spectral range: ' // &
               SPC(ISPC)%LAB
      RETURN
    END IF
!
    IILS = SPC(ISPC)%ILS
    SPC(ISPC)%WXL = SPC(ISPC)%WNL + ILS(IILS)%PT1
    SPC(ISPC)%WXU = SPC(ISPC)%WNU + ILS(IILS)%PT2
  END DO      ! end loop over spectral ranges
!
! Normal exit
  FAIL = .FALSE.
!
END SUBROUTINE ILSSPC
END MODULE ILSSPC_SUB
