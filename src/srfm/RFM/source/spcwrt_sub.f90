MODULE SPCWRT_SUB
CONTAINS
SUBROUTINE SPCWRT ( NAMTMP, TYP, ISPC, ITAN, IJAC, ILEV, NPT, IRREG, &
                    WNO, SPC, FAIL, ERRMSG )
! VERSION
!   02MAR23 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Write spectral data file
!   Called by SPCOUT
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
    USE PYOUT_DAT, ONLY: CAPTURE_ENABLED
!
! SUBROUTINES
    USE OPNOUT_SUB ! Open spectral output files for current spectral range
    USE WRTHDR_SUB ! Write header records of output files
    USE WRTSPC_SUB ! Write spectral data to file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMTMP ! Filename template
    CHARACTER(*),  INTENT(IN)  :: TYP    ! Type of spectrum 'ABS','COO',etc
    INTEGER(I4),   INTENT(IN)  :: ISPC   ! Spectral range index
    INTEGER(I4),   INTENT(IN)  :: ITAN   ! Index of output tan.path
    INTEGER(I4),   INTENT(IN)  :: IJAC   ! Index of Jacobian, or 0
    INTEGER(I4),   INTENT(IN)  :: ILEV   ! Index of output level, or 0
    INTEGER(I4),   INTENT(IN)  :: NPT    ! No. points to be written
    LOGICAL,       INTENT(IN)  :: IRREG  ! T=irreg grid, F=regular grid
    REAL(R8),      INTENT(IN)  :: WNO(:) ! Spectral axis [cm-1]
    REAL(R8),      INTENT(IN)  :: SPC(:) ! Output spectrum
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error occurs
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( CAPTURE_ENABLED ) THEN
    FAIL = .FALSE.
    ERRMSG = ''
    RETURN
  END IF
!
  CALL OPNOUT ( LUNTMP, NAMTMP, FAIL, ERRMSG, &
                ISPC=ISPC, ITAN=ITAN, IJAC=IJAC, ILEV=ILEV )
  IF ( FAIL ) RETURN
!
  CALL WRTHDR ( LUNTMP, TYP, IJAC, ILEV, ISPC, ITAN, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
  CALL WRTSPC ( LUNTMP, NPT, IRREG, WNO, SPC, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
  CLOSE ( LUNTMP ) 
!
END SUBROUTINE SPCWRT
END MODULE SPCWRT_SUB
