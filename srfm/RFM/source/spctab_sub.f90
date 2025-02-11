MODULE SPCTAB_SUB
CONTAINS
SUBROUTINE SPCTAB ( ISPC, IWID, NWID, FAIL, ERRMSG )
!
! VERSION
!   03APR23 AD Checked.
!   08NOV17 AD F90 original. Checked.
!
! DESCRIPTION
!   Write TABulated absorption coefficients
!   Called by RFMSPC for each widemesh interval if TAB flag enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE TABCOM_DAT ! Axes for tabulated absorption coefficients
    USE NAMCOM_DAT, ONLY: TABNAM ! Name of .tab output file
!
! SUBROUTINES
    USE OPNOUT_SUB ! Open spectral output files for current spectral range
    USE TABHDR_SUB ! Write TAB file header records
    USE TABWRT_SUB ! Write RFM tabulated absorption coefficients
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: ISPC   ! Spectral range number
    INTEGER(I4),   INTENT(IN)  :: IWID   ! Widemesh interval#
    INTEGER(I4),   INTENT(IN)  :: NWID   ! No.Widemesh intervals
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error occurs
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: ITAB ! Counter for .tab files/molecules
    INTEGER(I4) :: LUN  ! Logical Unit Number for .tab files
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! On first call in each spectral range, open .tab files and write headers
  IF ( IWID .EQ. 1 ) THEN
    LUN = LUNTAB
    DO ITAB = 1, NTAB
      CALL OPNOUT ( LUN, TABNAM, FAIL, ERRMSG, IGAS=ITAB, ISPC=ISPC )
      IF ( FAIL ) RETURN
      CALL TABHDR ( LUN, ITAB, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      LUN = LUN + 1
    END DO
  END IF
!
! Output abs.coeff. data for current widemesh interval 
  CALL TABWRT ( FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
! On last call in each spectral range, close .tab files
  IF ( IWID .EQ. NWID ) THEN
    LUN = LUNTAB
    DO ITAB = 1, NTAB
      CLOSE ( LUN )
    END DO
  END IF
!
END SUBROUTINE SPCTAB
END MODULE SPCTAB_SUB

