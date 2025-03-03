MODULE CHKFLX_SUB
CONTAINS
SUBROUTINE CHKFLX ( HGTTST, FAIL, ERRMSG )
!
! VERSION
!   19JAN24 AD Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Check output levels for Flux calculations
!   Called by TANCHK for each path with FLX flag enabled. 
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT, ONLY: HGTSFC, HGTTOA ! Surface, top-of-atmosphere altitudes
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
!
  IMPLICIT NONE 
!
! ARGUMENTS
    REAL(R4),      INTENT(IN)  :: HGTTST ! Value of level to be checked
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  FAIL = .TRUE.
!
  IF ( HGTTST .LT. HGTSFC ) THEN
    ERRMSG = 'F-CHKFLX: Level = ' // TRIM ( C9REAL(HGTTST) ) // &
         ' km below base of atmosphere (' // TRIM ( C9REAL(HGTSFC) ) // ' km)'
  ELSE IF ( HGTTST .GT. HGTTOA ) THEN
    ERRMSG = 'F-CHKFLX: Level = ' // TRIM ( C9REAL(HGTTST) ) // & 
         ' km above top of atmosphere (' // TRIM ( C9REAL(HGTTOA) ) //  ' km)'
  ELSE
    FAIL = .FALSE.
  END IF
!
END SUBROUTINE CHKFLX
END MODULE CHKFLX_SUB
