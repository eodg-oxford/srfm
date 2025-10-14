MODULE CHKHOM_SUB
CONTAINS
SUBROUTINE CHKHOM ( LENTST, FAIL, ERRMSG )
!
! VERSION
!   20JAN24 AD Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Check path length for HOMogeneous path calculation
!   Called by TANCHK if HOM flag enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
!
  IMPLICIT NONE 
!
! ARGUMENTS
    REAL(R4),      INTENT(IN)  :: LENTST ! Value of level to be checked
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  IF ( LENTST .LE. 0.0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-CHKHOM: Specified Homog. Path Length .LE. 0, value =' &
             // C9REAL(LENTST)
  ELSE
    FAIL = .FALSE.
  END IF
!
END SUBROUTINE CHKHOM
END MODULE CHKHOM_SUB
