MODULE ADDTAN_SUB
CONTAINS
SUBROUTINE ADDTAN ( ITAN, CLCTAN )
!
! VERSION
!   26NOV23 AD Checked.
!   31MAY17 AD Ensure %ITN always defined. Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Add tangent ray path
!   General purpose module.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE TANCOM_DAT ! Tangent path data
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: ITAN   ! Index of tan.path to be copied
    LOGICAL,     INTENT(IN) :: CLCTAN ! T=Add as Calc.tan pth
!
! LOCAL VARIABLES
    TYPE(TANTYP), ALLOCATABLE :: TANSAV(:) ! Saved TAN during reallocation
!
! EXECUTABLE CODE ------------------------------------------------------------
!
  CALL MOVE_ALLOC ( TAN, TANSAV ) 
  MTAN = MTAN + 1
  ALLOCATE ( TAN(MTAN) ) 
  TAN(1:MTAN-1) = TANSAV
  IF ( ITAN .GT. 0 ) THEN
    TAN(MTAN) = TAN(ITAN)
    TAN(MTAN)%ITN = ITAN
  ELSE
    TAN(MTAN)%ITN = 0
  END IF
  TAN(MTAN)%CLC = CLCTAN
!
END SUBROUTINE ADDTAN
END MODULE ADDTAN_SUB
