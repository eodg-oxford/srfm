MODULE ADDPTH_SUB
CONTAINS
SUBROUTINE ADDPTH ( PTHNEW, PTH, JDX )
!
! VERSION
!   03DEC21 AD Further modifications. Checked.
!   20SEP21 AD Rewritten/redefined
!
! DESCRIPTION
!   Add perturbed path segment
!   Called by JACPTH
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT, ONLY: CLCFLG ! T = explicit LBL calc for each path
    USE JACCOM_DAT, ONLY: JDXTEM, JDXPRE ! Indices of Tem,Pre Jacobians
    USE PTHCOM_DAT, ONLY: PTHTYP ! Path data type
!
! SUBROUTINES
    USE ICLPTH_FNC ! Index of corresponding (new) calculated path
!
  IMPLICIT NONE
!
! ARGUMENTS
    TYPE(PTHTYP), INTENT(IN)  :: PTHNEW ! Perturbed path
    TYPE(PTHTYP), ALLOCATABLE, &
                INTENT(INOUT) :: PTH(:) ! Set of unperturbed paths
    INTEGER(I4),  INTENT(IN)  :: JDX    ! Index of Jacobian type
!
! LOCAL VARIABLES
    INTEGER(I4)               :: NPTH      ! No. paths in PTH
    TYPE(PTHTYP), ALLOCATABLE :: PTHSAV(:) ! Saved ORG during reallocation
!
! EXECUTABLE CODE ------------------------------------------------------------
!
  NPTH = SIZE ( PTH ) 
  CALL MOVE_ALLOC ( PTH, PTHSAV ) 
  NPTH = NPTH + 1
  ALLOCATE ( PTH(NPTH) )
  PTH(1:NPTH-1) = PTHSAV
  PTH(NPTH) = PTHNEW
! Only create a new calculated path if new Tem or Pre Jacobian path, or NEWCLC
! otherwise keep %ICL pointing to same CLC path as unperturbed segment
  IF ( JDX .EQ. JDXTEM .OR. JDX .EQ. JDXPRE .OR. CLCFLG ) & 
    PTH(NPTH)%ICL = ICLPTH ( PTHNEW, .TRUE. )
!
END SUBROUTINE ADDPTH
END MODULE ADDPTH_SUB
