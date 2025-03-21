MODULE ICLPTH_FNC
CONTAINS
INTEGER(I4) FUNCTION ICLPTH ( PTH, NEWCLC ) 
!
! VERSION
!   20SEP21 AD Original. Checked.
!
! DESCRIPTION
!   Index of corresponding (new) calculated path
!   General purpose function.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated absorption coefficients
    USE FLGCOM_DAT, ONLY: CLCFLG ! T = explicit LBL calc for each path
    USE PTHCOM_DAT, ONLY: PTHTYP ! Path data type
    USE RFMCON_DAT, ONLY: PCGMAX, TCGMAX ! Tolerances for scaled paths
!
  IMPLICIT NONE
!
! ARGUMENTS
    TYPE(PTHTYP),      INTENT(IN) :: PTH    ! Path to be matched
    LOGICAL, OPTIONAL, INTENT(IN) :: NEWCLC ! T=force new calc path
!
! LOCAL VARIABLES
    INTEGER(I4)               :: ICLC      ! Counter for calculated paths
    LOGICAL                   :: ADDNEW    ! T=add new path    
    TYPE(CLCTYP), ALLOCATABLE :: CLCSAV(:) ! Saved CLC during reallocation
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  IF ( CLCFLG ) THEN   ! CLC flag overrides NEWCLC argument
    ADDNEW = .TRUE.
  ELSE IF ( PRESENT ( NEWCLC ) ) THEN
    ADDNEW = NEWCLC
  ELSE
    ADDNEW = .FALSE.
  END IF
!
  IF ( .NOT. ADDNEW ) THEN      ! Look for possible scaled path
    DO ICLC = 1, NCLC
      IF ( PTH%IGS .EQ. CLC(ICLC)%IGS .AND. &
           ABS ( PTH%PRE / CLC(ICLC)%PRE - 1.0 ) .LT. PCGMAX .AND. &
           ABS ( PTH%TEM - CLC(ICLC)%TEM )       .LT. TCGMAX  ) THEN
        ICLPTH = ICLC
        RETURN             ! Found a calc path for scaling
      END IF
    END DO
  END IF
!
! If this point reached, need to add another calculated path
  IF ( ALLOCATED ( CLC ) ) CALL MOVE_ALLOC ( CLC, CLCSAV )
  NCLC = NCLC + 1
  ALLOCATE ( CLC(NCLC) )
  IF ( ALLOCATED ( CLCSAV ) ) CLC(1:NCLC-1) = CLCSAV
!
  CLC(NCLC)%IGS = PTH%IGS
  CLC(NCLC)%AMT = 1.0        ! All calc paths have standard 1 kmol/cm2
  CLC(NCLC)%PPA = PTH%PPA
  CLC(NCLC)%PRE = PTH%PRE
  CLC(NCLC)%TEM = PTH%TEM
!
  ICLPTH = NCLC
!
END FUNCTION ICLPTH
END MODULE ICLPTH_FNC
