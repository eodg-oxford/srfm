MODULE INTRVL_SUB
CONTAINS
PURE SUBROUTINE INTRVL ( XTAB, X, IX, DX, LOGINT )
!
! VERSION
!   31DEC23 AD Checked.
!   01JUL17 AD Original. Checked.
! 
! DESCRIPTION    
!   Calculate interpolation interval
!   General purpose module
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE IBRAKT_GEN ! Lower index of array interpolation
!                  
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4),    INTENT(IN)  :: XTAB(:) ! Data Array
    REAL(R4),    INTENT(IN)  :: X       ! Value to be bracketted
    INTEGER(I4), INTENT(OUT) :: IX      ! Index of array .LE. X
    REAL(R4),    INTENT(OUT) :: DX      ! Fraction of interval of X in IX:IX+1
    LOGICAL, OPTIONAL, &
                 INTENT(IN)  :: LOGINT  ! T=Logarithmic interpolation
!
! LOCAL VARIABLES
    LOGICAL     :: LINT ! T=log interpolation, F=linear 
    INTEGER(I4) :: JX   ! IX+1
    INTEGER(I4) :: N    ! Size of array XTAB
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  N = SIZE ( XTAB ) 
  IF ( N .EQ. 1 ) THEN
    IX = 1
    DX = 0.0
    RETURN
  END IF
!
  IX = IBRAKT ( XTAB, X ) 
  IF ( IX .LT. 1 ) THEN
    IX = 1
    DX = 0.0
  ELSE IF ( IX .EQ. N ) THEN
    IX = N - 1
    DX = 1.0
  ELSE
    JX = IX + 1
    IF ( PRESENT ( LOGINT ) ) THEN
      LINT = LOGINT .AND. XTAB(IX) .GT. TINY(1.0) &
                    .AND. XTAB(JX) .GT. TINY(1.0) 
    ELSE
      LINT = .FALSE.
    END IF
!
    IF ( LINT ) THEN
      DX = ( LOG(X) - LOG(XTAB(IX)) ) / ( LOG(XTAB(JX)) - LOG(XTAB(IX)) ) 
    ELSE
      DX = ( X - XTAB(IX) ) / ( XTAB(JX) - XTAB(IX) ) 
    END IF
  END IF
!
END SUBROUTINE INTRVL
END MODULE INTRVL_SUB
