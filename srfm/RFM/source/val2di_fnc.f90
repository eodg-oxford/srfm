MODULE VAL2DI_FNC
CONTAINS
REAL(R4) PURE FUNCTION VAL2DI ( XTAB, XINT, YTAB, YINT, ZTAB, LOGINT )
!
! VERSION
!   06JAN24 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Interpolate value from 2D array
!   General purpose function
!   2D version of VAL1DI.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE INTRVL_SUB ! Calculate interpolation interval
    USE VAL1DI_GEN ! Interpolate value from 1D array
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: XTAB(:)   ! List of x-axis tabulation coordinates
    REAL(R4), INTENT(IN) :: XINT      ! x-axis interpolation coordinate
    REAL(R4), INTENT(IN) :: YTAB(:)   ! List of y-axis tabulation coordinates
    REAL(R4), INTENT(IN) :: YINT      ! y-axis interpolation coordinate
    REAL(R4), INTENT(IN) :: ZTAB(:,:) ! Interpolated array
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: LOGINT  ! TRUE=interpolate linearly in Log(YTAB)
!
! LOCAL VARIABLES
    LOGICAL     :: LINT   ! T=Log interpolation, F=linear interpolation
    INTEGER(I4) :: IX     ! Index of XTAB element below XINT
    INTEGER(I4) :: IY     ! Index of YTAB element below YINT
    INTEGER(I4) :: NXTAB  ! Size of XTAB axis
    INTEGER(I4) :: NYTAB  ! Size of YTAB axis
    REAL(R4)    :: DX     ! XINT as fraction of XTAB interval
    REAL(R4)    :: DY     ! YINT as fraction of YTAB interval
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NXTAB = SIZE ( XTAB ) 
  NYTAB = SIZE ( YTAB ) 
!
! Special case of just a single element in tabulated values
  IF ( NXTAB .EQ. 1 .AND. NYTAB .EQ. 1 ) THEN
    VAL2DI = ZTAB(1,1)
    RETURN
  END IF
!
  IF ( PRESENT ( LOGINT ) ) THEN
    LINT = LOGINT
  ELSE
    LINT = .FALSE.
  END IF
!
! If either X or Y dimension is size 1 treat as 1D linear interpolation
  IF ( NXTAB .EQ. 1 ) THEN
    VAL2DI = VAL1DI ( YTAB, YINT, ZTAB(1,:), LINT ) 
    RETURN
  ELSE IF ( NYTAB .EQ. 1 ) THEN
    VAL2DI = VAL1DI ( XTAB, XINT, ZTAB(:,1), LINT ) 
    RETURN
  END IF

  CALL INTRVL ( XTAB, XINT, IX, DX, LINT ) 
  CALL INTRVL ( YTAB, YINT, IY, DY, LINT ) 
!
  IF ( LINT ) THEN
    VAL2DI = EXP ( &
         (1.0-DX) * (1.0-DY) * LOG ( MAX ( TINY(1.0), ZTAB(IX,IY)     ) ) + &
         (1.0-DX) *      DY  * LOG ( MAX ( TINY(1.0), ZTAB(IX,IY+1)   ) ) + &
              DX  * (1.0-DY) * LOG ( MAX ( TINY(1.0), ZTAB(IX+1,IY)   ) ) + &
              DX  *      DY  * LOG ( MAX ( TINY(1.0), ZTAB(IX+1,IY+1) ) ) )
  ELSE 
    VAL2DI = (1.0-DX) * (1.0-DY) * ZTAB(IX,IY)   + &
             (1.0-DX) *      DY  * ZTAB(IX,IY+1) + &
                  DX  * (1.0-DY) * ZTAB(IX+1,IY) + &
                  DX  *      DY  * ZTAB(IX+1,IY+1) 
  END IF
!
END FUNCTION VAL2DI
END MODULE VAL2DI_FNC
