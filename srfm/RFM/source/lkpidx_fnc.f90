MODULE LKPIDX_FNC
CONTAINS
REAL(R4) PURE FUNCTION LKPIDX ( XIDX, YTAB )
!
! VERSION
!   10APR23 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Interpolate by real array index
!   General purpose function
!   Returns a value 0 if XIDX < 1 or XIDX > size(YTAB)
!   Unlike other interpolation functions, this one assumes that the 
!   interpolation interval has already been calculated
!
! VARIABLE KINDS
    USE KIND_DAT
! 
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: XIDX    ! Interpolation coordinates
    REAL(R4), INTENT(IN) :: YTAB(:) ! List of tabulated data values at XTAB
!
! LOCAL VARIABLES
    INTEGER(I4) :: IX   ! Index below XIDX
    INTEGER(I4) :: NTAB ! Size of YTAB array
    REAL(R4)    :: DX   ! Fraction of XINT above IX
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NTAB = SIZE ( YTAB )
!
  IF ( XIDX .LT. 1.0 .OR. XIDX .GT. NTAB ) THEN
    LKPIDX = 0.0
  ELSE
    IX = MIN ( INT ( XIDX ), NTAB - 1 ) 
    DX = XIDX - FLOAT ( IX ) 
    LKPIDX = (1.0-DX) * YTAB(IX) + DX * YTAB(IX+1) 
  END IF
!
END FUNCTION LKPIDX
END MODULE LKPIDX_FNC
