MODULE SNELL_SUB 
CONTAINS
PURE SUBROUTINE SNELL ( INPDCC, DRAD, DRFR, DSH, DCC, DSIN, DRAT )
!
! VERSION
!   12JAN24 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Apply Snell's Law
!   Called by RAYSUM
!   Calculate snells law constant R = - r/(n/n'), =0 if refractivity=0
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,  INTENT(IN)    :: INPDCC ! T=DCC specified, F=DSIN spec
    REAL(R8), INTENT(IN)    :: DRAD   ! Rad.of DZn from centre earth
    REAL(R8), INTENT(IN)    :: DRFR   ! Refractivity
    REAL(R8), INTENT(IN)    :: DSH    ! Dens or Refr Scale height [km]
    REAL(R8), INTENT(INOUT) :: DCC    ! n.r.sin(theta)
    REAL(R8), INTENT(INOUT) :: DSIN   ! Sine of zenith angle 
    REAL(R8), INTENT(OUT)   :: DRAT   ! R/(n/n') (R=DRAD)
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DRAT = DRAD * DRFR / ( DSH * (1.D0 + DRFR) ) 
  IF ( INPDCC ) THEN
    DSIN = DCC / ( (DRFR + 1.D0) * DRAD )
  ELSE
    DCC = DSIN * DRAD * ( DRFR + 1.D0 )
  END IF
!
END SUBROUTINE SNELL
END MODULE SNELL_SUB
