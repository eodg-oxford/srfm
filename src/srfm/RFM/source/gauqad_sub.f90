MODULE GAUQAD_SUB
CONTAINS
SUBROUTINE GAUQAD ( N, X, W )
!
! VERSION
!   09DEC23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Values and Weights for Gaussian First Moment Quadrature
!   Called by INIQAD
!   Integrating x.f(x).dx in interval x=0:1
!   Numbers taken from Abramowitz and Stegun, p.921.
!   Values NQAD=1:4 are identical to those listed in Clough et al 1992.
!   Digits here have been checked!
!
! REFERENCES
!     Handbook of Mathematical Functions 
!     M.Abramowitz and I.A.Stegun (Eds)
!     9th Dover printing, New York, 1972.
!
!     Clough et al, J.Geophys.Res. 97, 157610-15785 (1992).
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN)  :: N    ! Degree of quadrature  (1:5)
    REAL(R8),    INTENT(OUT) :: X(N) ! List of N x-values (abscissa)
    REAL(R8),    INTENT(OUT) :: W(N) ! List of N weights (0:1)
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  SELECT CASE ( N )
  CASE ( 1 ) 
    X(1) = 0.6666666667D0
    W(1) = 0.5D0
  CASE ( 2 ) 
    X(1) = 0.3550510257D0
    X(2) = 0.8449489743D0
    W(1) = 0.1819586183D0
    W(2) = 0.3180413817D0
  CASE ( 3 ) 
    X(1) = 0.2123405382D0
    X(2) = 0.5905331356D0
    X(3) = 0.9114120405D0
    W(1) = 0.0698269799D0
    W(2) = 0.2292411064D0
    W(3) = 0.2009319137D0
  CASE ( 4 ) 
    X(1) = 0.1397598643D0
    X(2) = 0.4164095676D0
    X(3) = 0.7231569864D0
    X(4) = 0.9428958039D0
    W(1) = 0.0311809710D0
    W(2) = 0.1298475476D0
    W(3) = 0.2034645680D0
    W(4) = 0.1355069134D0
  CASE ( 5 ) 
    X(1) = 0.0985350858D0
    X(2) = 0.3045357266D0
    X(3) = 0.5620251898D0
    X(4) = 0.8019865821D0
    X(5) = 0.9601901429D0
    W(1) = 0.0157479145D0
    W(2) = 0.0739088701D0
    W(3) = 0.1463869871D0
    W(4) = 0.1671746381D0
    W(5) = 0.0967815902D0
  CASE DEFAULT
    STOP 'F-GAUQAD: Argument N out of range'
  END SELECT
!
END SUBROUTINE GAUQAD
END MODULE GAUQAD_SUB
