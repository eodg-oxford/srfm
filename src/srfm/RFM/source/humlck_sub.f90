MODULE HUMLCK_SUB
CONTAINS
SUBROUTINE HUMLCK ( X, Y, V )
!
! VERSION
!   16MAY24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Calculate Humlicek complex prob.function
!   Called by MIXSHP or CHISHP.
!
!   The Voigt lineshape formulation:
!
!               g(X,Y) = S * g0 * K(X,Y)
!               g0 = 1/Ad * SQRT(ln2/pi)
!               X = (nu - nu0)/Ad *SQRT(ln2)
!               Y = Al/Ad *SQRT(ln2)
!               K(X,Y) = Y/pi * 
!               INT^(+infty)_(-infty){exp(-t**2)/[Y**2 + (X-t)**2]}dt
!
!   This routine calculates the complex probability function using a 
!   vectorized version of the Humlicek JQSRT V27 437 1982 paper. 
!
! REFERENCES
!   Humlicek, J.
!   Optimised computation of the Voigt and complex probability functiosn
!   J.Q.S.R.T, 27, 437-444 (1982)
!   doi:10.1016/0022-4073(82)90078-4
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4),    INTENT(IN)  :: X(:) ! X array 
    REAL(R4),    INTENT(IN)  :: Y    ! Y value 
    COMPLEX(R4), INTENT(OUT) :: V(:) ! Absorption 
!
! LOCAL VARIABLES
    INTEGER(I4) :: I
    INTEGER(I4) :: N 
    REAL(R4)    :: S1V, S2V
    COMPLEX(R4) :: U,T
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  N = SIZE ( X ) 
! Sort the (x,y) pairs into the 4 regions of the Humlicek expressions of the 
! Voigt line profile.
!
  DO I = 1, N 
    S1V = ABS(X(I)) + Y
    S2V = ( 0.195 * ABS(X(I)) ) - 0.176
    T = CMPLX ( Y, -X(I) )
!
! Region 1 of Humlicek
!
    IF ( S1V .GE. 15.0 ) THEN
      V(I) = T * 0.5641896 / ( 0.5 + (T*T) )
!
! Region 2 OF Humlicek
!
    ELSEIF ( S1V .GE. 5.5 ) THEN
      U = T * T
      V(I) = T * ( 1.410474 + U * 0.5641896 ) / ( 0.75 + U * ( 3.0 + U ) )
!
! Region 3 of Humlicek
!
    ELSE IF ( Y .GE. S2V ) THEN
      V(I) = ( 16.4955 + T * ( 20.20933 + T * ( 11.96482 + &
                T * ( 3.778987 + T * 0.5642236 ) ) ) ) / &
             ( 16.4955 + T * ( 38.82363 + T * ( 39.27121 + & 
                T * ( 21.69274 + T * ( 6.699398 + T ) ) ) ) )
!
! Region 4 of Humlicek
!
    ELSE
      U = T * T
      V(I) = CEXP ( U ) - T * ( 36183.31 - U * ( 3321.9905 - &
               U * ( 1540.787 - U * ( 219.0313 - U * ( 35.76683 - & 
               U * ( 1.320522 - U * 0.56419 ) ) ) ) ) ) /  &
               ( 32066.6 - U * ( 24322.84 - U * ( 9022.228 -  & 
               U * ( 2186.181 - U * ( 364.2191 - U * ( 61.57037 -  & 
               U * ( 1.841439-U ) ) ) ) ) ) )  
    END IF
  END DO
! 
END SUBROUTINE HUMLCK
END MODULE HUMLCK_SUB
