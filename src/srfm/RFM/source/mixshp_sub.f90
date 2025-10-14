MODULE MIXSHP_SUB
CONTAINS
SUBROUTINE MIXSHP ( DWNO, ABSORP )
!
! VERSION
!   08MAY24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Calculate Voigt Line shape allowing for line mixing
!   Called by LINSHP if MIX Flag enabled.
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
!   The calculation is performed for the array of x,y pairs for a given line 
!   over the fine mesh points of the current wide mesh. 
!
!   Uses path-adjusted line data in ADJCOM
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ADJCOM_DAT ! Path-adjusted line data
    USE PHYCON_DAT, ONLY: PI
!
! SUBROUTINES
    USE HUMLCK_SUB ! Calculate Humlicek complex prob.function
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(IN)  :: DWNO(:)   ! Array of Wavenumbers [cm-1]
    REAL(R4), INTENT(OUT) :: ABSORP(:) ! Absorption 
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: SQRLN2 = SQRT ( LOG ( 2.0 ) ) ! 0.83255461
    REAL(R4), PARAMETER :: RSQRPI = 1.0 / SQRT ( PI )    ! 0.5641895
!
! LOCAL VARIABLES
    REAL(R4)    :: H0
    REAL(R4)    :: REPWID   
    REAL(R4)    :: X(SIZE(DWNO))
    REAL(R4)    :: Y
    COMPLEX(R4) :: V(SIZE(DWNO)) 
!
! EXECUTABLE CODE -------------------------------------------------------------
!
   REPWID = SQRLN2 / DOPADJ
   H0 = REPWID * RSQRPI * STRADJ
   Y = WIDADJ * REPWID
   X = SNGL ( DWNO - WNOADJ ) * REPWID
   CALL HUMLCK ( X, Y, V )
!
! Compute absorption due to Voigt line shape
   ABSORP = H0 * ( REAL ( V ) + YMXADJ * AIMAG ( V ) )
!
END SUBROUTINE MIXSHP
END MODULE MIXSHP_SUB
