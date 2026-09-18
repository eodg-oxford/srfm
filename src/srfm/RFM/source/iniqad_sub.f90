MODULE INIQAD_SUB
CONTAINS
SUBROUTINE INIQAD ( NORD ) 
!
! VERSION
!   17JUL26 AD Checked.
!   01AUG25 AD Incorporate GAUQAD. Add NORD argument
!   26APR24 AD Checked.
!   05MAR19 AD Original. Checked.
!
! DESCRIPTION
!   Initialise Gaussian quadrature for flux calculations
!   Called once by CHKSFC or TANFLX.
!   Values and Weights for Gaussian First Moment Quadrature
!     Integrating x.f(x).dx in interval x=0:1
!     Numbers taken from Abramowitz and Stegun, p.921.
!     Values NQAD=1:4 are identical to those listed in Clough et al 1992.
!     Digits here have been checked!
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
! GLOBAL DATA
    USE QADCOM_DAT ! Gaussian quadrature data
    USE PHYCON_DAT, ONLY: PI
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), OPTIONAL, INTENT(IN) :: NORD ! No.points/weights for quadrature
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: NDEF = 4  ! Default value of NORD
    REAL(R8), PARAMETER :: X1(1) = (/ 0.6666666667D0 /)
    REAL(R8), PARAMETER :: X2(2) = (/ 0.3550510257D0, 0.8449489743D0 /)
    REAL(R8), PARAMETER :: X3(3) = (/ 0.2123405382D0, 0.5905331356D0, &
                                      0.9114120405D0 /)
    REAL(R8), PARAMETER :: X4(4) = (/ 0.1397598643D0, 0.4164095676D0, &
                                      0.7231569864D0, 0.9428958039D0 /)
    REAL(R8), PARAMETER :: X5(5) = (/ 0.0985350858D0, 0.3045357266D0, & 
                                      0.5620251898D0, 0.8019865821D0, &
                                      0.9601901429D0 /)
    REAL(R8), PARAMETER :: W1(1) = (/ 0.5D0 /)
    REAL(R8), PARAMETER :: W2(2) = (/ 0.1819586183D0, 0.3180413817D0 /)
    REAL(R8), PARAMETER :: W3(3) = (/ 0.0698269799D0, 0.2292411064D0, &
                                      0.2009319137D0 /)
    REAL(R8), PARAMETER :: W4(4) = (/ 0.0311809710D0, 0.1298475476D0, &
                                      0.2034645680D0, 0.1355069134D0 /)
    REAL(R8), PARAMETER :: W5(5) = (/ 0.0157479145D0, 0.0739088701D0, &
                                      0.1463869871D0, 0.1671746381D0, &
                                      0.0967815902D0 /)
!
! LOCAL VARIABLES
    INTEGER(I4) :: N  ! NORD or NDEF
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( PRESENT ( NORD ) ) THEN
    N = NORD
  ELSE
    N = NDEF
  END IF
!
  IF ( N .EQ. 0 ) THEN ! suppress integration over solid angle
    NQAD = 1            
    RPIQAD = 1.0D0
    ALLOCATE ( XQAD(1), WQAD(1) ) 
    XQAD(1) = 1.0D0
    WQAD(1) = 1.0D0
  ELSE 
    NQAD = N
    RPIQAD = 1.0D0 / PI
    ALLOCATE ( XQAD(NQAD), WQAD(NQAD) ) 
    SELECT CASE ( NQAD ) 
      CASE ( 1 ) ; XQAD = X1 ; WQAD = W1
      CASE ( 2 ) ; XQAD = X2 ; WQAD = W2
      CASE ( 3 ) ; XQAD = X3 ; WQAD = W3
      CASE ( 4 ) ; XQAD = X4 ; WQAD = W4
      CASE ( 5 ) ; XQAD = X5 ; WQAD = W5
      CASE DEFAULT
        STOP 'F-INIQAD: Argument NORD out of range 0:5'
    END SELECT
    WQAD = WQAD * 2 * PI
  END IF
!
END SUBROUTINE INIQAD
END MODULE INIQAD_SUB
