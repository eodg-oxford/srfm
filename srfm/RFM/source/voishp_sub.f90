MODULE VOISHP_SUB
CONTAINS
SUBROUTINE VOISHP ( DWNO, K )
!
! VERSION
!   10MAY24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Calculate Voigt Line shape
!   Called by LINSHP
!   Uses path-adjusted line data in ADJCOM
!
!   This modifies the GENLN2 VOIGT algorithm - see VOIGL2 for original.
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
!   modified version of the Humlicek algorithm (JQSRT V27 437 1982)
!   accepted for publication in JQSRT 1999.
!   The calculation is performed for the array of x,y pairs for a given line 
!   over the fine mesh points of the current wide mesh. 
!
!
! REFERENCES
!   Humlicek, J. (1982)
!   Optimised computation of the Voigt and complex probability functiosn 
!   J. Quant. Spectrosc. Radiat. Transfer, 27, 437-444. 
!   doi:10.1016/0022-4073(82)90078-4 
!
!   Wells, R. J. (1999)
!   Rapid approximation to the Voigt/Faddeeva function and its derivatives
!   J. Quant. Spectrosc. Radiat. Transfer,  62, 29-48
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE ADJCOM_DAT ! Path-adjusted line data
    USE PHYCON_DAT, ONLY: PI
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(IN)  :: DWNO(:) ! Array of Wavenumbers [/cm]
    REAL(R4), INTENT(OUT) :: K(:)    ! Absorption 
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: RSQRPI = 1.0 / SQRT ( PI ) 
    REAL(R4), PARAMETER :: SQRLN2 = SQRT ( LOG ( 2.0 ) ) !
    REAL(R4), PARAMETER :: C(0:5) = (/ 1.0117281, -0.75197147, 0.012557727, &
                      0.010022008, -0.00024206814, 0.00000050084806 /)
    REAL(R4), PARAMETER :: S(0:5) = (/ 1.393237, 0.23115241, -0.15535147, &
                      0.0062183662, 0.000091908299, -0.00000062752596 /)
    REAL(R4), PARAMETER :: T(0:5) = (/ 0.31424038, 0.94778839, 1.5976826, &
                                   2.2795071, 3.0206370, 3.8897249 /)
    REAL(R4), PARAMETER :: Y0 = 1.5                      ! for CPF12 algorithm
    REAL(R4), PARAMETER :: Y0PY0 = Y0+Y0
    REAL(R4), PARAMETER :: Y0Q = Y0 * Y0
!
! LOCAL VARIABLES
    INTEGER(I4) :: I, J                              ! Loop variables
    INTEGER(I4) :: N                                 ! Size of DWNO array
    INTEGER(I4) :: RG1, RG2, RG3                     ! y polynomial flags
    REAL(R4)    :: ABX, XQ, YQ, YRRTPI               ! |x|, x^2, y^2, y/SQRT(pi)
    REAL(R4)    :: XLIM0, XLIM1, XLIM2, XLIM3, XLIM4 ! |x| on region boundaries
    REAL(R4)    :: A0, D0, D2, E0, E2, E4, H0, H2, H4, H6 ! W4 temp. variables
    REAL(R4)    :: P0, P2, P4, P6, P8, Z0, Z2, Z4, Z6, Z8
    REAL(R4)    :: XP(0:5), XM(0:5), YP(0:5), YM(0:5) ! CPF12 temporary values
    REAL(R4)    :: MQ(0:5), PQ(0:5), MF(0:5), PF(0:5)
    REAL(R4)    :: D, YF, YPY0, YPY0Q, REPWID, XI, Y  
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  N = SIZE ( DWNO ) 
  REPWID = SQRLN2 / DOPADJ
  Y      = REPWID * WIDADJ
  YQ     = Y * Y                                ! y^2
!
  IF ( Y .GE. 70.55 ) THEN                    ! Region 0
    DO I = 1, N                               
      XI   = SNGL ( DWNO(I) - WNOADJ ) * REPWID            
      K(I) = STRADJ * REPWID * Y / ( PI * ( XI * XI + YQ ) )
    ENDDO
    RETURN
  ENDIF
!
  RG1 = 1                                     ! Set flags
  RG2 = 1
  RG3 = 1
!
  YRRTPI = Y * RSQRPI                           ! y/SQRT(pi)
  XLIM0  = SQRT ( 15100.0 + Y * (40.0 - Y * 3.6 ) )   ! y<70.55
  IF ( Y .GE. 8.425 ) THEN
    XLIM1 = 0.0
  ELSE
    XLIM1 = SQRT ( 164.0 - Y*(4.3 + Y*1.8) )
  ENDIF
  XLIM2 = 6.8 - Y
  XLIM3 = 2.4*Y
  XLIM4 = 18.1*Y + 1.65
!
  IF ( Y .LE. 0.000001 ) THEN                 ! When y<10^-6
    XLIM1 = XLIM0                              ! avoid W4 algorithm
    XLIM2 = XLIM0
  ENDIF
!
  DO I = 1, N
    XI = SNGL ( DWNO(I) - WNOADJ ) * REPWID        ! Loop over all points
    ABX = ABS ( XI )                           ! |x|
    XQ  = ABX * ABX                              ! x^2
    IF ( ABX .GE. XLIM0 ) THEN             ! Region 0 algorithm
      K(I) = YRRTPI / ( XQ + YQ )
!
    ELSE IF ( ABX .GE. XLIM1 ) THEN            ! Humlicek W4 Region 1
      IF ( RG1 .NE. 0 ) THEN                    ! First point in Region 1
        RG1 = 0
        A0 = YQ + 0.5                            ! Region 1 y-dependents
        D0 = A0 * A0
        D2 = YQ + YQ - 1.0
      ENDIF
      D = RSQRPI / ( D0 + XQ * ( D2 + XQ ) )
      K(I) = D * Y * ( A0 + XQ )
!
    ELSE IF ( ABX .GE. XLIM2 ) THEN            ! Humlicek W4 Region 2 
      IF ( RG2 .NE. 0 ) THEN                    ! First point in Region 2
        RG2 = 0
        H0 =  0.5625 + YQ*(4.5 + YQ*(10.5 + YQ*(6.0 + YQ))) ! Region 2 y-deps
        H2 = -4.5    + YQ*(9.0 + YQ*( 6.0 + YQ* 4.0))
        H4 = 10.5    - YQ*(6.0 - YQ*  6.0)
        H6 = -6.0    + YQ* 4.0
        E0 =  1.875  + YQ*(8.25 + YQ*(5.5 + YQ))
        E2 =  5.25   + YQ*(1.0  + YQ* 3.0)
        E4 =  0.75*H6
      ENDIF
      D = RSQRPI / (H0 + XQ*(H2 + XQ*(H4 + XQ*(H6 + XQ))))
      K(I) = D*Y   *(E0 + XQ*(E2 + XQ*(E4 + XQ)))
!
    ELSE IF ( ABX .LT. XLIM3 ) THEN             ! Humlicek W4 Region 3
      IF ( RG3 .NE. 0 ) THEN                     ! First point in Region 3
        RG3 = 0
        Z0 = 272.1014     + Y*(1280.829 + Y*(2802.870 + Y*(3764.966 & ! y-deps
             + Y*(3447.629 + Y*(2256.981 + Y*(1074.409 + Y*(369.1989 &
             + Y*(88.26741 + Y*(13.39880 + Y)))))))))
         Z2 = 211.678      + Y*(902.3066 + Y*(1758.336 + Y*(2037.310 &
                           + Y*(1549.675 + Y*(793.4273 + Y*(266.2987 &
                           + Y*(53.59518 + Y*5.0)))))))
         Z4 = 78.86585     + Y*(308.1852 + Y*(497.3014 + Y*(479.2576 &
                           + Y*(269.2916 + Y*(80.39278 + Y*10.0))))) 
         Z6 = 22.03523     + Y*(55.02933 + Y*(92.75679 + Y*(53.59518 &
                           + Y*10.0)))
         Z8 = 1.496460     + Y*(13.39880 + Y*5.0)
         P0 = 153.5168     + Y*(549.3954 + Y*(919.4955 + Y*(946.8970 &
             + Y*(662.8097 + Y*(328.2151 + Y*(115.3772 + Y*(27.93941 &
                           + Y*(4.264678 + Y*0.3183291))))))))
         P2 = -34.16955    + Y*(-1.322256+ Y*(124.5975 + Y*(189.7730 &
                           + Y*(139.4665 + Y*(56.81652 + Y*(12.79458 &
                           + Y*1.2733163))))))
         P4 = 2.584042     + Y*(10.46332 + Y*(24.01655 + Y*(29.81482 &
                           + Y*(12.79568 + Y*1.9099744))))
         P6 = -0.07272979  + Y*(0.9377051+ Y*(4.266322 + Y*1.273316))
         P8 = 0.0005480304 + Y*0.3183291
       ENDIF
       D = 1.7724538 / (Z0 + XQ*(Z2 + XQ*(Z4 + XQ*(Z6 + XQ*(Z8+XQ)))))
       K(I) = D*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))))
     ELSE                                         ! Humlicek CPF12 algorithm
       YPY0 = Y + Y0
       YPY0Q = YPY0*YPY0
       K(I) = 0.0
       DO J = 0, 5
         D = XI - T(J)
         MQ(J) = D*D
         MF(J) = 1.0 / (MQ(J) + YPY0Q)
         XM(J) = MF(J)*D
         YM(J) = MF(J)*YPY0
         D = XI + T(J)
         PQ(J) = D*D
         PF(J) = 1.0 / (PQ(J) + YPY0Q)
         XP(J) = PF(J)*D
         YP(J) = PF(J)*YPY0
       ENDDO
!
       IF ( ABX .LE. XLIM4 ) THEN                   ! Humlicek CPF12 Region I
         DO J = 0, 5
            K(I) = K(I) + C(J)*(YM(J)+YP(J)) - S(J)*(XM(J)-XP(J))
         ENDDO

       ELSE                                        ! Humlicek CPF12 Region II
         YF   = Y + Y0PY0
         DO J = 0, 5
          K(I) = K(I) &
           + (C(J)*(MQ(J)*MF(J)-Y0*YM(J)) + S(J)*YF*XM(J)) / (MQ(J)+Y0Q) &
           + (C(J)*(PQ(J)*PF(J)-Y0*YP(J)) - S(J)*YF*XP(J)) / (PQ(J)+Y0Q)
         ENDDO
         K(I) = Y * K(I) + EXP ( -XQ )
       ENDIF
     ENDIF
     K(I) = STRADJ * RSQRPI * REPWID * K(I)
   ENDDO
!
END SUBROUTINE VOISHP
END MODULE VOISHP_SUB
