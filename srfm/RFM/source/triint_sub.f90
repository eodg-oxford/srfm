MODULE TRIINT_SUB
CONTAINS
SUBROUTINE TRIINT ( XO, YO, N, X, Y, NTRI, IDXTRI, IFIT, WFIT, XFIT, YFIT ) 
!
! VERSION
!   05APR24 AD Checked.
!   08FEB19 AD Bug#16 Add NTRI argument
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   2-D interpolation of irregular grid using triangulation
!   Called by SPCXSC.
!
!   A plane is fitted to the three triangulation points and the interpolated
!   value on this plane is expressed as a set of weights applied to the 
!   value of the function at the triangulation points.
!   If the interpolation point lies outside the triangulated field, a value
!   is returned from the closest point along the boundary. In this case 
!   the arguments XFIT, YFIT will be different from the input coordinates 
!   XO,YO.
!
!   NB: the output indices of the interpolated points may be repeated 
!   (eg if only one point is chosen, the index will be repeated three times
!   and the corresponding weights will be 1.0, 0.0 and 0.0).
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4),    INTENT(IN)  :: XO          ! x-coordinate of interpolated point
    REAL(R4),    INTENT(IN)  :: YO          ! y-coordinate of interpolated point
    INTEGER(I4), INTENT(IN)  :: N           ! No.(x,y) coords in domain
    REAL(R4),    INTENT(IN)  :: X(N)        ! List of x-coords of domain points
    REAL(R4),    INTENT(IN)  :: Y(N)        ! List of y-coords of domain points
    INTEGER(I4), INTENT(IN)  :: NTRI        ! Number of triangles
    INTEGER(I4), INTENT(IN)  :: IDXTRI(:,:) ! Indices of triangle vertices
    INTEGER(I4), INTENT(OUT) :: IFIT(3)     ! Indices of points used for interp.
    REAL(R4),    INTENT(OUT) :: WFIT(3)     ! Weights of points used for interp.
    REAL(R4),    INTENT(OUT) :: XFIT        ! x-coord of point actually fitted 
    REAL(R4),    INTENT(OUT) :: YFIT        ! y-coord of point actually fitted
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: SMALL = 1.0E-5  !  Min. significant diff in cosine.
!
! LOCAL VARIABLES
    INTEGER(I4) :: I,J,K ! Indices of original points X,Y
    INTEGER(I4) :: ISID  ! Side counter (of triangle)
    INTEGER(I4) :: ITRI  ! Triangle counters
    REAL(R4)    :: CIJO,CJKO,CKIO         ! Cross products of vectors
    REAL(R4)    :: CIJTOL, CJKTOL, CKITOL ! Cross-product tolerances
    REAL(R4)    :: DIJO, DJIO             ! Dot products of vectors
    REAL(R4)    :: FACTX,FACTY,FACTO      ! Factors in comp. of triang. wgts
    REAL(R4)    :: R     ! Distance between points
    REAL(R4)    :: RMIN  ! Minimum distance from interp.point 
    REAL(R4)    :: RTOL  ! Tolerance in comparing R
    REAL(R4)    :: S     ! Length of vector
    REAL(R4)    :: SMIN  ! Minimum length of vector
    REAL(R4)    :: SJP2  ! Square of length of vector JP 
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Deal with special cases of N=1 and N=2 (ignoring IDXTRI)
  IF ( N .EQ. 1 ) THEN  ! Simply take value of single grid point
    IFIT(1) = 1
    IFIT(2) = 1
    IFIT(3) = 1
    WFIT(1) = 1.0
    WFIT(2) = 0.0
    WFIT(3) = 0.0
    XFIT    = X(1)
    YFIT    = Y(1)
    RETURN
  ELSE IF ( N .EQ. 2 ) THEN ! Interpolate normally from line IJ
    I = 1
    J = 2
    IFIT(3) = 1
    WFIT(3) = 0.0
    DIJO = (X(I)-X(J)) * (XO-X(J)) + (Y(I)-Y(J)) * (YO-Y(J))
    DJIO = (X(J)-X(I)) * (XO-X(I)) + (Y(J)-Y(I)) * (YO-Y(I))
    IF ( DIJO .GT. 0.0 .AND. DJIO .GT. 0.0 ) THEN  ! Interpolate
      SJP2 = DIJO**2 / ( (X(J)-X(I))**2 + (Y(J)-Y(I))**2 )
      IFIT(1) = I
      IFIT(2) = J
      WFIT(1) = SJP2 / DIJO
      WFIT(2) = 1.0 - WFIT(1)
      XFIT = X(I)*WFIT(1) + X(J)*WFIT(2)
      YFIT = Y(I)*WFIT(1) + Y(J)*WFIT(2)
    ELSE                                 ! Duplicate closest point
      IFIT(2) = 1
      WFIT(2) = 0.0
      WFIT(1) = 1.0
      IF ( (X(I)-XO)**2 + (Y(I)-YO)**2 .LT. &
           (X(J)-XO)**2 + (Y(J)-YO)**2      ) THEN  ! I is closest point
        IFIT(1) = I
        XFIT    = X(I)
        YFIT    = Y(I)
      ELSE                                          ! J is closest point
        IFIT(1) = J
        XFIT    = X(J)
        YFIT    = Y(J)
      END IF
    END IF
    RETURN
  END IF            
!
! Continue from here for N > 2, ie if triangulation possible
!  NTRI = N - 2    Bug#16 - note that NTRI can be greater than N-2 
! 
! First, find closest grid point to interpolation point and set to return value
! at this point if no other location is found.
  RMIN = MINVAL ( ( XO - X )**2 + ( YO - Y )**2 )
  I    = MINLOC ( ( XO - X )**2 + ( YO - Y )**2, 1 )
  IFIT(1) = I
  IFIT(2) = 1
  IFIT(3) = 1
  WFIT(1) = 1.0
  WFIT(2) = 0.0
  WFIT(3) = 0.0
  XFIT = X(I)
  YFIT = Y(I)
  IF ( RMIN .EQ. 0.0 ) RETURN
  SMIN = 0.0
!
! Look for existing triangle enclosing interpolated point.
! Enclosing defined by cross-products, as before, but normalise by square of
! length of side to avoid choosing "flat" triangles
  DO ITRI = 1, NTRI
!    IF ( FINITE ) THEN
    I = IDXTRI(ITRI,1)
    J = IDXTRI(ITRI,2)
    K = IDXTRI(ITRI,3)
    CIJO = ( X(J)-X(I))*(YO-Y(J)) - (Y(J)-Y(I))*(XO-X(J) ) 
    CIJTOL = SMALL * SQRT (( (X(J)-X(I))**2 + (Y(J)-Y(I))**2 )* &
                           ( ( XO -X(J))**2 + ( YO -Y(J))**2 ))
    CJKO = ( X(K)-X(J))*(YO-Y(K)) - (Y(K)-Y(J))*(XO-X(K) ) 
    CJKTOL = SMALL * SQRT (( (X(K)-X(J))**2 + (Y(K)-Y(J))**2 )* &
                           ( ( XO -X(K))**2 + ( YO -Y(K))**2 ))
    CKIO = ( X(I)-X(K))*(YO-Y(I)) - (Y(I)-Y(K))*(XO-X(I) ) 
    CKITOL = SMALL * SQRT (( (X(I)-X(K))**2 + (Y(I)-Y(K))**2 )* &
                           ( ( XO -X(I))**2 + ( YO -Y(I))**2 ))
    IF ( ( CIJO .GT. CIJTOL .AND. &
           CJKO .GT. CJKTOL .AND. &
           CKIO .GT. CKITOL          ) .OR. &
         ( CIJO .LT. -CIJTOL .AND. &
           CJKO .LT. -CJKTOL .AND. &
           CKIO .LT. -CKITOL         ) ) THEN
!
! Fit plane to points I,J,K and interpolate plane to value at (XO,YO).
! This next bit of code effectively does both by inverting the 3x3 matrix to
! establish the weights (WFIT) for I,J,K in the final interpolation directly
      IFIT(1) = I
      IFIT(2) = J
      IFIT(3) = K
      FACTX = 1.0 / &
        ( ( X(I)-X(J) ) * ( Y(J)-Y(K) ) - ( X(J)-X(K) ) * ( Y(I)-Y(J) ) )
      FACTY = 1.0 / &
        ( ( Y(I)-Y(J) ) * ( X(J)-X(K) ) - ( Y(J)-Y(K) ) * ( X(I)-X(J) ) )
      FACTO = 1.0 / &
        ( X(K)*Y(J)-X(J)*Y(K) + X(I)*Y(K)-X(K)*Y(I) + X(J)*Y(I)-X(I)*Y(J) )
      WFIT(1) = XO * ( Y(J)-Y(K) ) * FACTX + &
                YO * ( X(J)-X(K) ) * FACTY + &
                     ( Y(J)*X(K) - X(J)*Y(K) ) * FACTO
      WFIT(2) = XO * ( Y(K)-Y(I) ) * FACTX + &
                YO * ( X(K)-X(I) ) * FACTY + &
                     ( Y(K)*X(I) - X(K)*Y(I) ) * FACTO
      WFIT(3) = XO * ( Y(I)-Y(J) ) * FACTX + &
                YO * ( X(I)-X(J) ) * FACTY + &
                     ( Y(I)*X(J) - X(I)*Y(J) ) * FACTO
      XFIT = XO
      YFIT = YO
      RETURN             ! Exit with triangulation completed
!      END IF            if finite
    END IF
! For each side of triangle in turn, find if the interpolation point lies 
! perpendicular ie if the triangle OIJ has acute angles at corners I and J
! so that both dot products are > 0. Note that the previous point-by-point
! test will have already found the case of dot product=0.
    DO ISID = 1, 3
      I = IDXTRI(ITRI,ISID)                ! I = 1, 2, 3
      J = IDXTRI(ITRI,1+MOD(ISID,3))       ! J = 2, 3, 1
      DIJO = (X(I)-X(J)) * (XO-X(J)) + (Y(I)-Y(J)) * (YO-Y(J))
      DJIO = (X(J)-X(I)) * (XO-X(I)) + (Y(J)-Y(I)) * (YO-Y(I))
      IF (DIJO .GT. 0.0 .AND. DJIO .GT. 0.0 )  THEN
! Project interpolation point O to point P along line IJ. 
! IF OP < RMIN, save this point (P) as potentially closest point.
! If OP = RMIN, only save this point if length of IJ (=S) is shorter than
! previous (representing a linear interpolation between two closer points).
! The distance JP = JO.cos(IJO) = DIJO/IJ since DIJO = IJ.JO.cos(IJO)
! So from right-angle triang JPO, distance OP is given by OP^2 = JO^2 - JP^2
! But DIJO = IJ.JO.cos(IJO) so cos(IJO) = DIJO/(IJ.JO)
! Also: DIJO + DJIO = IJ.(JO.cos(IJO)+IO.cos(JIO)) = IJ^2. 
        SJP2 = DIJO**2 / ( (X(J)-X(I))**2 + (Y(J)-Y(I))**2 )  
        R = (X(J)-XO)**2 + (Y(J)-YO)**2 - SJP2
        S = DIJO + DJIO
        RTOL = SJP2*SMALL + S*SMALL
        IF ( R .LT. RMIN-RTOL .OR. ( R .LT. RMIN+RTOL .AND. S.LT.SMIN) ) THEN
!
! Weight of point I is then JP/IJ. But JP = DIJO/IJ, so Wgt=JP^2/DIJO
! But DIJO = IJ.JO.cos(IJO), so PJ/IJ=DIJO/(IJ)**2 
! NB: DIJO=0 shouldn't arise since then R=RMIN from previous loop looking for
! closest point (J) and SMIN was set to 0.0
          IFIT(1) = I
          IFIT(2) = J
          IFIT(3) = 1
          WFIT(1) = SJP2 / DIJO 
          WFIT(2) = 1.0 - WFIT(1)
          WFIT(3) = 0.0
          XFIT = X(I)*WFIT(1) + X(J)*WFIT(2)
          YFIT = Y(I)*WFIT(1) + Y(J)*WFIT(2)
          RMIN = R
          SMIN = DIJO + DIJO
        END IF
      END IF
    END DO
  END DO
!
END SUBROUTINE TRIINT
END MODULE TRIINT_SUB

