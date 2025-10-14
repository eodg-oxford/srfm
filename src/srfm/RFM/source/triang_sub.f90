MODULE TRIANG_SUB
CONTAINS
SUBROUTINE TRIANG ( N, X, Y, NTRI, IDXTRI ) 
!
! VERSION
!   19FEB20 AD Bug#24 Amend perimeter check. Checked.
!   10FEB19 AD Bug#16 Add NTRI argument
!   01MAY17 AD F90 conversion of part of triang.for. Checked.
!
! DESCRIPTION
!   2-D interpolation of irregular grid using triangulation
!   General purpose module.
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
    INTEGER(I4), INTENT(IN)  :: N           ! No. pairs of (x,y) coordinates
    REAL(R4),    INTENT(IN)  :: X(:)        ! x-coordinates of points in field
    REAL(R4),    INTENT(IN)  :: Y(:)        ! y-coordinates of points in field
    INTEGER(I4), INTENT(OUT) :: NTRI        ! Number of triangles constructed
    INTEGER(I4), ALLOCATABLE, &
                 INTENT(OUT) :: IDXTRI(:,:) ! Indices of triangle vertices
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: SMALL = 1.0E-5 !  Min sig diff in cosine. 
!
! LOCAL VARIABLES
    LOGICAL     :: FINITE       ! T=space has finite area, F=straight line
    LOGICAL     :: INSIDE       ! T=inside triangle, F=outside
    LOGICAL     :: PERIM(N)     ! T=perimeter point, F=internal point
    INTEGER(I4) :: I,J,K,L      ! Indices of original points X,Y
    INTEGER(I4) :: II,JJ,KK,LL  ! Saved values of indices
    INTEGER(I4) :: IDXPER(N)    ! Index of perimeter points
    INTEGER(I4) :: ISID         ! Side counter (of triangle)
    INTEGER(I4) :: ITRI,JTRI,KTRI,LTRI ! Triangle counters
    INTEGER(I4) :: IXMAX        ! (Original) Index of max value of x-coordinate
    INTEGER(I4) :: NPER         ! Number of perimeter points
    REAL(R4)    :: CIJK,CIJL,CJKL,CKIL ! Cross products of vectors
    REAL(R4)    :: CKLI,CLIJ    ! Cross products of vectors 
    REAL(R4)    :: COSIJK, COSKLI ! Cosine of angles IJK and KLI
    REAL(R4)    :: COSINE       ! Cosine of angle between two vectors
    REAL(R4)    :: COSMAX       ! Maximum value of cosine
    REAL(R4)    :: DIJK, DKLI   ! Dot products of vectors
    REAL(R4)    :: DLIJ, DLJK, DLKI  ! Dot products of vectors
    REAL(R4)    :: FLIP         ! Benefit (+ve reduced dist) from diagonal flip
    REAL(R4)    :: FLIPMX       ! Maximum value of FLIP found so far
    REAL(R4)    :: R            ! Distance between points
    REAL(R4)    :: RMIN         ! Minimum distance from interp.point 
    REAL(R4)    :: SIJ, SJK, SKL, SLI ! Distances between points
    REAL(R4)    :: SINE         ! Sine of angle
    REAL(R4)    :: UIJ,UIL,UKI  ! x-coordinate of unit vectors along perimeter
    REAL(R4)    :: VIJ,VIL,VKI  ! y-coordinate of unit vectors along perimeter
    REAL(R4)    :: XMAX         ! Maximum value of x-coordinate
    INTEGER(I4), ALLOCATABLE :: IDXSAV(:,:) ! Saved IDXTRI during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Deal with special cases of N=1 and N=2
  IF ( N .EQ. 1 ) THEN  ! replicate (x,y) coordinate 3 times
    ALLOCATE ( IDXTRI(1,3) )
    IDXTRI(1,:) = 1
    NTRI = 1
    RETURN
  ELSE IF ( N .EQ. 2 ) THEN ! replicate 1st (x,y) coordinate twice
    ALLOCATE ( IDXTRI(1,3) ) 
    IDXTRI(1,1:2) = 1
    IDXTRI(1,3) = 2
    NTRI = 1
    RETURN
  END IF
!
! Step#1: establish point with the maximum x-coordinate 
  XMAX = MAXVAL ( X )
  IXMAX = MAXLOC ( X, 1 ) 
!
! Step#2: Determine the array of points that form the perimeter
! Given a known point along the perimeter I ...
! for each unit vector KI along the perimeter (leading to point (I)), we want 
! to find the next point (J) whose unit vector IJ lies at the smallest angle 
! relative to KI (because we are proceeding anti-c/w around the perimeter, all
! angles are between 0 and 180 degrees by definition). 
! Cos(angle) between vectors is a suitable parameter to maximise since this
! decreases monotonically from 1 (0 degrees) to -1 (180 degrees).
!
  PERIM = .FALSE.
  NPER = 0
  PERIM(IXMAX) = .TRUE.      ! Prevent IXMAX selected first time
  I = IXMAX       ! Start at a known perimeter point, (X(IXMAX),Y(IXMAX))
  UKI = 0.0       ! Start with unit vector KI in +y direction
  VKI = 1.0  
  J = 0
  DO WHILE ( J .NE. IXMAX )        ! Stop when returns to point IXMAX
    RMIN = 0.0
    COSMAX = -2.0                  ! Initialise search for next perim.point
    DO L = 1, N                    ! Try all points in field (L)
      IF ( .NOT. PERIM(L) ) THEN    ! avoiding previous perimeter points
        R  = SQRT ( (X(L)-X(I))**2 + (Y(L)-Y(I))**2 )
        UIL = ( X(L) - X(I) ) / R       ! Unit vector in dir. (I) to (J)
        VIL = ( Y(L) - Y(I) ) / R
        COSINE = UKI * UIL + VKI * VIL  ! Dot-product of KI.IL = cos(angle)
        IF ( COSINE .GT. COSMAX+SMALL ) THEN  ! New direction
          RMIN   = R
          J      = L
          UIJ    = UIL
          VIJ    = VIL
          COSMAX = COSINE
        ELSE IF ( COSINE .GT. COSMAX-SMALL .AND. R .LT. RMIN ) THEN
          J = L                                   ! Same direction, but closer
          RMIN = R
        END IF
      END IF
    END DO   
    NPER = NPER + 1                ! Add point I to perimeter list
    IDXPER(NPER) = J
    PERIM(J) = .TRUE.
    I   = J                        ! Set current perimeter point I
    UKI = UIJ                      ! Unit vector along last perim.section
    VKI = VIJ
! Allow point IXMAX to be selected to complete perimeter but only after two 
! other perimeter points have been identified (Bug#24)
     IF ( NPER .EQ. 2 ) PERIM(IXMAX) = .FALSE.  
  END DO
  PERIM(IXMAX) = .TRUE.
!
! Step#3: Form triangles by drawing lines from IXMAX to every other perimeter
! point. Since there are NPER perimeter points, there must be NPER-2 triangles.
! Also check if perimeter section subtends a finite angle at IXMAX: if it 
! does then the grid has a finite area (otherwise it is a straight line)
  NTRI = NPER - 2
  ALLOCATE ( IDXTRI(NTRI,3) )
  FINITE = .FALSE.
  DO ITRI = 1, NTRI
    I = IXMAX
    J = IDXPER(ITRI)
    K = IDXPER(ITRI+1)
    IDXTRI(ITRI,1) = I
    IDXTRI(ITRI,2) = J
    IDXTRI(ITRI,3) = K
    IF ( .NOT. FINITE ) THEN
      COSINE = ((X(J)-X(I))*(X(K)-X(I)) + (Y(J)-Y(I))*(Y(K)-Y(I))) &
               / SQRT ( (X(J)-X(I))**2 + (Y(J)-Y(I))**2 )      &
               / SQRT ( (X(K)-X(I))**2 + (Y(K)-Y(I))**2 )
      FINITE = ( 1.0 - ABS ( COSINE ) ) .GT. SMALL
    END IF
  END DO
!
! Step#4: For each non-perimeter point L, find which existing triangle IJK it 
! lies inside, and sub-divide that triangle into three smaller triangles by 
! drawing lines from point to each vertex: IJL, JKL, KIL
  DO L = 1, N
    IF ( .NOT. PERIM(L) ) THEN
      DO ITRI = 1, NTRI
        I = IDXTRI(ITRI,1)
        J = IDXTRI(ITRI,2)
        K = IDXTRI(ITRI,3)
! If the point L is enclosed, then the cross-products of the vectors formed
! by each of the sides in turn (IJ, JK, KI) and the lines from each vertex to L
! (JL, KL, IL) must all be of the same sign. 
        CIJL = (X(J)-X(I))*(Y(L)-Y(J)) - (Y(J)-Y(I))*(X(L)-X(J))
        CJKL = (X(K)-X(J))*(Y(L)-Y(K)) - (Y(K)-Y(J))*(X(L)-X(K))
        CKIL = (X(I)-X(K))*(Y(L)-Y(I)) - (Y(I)-Y(K))*(X(L)-X(I))
        INSIDE = (CIJL.GE.0.0 .AND. CJKL.GE.0.0 .AND. CKIL.GE.0.0) .OR. &
                 (CIJL.LE.0.0 .AND. CJKL.LE.0.0 .AND. CKIL.LE.0.0) 
! However, if all three cross-products are zero, point lies on same straight
! line as the other three points, so need to determine if interpolation 
! ("inside") or extrapolation ("outside"). Do this by looking for at least one 
! negative dot product between LI.LJ, LJ.LK and LK.LI
        IF ( CIJL.EQ.0.0 .AND. CJKL.EQ.0.0 .AND. CKIL.EQ.0.0 ) THEN
          DLIJ = (X(I)-X(L))*(X(J)-X(L)) + (Y(I)-Y(L))*(Y(J)-Y(L))
          DLJK = (X(J)-X(L))*(X(K)-X(L)) + (Y(J)-Y(L))*(Y(K)-Y(L))
          DLKI = (X(K)-X(L))*(X(I)-X(L)) + (Y(K)-Y(L))*(Y(I)-Y(L))
          INSIDE = DLIJ.LT.0.0 .OR. DLJK.LT.0.0 .OR. DLKI.LT.0.0
        END IF
        IF ( INSIDE ) THEN
          CALL MOVE_ALLOC ( IDXTRI, IDXSAV ) 
          ALLOCATE ( IDXTRI(NTRI+2,3) )
          IDXTRI(1:NTRI,:) = IDXSAV
          IDXTRI(ITRI,3) = L        ! ITRI becomes IJL
          IDXTRI(NTRI+1,1) = I      ! NTRI+1 becomes IKL
          IDXTRI(NTRI+1,2) = K
          IDXTRI(NTRI+1,3) = L
          IDXTRI(NTRI+2,1) = J      ! NTRI+1 becomes JKL
          IDXTRI(NTRI+2,2) = K
          IDXTRI(NTRI+2,3) = L
          NTRI = NTRI + 2           ! 2 more triangles added (net effect)
          GOTO 100
        END IF
      END DO
      STOP 'F-TRIANG: logical error#1'  ! All points must lie in a triangle
 100  CONTINUE
    END IF
  END DO
!
! At this point the entire 2D field has been triangulated 
!
! Step#5: Optimising the triangulation (try to make triangles more equilateral)
! Look for pairs of triangles IJK, KLI which share a side IK
  DO
    FLIPMX = 0.0
    DO JTRI = 1, NTRI-1                       ! Loop over all triangles
      DO ISID = 1, 3                          ! Loop over each side of JTRI
        I = IDXTRI(JTRI,ISID)                 ! I = 1, 2, 3
        J = IDXTRI(JTRI,1+MOD(ISID,3))        ! J = 2, 3, 1
        K = IDXTRI(JTRI,1+MOD(ISID+1,3))      ! K = 3, 1, 2
        DO LTRI = JTRI+1, NTRI                ! Compare with other triangles
          IF ( ( I .EQ. IDXTRI(LTRI,1) .OR.  &
                 I .EQ. IDXTRI(LTRI,2) .OR.  &
                 I .EQ. IDXTRI(LTRI,3)      ) .AND.  &
               ( K .EQ. IDXTRI(LTRI,1) .OR.  &
                 K .EQ. IDXTRI(LTRI,2) .OR.  &
                 K .EQ. IDXTRI(LTRI,3)      ) ) THEN  ! Found shared boundary
            L = IDXTRI(LTRI,1)                        ! Establish 4th point, L
            IF ( L .EQ. I .OR. L .EQ. K ) THEN
              L = IDXTRI(LTRI,2)
              IF ( L .EQ. I .OR. L .EQ. K ) L = IDXTRI(LTRI,3)
            END IF
! We now have a quadrilateral IJKL. To allow the diagonal IK to be flipped to 
! JL it is necessary that each corner is convex (inside angle up to 180deg)
! For this, the sine (ie cross-product) between vectors representing 
! successive sides IJ, JK, JL, LI must all be the same sign.
            CIJK = (X(J)-X(I))*(Y(K)-Y(J))-(X(K)-X(J))*(Y(J)-Y(I))
            CJKL = (X(K)-X(J))*(Y(L)-Y(K))-(X(L)-X(K))*(Y(K)-Y(J))
            CKLI = (X(L)-X(K))*(Y(I)-Y(L))-(X(I)-X(L))*(Y(L)-Y(K))
            CLIJ = (X(I)-X(L))*(Y(J)-Y(I))-(X(J)-X(I))*(Y(I)-Y(L))
            IF ( ( CIJK .GE. 0.0 .AND. CJKL .GE. 0.0 .AND. &
                   CKLI .GE. 0.0 .AND. CLIJ .GE. 0.0       ) .OR. &
                 ( CIJK .LE. 0.0 .AND. CJKL .LE. 0.0 .AND. &
                   CKLI .LE. 0.0 .AND. CLIJ .LE. 0.0       ) ) THEN
!
! Flip if the circle through one triangle includes the fourth point.
! This is true if the sum of angles IJK+KLI > 180. If this is true, then it
! will not be true for the other pair of corners (since sum all angles=360),
! so "benefit" is increase in sin(sum opposite corners).
              DIJK = (X(I)-X(J))*(X(K)-X(J))+(Y(I)-Y(J))*(Y(K)-Y(J))
              DKLI = (X(K)-X(L))*(X(I)-X(L))+(Y(K)-Y(L))*(Y(I)-Y(L))
              SIJ = SQRT ( (X(J)-X(I))**2 + (Y(J)-Y(I))**2 )
              SJK = SQRT ( (X(K)-X(J))**2 + (Y(K)-Y(J))**2 )
              SKL = SQRT ( (X(L)-X(K))**2 + (Y(L)-Y(K))**2 )             
              SLI = SQRT ( (X(I)-X(L))**2 + (Y(I)-Y(L))**2 )
              COSIJK = DIJK / ( SIJ * SJK )
              COSKLI = DKLI / ( SKL * SLI )
              SINE = COSIJK * SQRT ( MAX ( 0.0, 1.0 - COSKLI**2 ) ) + &
                     COSKLI * SQRT ( MAX ( 0.0, 1.0 - COSIJK**2 ) )
              FLIP = - SINE
!
! If two "flat" triangles (so cos=-1), swap diagonal to shorten boundary
              IF (COSIJK+1.0 .LE. SMALL .AND. COSKLI+1.0 .LE. SMALL) FLIP = 2.0
!
              IF ( FLIP .GT. FLIPMX ) THEN      ! Found best so far
                FLIPMX = FLIP
                II = I
                JJ = J
                KK = K
                LL = L
                KTRI = JTRI
                ITRI = LTRI
              END IF     ! end test for best flip so far
            END IF       ! end test for flippable diagonals
          END IF         ! end test for shared boundary
        END DO           ! end loop over LTRI
      END DO             ! end loop over sides of JTRI
    END DO               ! end loop over JTRI
!
! Flip diagonal for pair of triangles which gives most improvement
! So triangles IJK, KLI become LIJ, JKL
    IF ( FLIPMX .GT. SMALL ) THEN
      IDXTRI(ITRI,1) = LL
      IDXTRI(ITRI,2) = II
      IDXTRI(ITRI,3) = JJ
      IDXTRI(KTRI,1) = JJ
      IDXTRI(KTRI,2) = KK
      IDXTRI(KTRI,3) = LL
    ELSE
      EXIT                 ! Repeat until no more improvements found
    END IF
  END DO
!
END SUBROUTINE TRIANG
END MODULE TRIANG_SUB
