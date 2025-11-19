MODULE RAYGRA_SUB
CONTAINS
SUBROUTINE RAYGRA ( HGTINI, PSIINI, ZENINI, IFIN, S,  HGTFIN, PSIFIN, ZENFIN )
!
! VERSION
!   03JUN24 AD Checked.
!   04APR19 AD Bug#20: Ensure any surface grazing ray avoids surface. Checked.
!   01MAY17 AD F90 conversion. Checked
!
! DESCRIPTION
!   Ray-tracing in a 2-D atmosphere (z,psi)
!   Called by GRASUM, GRACNV
!   If IFIN=0, terminates after one step DS
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT, ONLY: HGTSFC ! Surface altitude [km]
    USE PHYADJ_DAT, ONLY: RADCRV ! Local radius of curvature [km]
    USE PHYCON_DAT, ONLY: DG2RAD ! [rad]/[deg] conversion factor
!
! SUBROUTINES
    USE GRADVS_SUB ! Derivatives for 2-D ray-tracing
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4),    INTENT(IN)    :: HGTINI ! Initial altitude [km]
    REAL(R4),    INTENT(IN)    :: PSIINI ! Initial LOS angle [deg]
    REAL(R4),    INTENT(IN)    :: ZENINI ! Initial Zenith Angle [deg]
    INTEGER(I4), INTENT(IN)    :: IFIN   ! Termin. parameter (1=HGT,2=PSI,3=ZEN)
    REAL(R4),    INTENT(INOUT) :: S      ! Init Path incr. (I) or Cum. path (O)
    REAL(R4),    INTENT(INOUT) :: HGTFIN ! Final height [km]
    REAL(R4),    INTENT(INOUT) :: PSIFIN ! Final LOS angle [deg]
    REAL(R4),    INTENT(INOUT) :: ZENFIN ! Final Zenith Angle [deg]
!
! LOCAL CONSTANTS
    REAL(R8), PARAMETER :: PTOL = 1.0D-6 ! Toler. [rads] to match PFIN  
    REAL(R8), PARAMETER :: RTOL = 1.0D-3 ! Toler. [km] to match RFIN  
    REAL(R8), PARAMETER :: TTOL = 1.0D-6 ! Toler. [rads] to match TFIN  
!
! LOCAL VARIABLES
    LOGICAL  :: SURFCE ! T=stopped at surface, F=stopped at reqd.limit
    REAL(R8) :: DFIN   ! Difference from target value of P,R,T
    REAL(R8) :: DPDS, DPDS1, DPDS2, DPDS3, DPDS4 ! d(P)/d(S)
    REAL(R8) :: DRDS, DRDS1, DRDS2, DRDS3, DRDS4 ! d(R)/d(S)
    REAL(R8) :: DS     ! Path length increment [km]
    REAL(R8) :: DS2    ! DS/2
    REAL(R8) :: DSFC   ! Distance above surface [km]
    REAL(R8) :: DSSUM  ! Cumulatived path [km]
    REAL(R8) :: DTDS, DTDS1, DTDS2, DTDS3, DTDS4 ! d(T)/d(S)
    REAL(R8) :: P      ! D.P. horizontal angle [rad]
    REAL(R8) :: PFIN   ! Termination value of P (if IFIN=2)
    REAL(R8) :: PINT   ! P evaluated at intermediate points
    REAL(R8) :: R      ! Radius coordinate of path [km]
    REAL(R8) :: RFIN   ! Termination value of R (if IFIN=1)
    REAL(R8) :: RINT   ! R evaluated at intermediate points
    REAL(R8) :: RSFC   ! Radius of curvature at earth's surface
    REAL(R8) :: T      ! D.P. zenith angle [rad]
    REAL(R8) :: TFIN   ! Termination value of T (if IFIN=3)
    REAL(R8) :: TINT   ! T evaluated at intermediate points
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  SURFCE = .FALSE.
! Ensure that any ray starting at tangent point or upwards never meets surface
  IF ( ABS ( ZENINI ) .LE. 90.0 ) THEN
    RSFC = 0.0D0
  ELSE
    RSFC = DBLE ( HGTSFC ) + RADCRV 
  END IF
  DS = DBLE ( S ) 
  DSSUM = 0.0D0
  DFIN  = 0.0D0            ! Not necessary but avoids warnings
!
! Look up atmospheric parameters for (z,psi)_INI
  R = RADCRV + DBLE ( HGTINI )
  P = DBLE ( PSIINI * DG2RAD )
  T = DBLE ( ZENINI * DG2RAD )
!
! Repeat from here for each path increment
  DO                      
    DSFC = R - RSFC 
    IF ( IFIN .NE. 0 .AND. ABS ( DSFC ) .LE. RTOL ) THEN
      SURFCE = .TRUE.
      EXIT                
    END IF
    SELECT CASE ( IFIN ) 
    CASE ( 1 ) 
      RFIN = RADCRV + DBLE ( HGTFIN )
      DFIN = R - RFIN 
      IF ( ABS ( DFIN ) .LE. RTOL ) EXIT
    CASE ( 2 ) 
      PFIN = DBLE ( PSIFIN * DG2RAD )
      DFIN = P - PFIN 
      IF ( ABS ( DFIN ) .LE. PTOL ) EXIT
    CASE ( 3 ) 
      TFIN = DBLE ( ZENFIN * DG2RAD )
      DFIN = T - TFIN 
      IF ( ABS ( DFIN ) .LE. TTOL ) EXIT
    CASE ( 0 ) 
      CONTINUE
    CASE DEFAULT 
      STOP 'F-RAYGRA: Logical error'
    END SELECT
! 
! Repeat from here if modifying step size
    DO 
      DS2 = DS * 0.5D0
! Runge-Kutta 4th order numerical integration increments y according to 
!  y_(i+1) = y_i + ( dy/dx_1 + 2*dy/dx_2 + 2*dy/dx_3 + dy/dx_4 )*dx/6
! where derivatives dy/dx_i are evaluated at: 
!  dy/dx_1: y_i 
!  dy/dx_2: y_i + (dy/dx_1)*dx/2
!  dy/dx_3: y_i + (dy/dx_2)*dx/2 
!  dy/dx_4: y_i + (dy/dx_3)*dx
      CALL GRADVS ( P, R, T, DPDS1, DRDS1, DTDS1 )      
      PINT = P + DPDS1 * DS2
      RINT = R + DRDS1 * DS2
      TINT = T + DTDS1 * DS2
      CALL GRADVS ( PINT, RINT, TINT, DPDS2, DRDS2, DTDS2 )   
      PINT = P + DPDS2 * DS2              
      RINT = R + DRDS2 * DS2 
      TINT = T + DTDS2 * DS2 
      CALL GRADVS ( PINT, RINT, TINT, DPDS3, DRDS3, DTDS3 )   
      PINT = P + DPDS3 * DS
      RINT = R + DRDS3 * DS
      TINT = T + DTDS3 * DS
      CALL GRADVS ( PINT, RINT, TINT, DPDS4, DRDS4, DTDS4 )      
      DPDS = ( DPDS1 + DPDS4 + 2.0D0 * ( DPDS2 + DPDS3 ) ) / 6.0D0
      DRDS = ( DRDS1 + DRDS4 + 2.0D0 * ( DRDS2 + DRDS3 ) ) / 6.0D0
      DTDS = ( DTDS1 + DTDS4 + 2.0D0 * ( DTDS2 + DTDS3 ) ) / 6.0D0
!      
! Putting these lines in gives 1st order (Euler) integration
!      dpds = dpds1
!      drds = drds1
!      dtds = dtds1
!
! If moving further under surface or further away from termination criterion, 
! set smaller stepsize DS
      IF ( IFIN .EQ. 0 ) EXIT  ! no iterating
      IF ( DSFC + DRDS * DS .LT. - RTOL ) THEN       ! Below surface
        DS = - DSFC / DRDS 
      ELSE IF ( IFIN .EQ. 1 .AND. &
              ABS ( DFIN + DRDS * DS ) .GT. ABS ( DFIN ) ) THEN
        DS = -DS * DFIN / ( DFIN + DRDS * DS )
      ELSE IF ( IFIN .EQ. 2 .AND. &
              ABS ( DFIN + DPDS * DS ) .GT. ABS ( DFIN ) ) THEN
        DS = -DS * DFIN / ( DFIN + DPDS * DS )
      ELSE IF ( IFIN .EQ. 3 .AND. &
                ABS ( DFIN + DTDS * DS ) .GT. ABS ( DFIN ) ) THEN
        DS = -DS * DFIN / ( DFIN + DTDS * DS )
      ELSE
        EXIT   ! Increment OK
      END IF
    END DO
! Add increment with same step size
    P = P + DPDS * DS
    R = R + DRDS * DS
    T = T + DTDS * DS
    DSSUM = DSSUM + DS
    IF ( IFIN .EQ. 0 ) EXIT   ! IFIN=0 means one iteration only
  END DO                      ! Loop for next increment
!
! Set final values to last iteration unless original term.criterion reached.
  IF ( IFIN .EQ. 0 ) THEN
    HGTFIN = SNGL ( R - RADCRV )
    PSIFIN = SNGL ( P ) / DG2RAD 
    ZENFIN = SNGL ( T ) / DG2RAD
  ELSE IF ( SURFCE ) THEN
    HGTFIN = HGTSFC
    PSIFIN = SNGL ( P ) / DG2RAD 
    ZENFIN = SNGL ( T ) / DG2RAD
  ELSE 
    IF ( IFIN .NE. 1 ) HGTFIN = SNGL ( R - RADCRV )
    IF ( IFIN .NE. 2 ) PSIFIN = SNGL ( P ) / DG2RAD 
    IF ( IFIN .NE. 3 ) ZENFIN = SNGL ( T ) / DG2RAD
  END IF
  S = SNGL ( DSSUM )
!
END SUBROUTINE RAYGRA
END MODULE RAYGRA_SUB

