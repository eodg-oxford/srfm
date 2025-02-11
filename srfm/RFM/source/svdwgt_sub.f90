MODULE SVDWGT_SUB
CONTAINS
SUBROUTINE SVDWGT ( FINISH, IGAS, NLNP, LNP1, DLNP, NTEM, TEM1, DTEM, &
                    PRE, TEM, IDXLKP, WGTLKP )
! VERSION
!   02JAN24 AD Checked.
!   01JUL17 AD F90 conversion of part of svdpth.for. Checked.
! 
! DESCRIPTION    
!   Calculate SVD-LUT p,T-axis interpolation weights for path
!   Called by SVDPTH
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE PHYCON_DAT, ONLY: ATMB ! Std Atmos. Pressure [mb]
!
! SUBROUTINES
    USE SVDWRN_SUB ! Warn if SVD-LUT p,T axis limits exceeded
!                  
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,     INTENT(IN)  :: FINISH    ! T=last path, F=continue
    INTEGER(I4), INTENT(IN)  :: IGAS      ! Index of absorber
    INTEGER(I4), INTENT(IN)  :: NLNP      ! No. -lnp axis points
    REAL(R4),    INTENT(IN)  :: LNP1      ! 1st -lnp axis point
    REAL(R4),    INTENT(IN)  :: DLNP      ! -lnp axis interval
    INTEGER(I4), INTENT(IN)  :: NTEM      ! No. Tem axis points
    REAL(R4),    INTENT(IN)  :: TEM1      ! 1st Tem axis point
    REAL(R4),    INTENT(IN)  :: DTEM      ! Tem axis interval
    REAL(R4),    INTENT(IN)  :: PRE       ! Pressure [mb]
    REAL(R4),    INTENT(IN)  :: TEM       ! Temperature [K]
    INTEGER(I4), INTENT(OUT) :: IDXLKP(:) ! Interpolation points
    REAL(I4),    INTENT(OUT) :: WGTLKP(:) ! Interpolation weights
!
! LOCAL VARIABLES
    INTEGER(I4) :: ILKP       ! Counter for interpolation points
    INTEGER(I4) :: IP,IT      ! Indices of lower points in p,T-axes
    INTEGER(I4) :: JP,JT      ! Indices of upper points in p,T-axes
    INTEGER(I4) :: KP,KT      ! Counter for lower/upper points in p,T-axes
    INTEGER(I4) :: NP1LIM = 0 ! Counter: no.times p outside PAXTAB lower index 
    INTEGER(I4) :: NP2LIM = 0 ! Counter: no.times p outside PAXTAB upper index
    INTEGER(I4) :: NT1LIM = 0 ! Counter: no.times T outside TAXTAB lower index
    INTEGER(I4) :: NT2LIM = 0 ! Counter: no.times T outside TAXTAB upper index
    INTEGER(I4) :: NTOT   = 0 ! Total number of calculated paths for IGAS
    REAL(R4)    :: DP,DT      ! Fraction of p,T-axis interval above IP,IT
    REAL(R4)    :: LNP        ! -ln(Pressure [mb]) of path segment 
    REAL(R4)    :: PFACT      ! Interpolation weight for p-axis (DP or 1-DP)
    REAL(R4)    :: TFACT      ! Interpolation weight for T-axis (DT or 1-DT)
    REAL(R4)    :: XFACT      ! Product of interpolation weights PFACT,TFACT
    REAL(R4)    :: XP,XT      ! Position of path p,T in p,T-axis units
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  NTOT = 1
!
  LNP = - LOG ( PRE * ATMB )          ! convert [atm] to [mb]
  XP = ( LNP - LNP1 ) / DLNP + 1.0
  IP = INT ( XP )
  IF ( XP .LT. 1 ) THEN
    IP = 1
    DP = 0.0
    NP1LIM = NP1LIM + 1
  ELSE IF ( XP .GT. NLNP ) THEN
    IP = NLNP
    DP = 0.0
    NP2LIM = NP2LIM + 1
  ELSE
    DP = XP - IP
  END IF
  JP = MIN ( NLNP, IP + 1 )     ! MIN to allow for IP=NP
!
  XT = ( TEM - TEM1 ) / DTEM + 1.0
  IT = INT ( XT )
  IF ( XT .LT. 1 ) THEN
    IT = 1
    DT = 0.0
    NT1LIM = NT1LIM + 1
  ELSE IF ( XT .GT. NTEM ) THEN
    IT = NTEM 
    DT = 0.0
    NT2LIM = NT2LIM + 1
  ELSE
    DT = XT - IT
  END IF
  JT = MIN ( NTEM, IT + 1 )     ! MIN to allow for IT=NTSVD
!
! Calculate weights for lower and upper points in p,T axes
  IDXLKP = 1
  WGTLKP = 0.0
  ILKP = 0
  TFACT = 1.0 - DT
  DO KT = IT, JT
    PFACT = 1.0 - DP
    DO KP = IP, JP
      XFACT = TFACT * PFACT 
      ILKP = ILKP + 1
      WGTLKP(ILKP) = XFACT
      IDXLKP(ILKP) = (KT-1) * NLNP + KP
      PFACT = DP
    END DO
    TFACT = DT
  END DO                                
!
  IF ( .NOT. FINISH ) RETURN
!
  CALL SVDWRN ( IGAS, NLNP, DLNP, NTEM, DTEM, &
                NTOT, NP1LIM, NP2LIM, NT1LIM, NT2LIM ) 

! Reset counters
  NP1LIM = 0
  NP2LIM = 0
  NT1LIM = 0
  NT2LIM = 0
  NTOT = 0
!
END SUBROUTINE SVDWGT
END MODULE SVDWGT_SUB


