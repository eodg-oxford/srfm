MODULE FLXATM_SUB
CONTAINS
SUBROUTINE FLXATM ( DOWN, RQAD, IPTB, ITNOFF ) 
!
! VERSION
!   29MAY26 AD Bug#57: Integrate from atm level IATSFC rather than 1
!   15MAY26 AD Bug#56: Ensure FLXRAD defined at start of upward path
!   01AUG25 AD Remove LEVCOM. Use IDXARR to determine if output level.
!   24APR24 AD Checked.
!   05MAR19 AD Original. Extracted from SPCFLX. Checked.
!
! DESCRIPTION
!   Atmospheric flux calculation
!   Called by SPCFLX, RADMTX.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FULCOM_DAT ! Full grid data
    USE QADCOM_DAT ! Gaussian quadrature data
    USE TANCOM_DAT ! Tangent path data
    USE ATMCOM_DAT, ONLY: NATM, iatsfc ! No. Atmospheric profile levels
    USE FINCOM_DAT, ONLY: NFIN ! No. finemesh grid points
    USE FLGCOM_DAT, ONLY: COOFLG, NADFLG ! Option flags
!
! SUBROUTINES
    USE COOLAY_SUB ! Cooling rate calculation for atmospheric level
    USE FLXLAY_SUB ! Radiative flux calculation through a layer
    USE IDXARR_GEN ! Index of value within an array, else 0
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,     INTENT(IN)    :: DOWN            ! T=downward path, F=upward
    REAL(R8),    INTENT(INOUT) :: RQAD(NFIN,NQAD) ! Radiances for quad. paths
    INTEGER(I4), INTENT(IN)    :: IPTB            ! Perturbation level
    INTEGER(I4), INTENT(IN)    :: ITNOFF          
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: MAXWGT = 5 !  Max.no o/p levels affected by atm#
!                                           Quadratic fit implies a maximum=5
! LOCAL VARIABLES
    LOGICAL     :: LSAVE           ! Save intermediate outputs
    INTEGER(I4) :: IATM            ! Counter for atmospheric levels
    INTEGER(I4) :: IATM1, IATM2    ! Start/end levels for path
    INTEGER(I4) :: IDIR            ! Direction of integration
    INTEGER(I4) :: IQAD            ! Counter for quadrature points
    INTEGER(I4) :: ITAN            ! Counter for output levels
    INTEGER(I4) :: JTAN            ! Index for secondary ray path
    REAL(R8)    :: OPT(NFIN)       ! Cumulative optical path
    REAL(R8)    :: OPTLAY(NFIN)    ! Optical depth of single atm layer
    REAL(R8)    :: TQAD(NFIN,NQAD) ! Transmittances for quadrature paths
    REAL(R8)    :: FLXRAD(NFIN)    ! Integrated radiance flux
!       
! EXECUTABLE CODE -------------------------------------------------------------
!
! Downward path 
  IF ( DOWN ) THEN
    IATM1 = NATM 
    IATM2 = IATSFC
    IDIR = -1
  ELSE
    IATM1 = IATSFC
    IATM2 = NATM
    IDIR = +1
  END IF
!
  OPT = 0.0D0
!
! Don't store outputs if NADFLG on downward path
  LSAVE = .NOT. ( NADFLG .AND. IDIR .EQ. -1 ) 

  DO IATM = IATM1, IATM2, IDIR               ! Loop over levels
!
! Lev=bottom of layer so on downward path calc layer contribution first
    IF ( IDIR .EQ. -1 ) THEN 
      CALL FLXLAY ( IATM, IDIR, (IATM.EQ.IPTB), XQAD, RQAD, OPTLAY ) 
      OPT = OPT + OPTLAY
    END IF
    FLXRAD = IDIR * MATMUL(RQAD,WQAD)    ! -ve rad if down. Bug#56

    ITAN = IDXARR ( TAN(1:NTAN)%IAT, IATM ) 
    IF ( ITAN .GT. 0 .AND. LSAVE ) THEN
      JTAN = ITAN + ITNOFF
      RADFUL(IFUL1:IFUL2,JTAN) = RADFUL(IFUL1:IFUL2,JTAN) + FLXRAD
 ! Only calculate transmittance for unperturbed case
      IF ( IPTB .EQ. 0 ) THEN
        OPTFUL(IFUL1:IFUL2,ITAN) = OPT  
       DO IQAD = 1, NQAD
          TQAD(:,IQAD) = EXP ( - OPT / XQAD(IQAD) ) * RPIQAD
        END DO
        TRAFUL(IFUL1:IFUL2,ITAN) = MATMUL ( TQAD, WQAD ) 
      END IF
    END IF
! Add contribution of this layer to cooling rates
    IF ( COOFLG .AND. LSAVE ) CALL COOLAY ( IATM, FLXRAD, ITNOFF )
!
! Lev=bottom of layer so on upward path calc layer contribution last
    IF ( IDIR .EQ. 1 ) THEN 
      CALL FLXLAY ( IATM, IDIR, (IATM.EQ.IPTB), XQAD, RQAD, OPTLAY ) 
      OPT = OPT + OPTLAY
      FLXRAD = IDIR * MATMUL(RQAD,WQAD)   
    END IF

  END DO
!
END SUBROUTINE FLXATM
END MODULE FLXATM_SUB
