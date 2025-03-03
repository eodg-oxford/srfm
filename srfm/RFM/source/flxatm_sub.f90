MODULE FLXATM_SUB
CONTAINS
SUBROUTINE FLXATM ( DOWN, RQAD, ILEV ) 
!
! VERSION
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
    USE ATMCOM_DAT ! Atmospheric profile data
    USE FULCOM_DAT ! Full grid data
    USE LEVCOM_DAT ! Intermediate output levels
    USE QADCOM_DAT ! Gaussian quadrature data
    USE FINCOM_DAT, ONLY: NFIN           ! Finemesh data
    USE FLGCOM_DAT, ONLY: COOFLG, NADFLG ! Option flags
!
! SUBROUTINES
    USE FLXLAY_SUB ! Radiative flux calculation through a layer
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,     INTENT(IN)    :: DOWN            ! T=downward path, F=upward
    REAL(R8),    INTENT(INOUT) :: RQAD(NFIN,NQAD) ! Radiances for quad. paths
    INTEGER(I4), &
       OPTIONAL, INTENT(IN)    :: ILEV            ! Perturbation level
!
! LOCAL VARIABLES
    LOGICAL     :: LSAVE           ! Save intermediate outputs
    INTEGER(I4) :: IATM            ! Counter for atmospheric levels
    INTEGER(I4) :: IATM1, IATM2    ! Start/end levels for path
    INTEGER(I4) :: IDIR            ! Direction of integration
    INTEGER(I4) :: IPTB            ! Atmospheric level for ILEV
    INTEGER(I4) :: IQAD            ! Counter for quadrature points
    INTEGER(I4) :: ITAN            ! Counter for output levels
    INTEGER(I4) :: IWGT            ! Counter for cooling rate weights
    REAL(R8)    :: OPT(NFIN)       ! Cumulative optical path
    REAL(R8)    :: OPTLAY(NFIN)    ! Optical depth of single atm layer
    REAL(R8)    :: TQAD(NFIN,NQAD) ! Transmittances for quadrature paths
!       
! EXECUTABLE CODE -------------------------------------------------------------
!
! Downward path 
  IF ( DOWN ) THEN
    IATM1 = NATM 
    IATM2 = 1
    IDIR = -1
  ELSE
    IATM1 = 1
    IATM2 = NATM
    IDIR = +1
  END IF
!
  IF ( PRESENT ( ILEV ) ) THEN
    IPTB = LEV(ILEV)%IAT
  ELSE
    IPTB = 0
  END IF
  OPT = 0.0D0

! Don't store outputs if NADFLG on downward path
  LSAVE = .NOT. ( NADFLG .AND. IDIR .EQ. -1 ) ! 

  DO IATM = IATM1, IATM2, IDIR               ! Loop over levels
!
! Lev=bottom of layer so on downward path calc layer contribution first
    IF ( IDIR .EQ. -1 ) THEN 
      CALL FLXLAY ( IATM, IDIR, (IATM.EQ.IPTB), XQAD, RQAD, OPTLAY ) 
      OPT = OPT + OPTLAY
    END IF

    IF ( ITNATM(IATM) .GT. 0 .AND. LSAVE ) THEN
      ITAN = ITNATM(IATM)
      IF ( IPTB .GT. 0 ) ITAN = ITNLEV(ITAN,ILEV) 
      RADFUL(IFUL1:IFUL2,ITAN) = RADFUL(IFUL1:IFUL2,ITAN) + &
                                 IDIR *  MATMUL(RQAD,WQAD) ! -ve rad if down.
      OPTFUL(IFUL1:IFUL2,ITAN) = OPT
! Only calculate transmittance for unperturbed case
      IF ( IPTB .EQ. 0 ) THEN
        DO IQAD = 1, NQAD
          TQAD(:,IQAD) = EXP ( - OPT / XQAD(IQAD) ) * RPIQAD
        END DO
        TRAFUL(IFUL1:IFUL2,ITAN) = MATMUL ( TQAD, WQAD ) 
      END IF
    END IF
!
    IF ( COOFLG .AND. LSAVE ) THEN
      DO IWGT = 1, NCOATM(IATM)
        ITAN = ICOATM(IWGT,IATM) 
        IF ( IPTB .GT. 0 ) ITAN = ITNLEV(ITAN,ILEV)
        COOFUL(IFUL1:IFUL2,ITAN) = COOFUL(IFUL1:IFUL2,ITAN) + &
                                 IDIR * MATMUL(RQAD,WQAD) * WCOATM(IWGT,IATM)
      END DO
    END IF
!
! Lev=bottom of layer so on upward path calc layer contribution last
    IF ( IDIR .EQ. 1 ) THEN 
      CALL FLXLAY ( IATM, IDIR, (IATM.EQ.IPTB), XQAD, RQAD, OPTLAY ) 
      OPT = OPT + OPTLAY
    END IF

  END DO
!
END SUBROUTINE FLXATM
END MODULE FLXATM_SUB
