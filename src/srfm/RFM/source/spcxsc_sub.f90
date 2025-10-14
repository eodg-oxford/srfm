MODULE SPCXSC_SUB
CONTAINS
SUBROUTINE SPCXSC ( IWID ) 
!
! VERSION
!   06APR24 AD Checked.
!   08FEB19 AD Bug#16: Add NTRI argument to TRIINT. Checked.
!   01MAY17 AD F90 conversion of rfmxsc.for. Checked.
!
! DESCRIPTION    
!   Absorption of x/s molecules
!   Called by RFMSPC for each widemesh interval
!
! VARIABLE KINDS
    USE KIND_DAT  
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE FINCOM_DAT ! Finemesh data
    USE GASCOM_DAT ! Molecule and isotope data
    USE STTCOM_DAT ! Widemesh statistics
    USE XSCCOM_DAT ! Tabulated Cross-Section data
    USE FLGCOM_DAT, ONLY: WIDFLG ! T=collect widemesh statistics
    USE PHYCON_DAT, ONLY: AVOG ! Avogadro's constant [kmole-1]
!
! SUBROUTINES
    USE TRIINT_SUB ! 2-D interpolation of irregular grid using triangulation
    USE XSCINT_FNC ! Spectral interpolation of XSC dataset
!                  
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IWID ! Widemesh interval
!
! LOCAL VARIABLES
    INTEGER(I4) :: ICLC      ! Index of calculated path
    INTEGER(I4) :: IDXLKP(3) ! Indices of tabulations use for (p,T) interp
    INTEGER(I4) :: IGAS      ! Index of absorber in GASCOM
    INTEGER(I4) :: IXSC      ! Counter for X/S data sets
    REAL(R4)    :: KXS(NFIN) ! x/s data interpolated to fine grid
    REAL(R4)    :: PFIT      ! Fitted pressure coordinate value (dummy)
    REAL(R4)    :: TFIT      ! Fitted temperature coordinate value (dummy)
    REAL(R4)    :: WGTLKP(3) ! Weights of tabulations use for (p,T) interp
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( NXSC .EQ. 0 ) RETURN
!
! Loop over cross-section data file
  DO IXSC = 1, NXSC
    IF ( XSC(IXSC)%WNL .GT. WN2FIN .OR. XSC(IXSC)%WNU .LE. WN1FIN ) CYCLE
    IGAS = XSC(IXSC)%IGS
    IF ( .NOT. GAS(IGAS)%XSC ) CYCLE  ! Set FALSE if LUT used for this spc.range
    DO ICLC = 1, NCLC
      IF ( CLC(ICLC)%IGS .NE. IGAS ) CYCLE
      CALL TRIINT ( CLC(ICLC)%TEM, CLC(ICLC)%PRE*760.0, XSC(IXSC)%NXT, & 
                    XSC(IXSC)%TEM, XSC(IXSC)%PRE, &
                    XSC(IXSC)%NTRI, XSC(IXSC)%ITRI, &
                    IDXLKP, WGTLKP, TFIT, PFIT )
      KXS = XSCINT ( IXSC, IDXLKP(1), NFIN, WNOFIN ) * WGTLKP(1) + & 
            XSCINT ( IXSC, IDXLKP(2), NFIN, WNOFIN ) * WGTLKP(2) + &
            XSCINT ( IXSC, IDXLKP(3), NFIN, WNOFIN ) * WGTLKP(3)
      ABSFIN(:,ICLC) = CLC(ICLC)%AMT * AVOG * KXS
      IF ( WIDFLG ) STT(IWID,IGAS)%IXS = 1
    END DO
  END DO
!
END SUBROUTINE SPCXSC
END MODULE SPCXSC_SUB
