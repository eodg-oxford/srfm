MODULE XSCINT_FNC
CONTAINS
FUNCTION XSCINT ( IXSC, IXT, NWNO, WNOLST ) 
!
! VERSION
!   24JUN24 AD Checked.
!   01MAY17 AD F90 version. Checked.
! 
! DESCRIPTION    
!   Spectral interpolation of XSC dataset
!   Called by SPCXSC.
!
! VARIABLE KINDS
    USE KIND_DAT
!                  
! GLOBAL DATA
    USE XSCCOM_DAT ! Tabulated Cross-Section data
!
  IMPLICIT NONE
!
! ARGUMENTS      
    INTEGER(I4), INTENT(IN) :: IXSC      ! Index of x/s dataset
    INTEGER(I4), INTENT(IN) :: IXT       ! Index of (p,T) within dataset
    INTEGER(I4), INTENT(IN) :: NWNO      ! No.points in WNOLST
    REAL(R8),    INTENT(IN) :: WNOLST(:) ! List of reqd wavenumbers [cm-1]
!
! FUNCTION TYPE
    REAL(R4) :: XSCINT(NWNO) ! Function returns array same size as WNOLST
!
! LOCAL VARIABLES
    INTEGER(I4) :: IOF       ! Offset of data in XSC%ABS array
    INTEGER(I4) :: IPT(NWNO) ! Index of datapoint below reqd wavenumbers
    REAL(R4)    :: DPT(NWNO) ! Fraction of interval between tabulated wnos
    REAL(R4)    :: XPT(NWNO) ! Fractional index of WNOLST in tabulated data
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IOF = XSC(IXSC)%IOF(IXT)
  XPT = 1.0 + SNGL ( ( WNOLST - XSC(IXSC)%WN1(IXT) ) / XSC(IXSC)%DWN(IXT) )
  WHERE ( XPT .GE. 1.0 .AND. XPT .LT. XSC(IXSC)%NPT(IXT) ) 
    IPT = INT ( XPT )
    DPT = XPT - FLOAT ( IPT ) 
    XSCINT = XSC(IXSC)%ABS(IOF+IPT) * ( 1.0 - DPT ) + &
             XSC(IXSC)%ABS(IOF+IPT+1) * DPT         
  ELSEWHERE
    XSCINT = 0.0
  END WHERE
!
END FUNCTION XSCINT
END MODULE XSCINT_FNC
