MODULE CHISHP_SUB
CONTAINS
SUBROUTINE CHISHP ( DWNO, ABSORP )
!
! VERSION
!   20APR23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Calculate Voigt Line shape allowing for chi-factor
!   Called by LINSHP if CHI shape selected.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ADJCOM_DAT ! Path-adjusted line data
    USE FLGCOM_DAT, ONLY: MIXFLG ! T = use line-mixing
    USE PHYCON_DAT, ONLY: PI
!
! SUBROUTINES
    USE CHICO2_FNC ! Chi-factor for CO2 lines
    USE HUMLCK_SUB ! Calculate Humlicek complex prob.function
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(IN)  :: DWNO(:)   ! Array of Wavenumbers [/cm]
    REAL(R4), INTENT(OUT) :: ABSORP(:) ! Absorption 
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: SQRLN2 = SQRT ( LOG ( 2.0 ) ) !
    REAL(R4), PARAMETER :: SL2RPI = SQRT ( LOG ( 2.0 ) / PI )
    REAL(R4), PARAMETER :: CHILIM = 10.0 ! Inner lim. [cm-1] chi-fac only calc.
    REAL(R4), PARAMETER :: MIXLIM =  8.0 ! Outer lim. [cm-1] line-mix only calc.
!
! LOCAL VARIABLES
    INTEGER(I4) :: IWNO          ! Wavenumber array counter
    INTEGER(I4) :: NWNO          ! No. Wno points to be evaluted
    REAL(R4)    :: CHIFAC(SIZE(DWNO)) ! Chi factor
    REAL(R4)    :: CHIMUL        ! Weighting for Chi-factor term
    REAL(R4)    :: DELWNO        ! Wno difference from line centre
    REAL(R4)    :: H0            ! REPWID * Line strength
    REAL(R4)    :: REPWID        ! sqrt(ln2) / Doppler width
    REAL(R4)    :: X(SIZE(DWNO)) ! scaled wno from line ctr (x-coord for Hum.)
    REAL(R4)    :: Y             ! REPWID * Lorentz width (y-coord for Humlicek)
    COMPLEX(R4) :: V(SIZE(DWNO)) ! Complex Voigt vector
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NWNO = SIZE ( DWNO )
!
! Currently, Chi-factor only tabulated for CO2
  CHIFAC = CHICO2 ( TEMADJ, PREADJ, PPAADJ, WNOADJ, DWNO )
!
! Calculate Voigt line shape
  REPWID = SQRLN2 / DOPADJ
  H0 = SL2RPI * STRADJ / DOPADJ
  Y = WIDADJ * REPWID
  X =  SNGL ( DWNO - WNOADJ ) * REPWID
  CALL HUMLCK ( X, Y, V )
!
  IF ( MIXFLG .AND. YMXADJ .NE. 0.0 ) THEN
!
! Compute absorption due to line mixing combined with chi-factor
    DO IWNO = 1, NWNO 
      DELWNO = SNGL ( ABS ( DWNO(IWNO) - WNOADJ ) )
      IF ( DELWNO .LE. MIXLIM ) THEN    ! Line-mixing only
        ABSORP(IWNO) = H0 * ( REAL(V(IWNO)) + YMXADJ * AIMAG(V(IWNO)) )
      ELSE IF ( DELWNO .GE. CHILIM ) THEN  ! Chi-factor only
        ABSORP(IWNO) = H0 * REAL(V(IWNO)) * CHIFAC(IWNO) 
      ELSE
        CHIMUL =  ( DELWNO - MIXLIM ) / ( CHILIM - MIXLIM )           
        ABSORP(IWNO) = CHIMUL * H0 *  REAL(V(IWNO)) * CHIFAC(IWNO) &
           + (1.0-CHIMUL) * H0 * ( REAL(V(IWNO)) +  YMXADJ * AIMAG(V(IWNO)) )
      END IF
    END DO
  ELSE                            
!
! Compute absorption due chi-factor alone
    ABSORP = H0 * REAL(V) * CHIFAC
  END IF
!
END SUBROUTINE CHISHP
END MODULE CHISHP_SUB

