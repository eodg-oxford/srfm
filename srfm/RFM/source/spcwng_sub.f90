MODULE SPCWNG_SUB
CONTAINS
SUBROUTINE SPCWNG ( IWID )
!
! VERSION
!   04SEP24 AD Checked
!   01MAY17 AD F90 conversion of rfmwng.for. Checked.
!
! DESCRIPTION    
!   Interpolate wide mesh absorption across fine mesh
!   Called by RFMSPC for each widemesh interval.
!   The accumulated line-wing contributions for each wide-mesh interval are
!   interpolated quadratically or inverse-quadratically to the fine-mesh grid.
!
!   Fitted quadratic is Y = Y0 + B*U + A*U*(U-1)
!   where B = Y2 - Y0 and A = 2 * ( Y0 + Y2 - 2 * Y1 )
!
!   The default option is to try to fit an `inverse quadratic', ie fit a 
!   quadratic to the reciprocals of the 3 interpolation points (which gives a
!   better fit for Lorentz line wings). 
!   But this inversion can only be used safely if certain criteria are met:
!   (1) all three fitted points must have values > ABSMIN
!   (2) there has to be a significant curvature component ( U**2 coeff.)
!   (3) the max/min has to lie outside the widemesh interval (this is to 
!       prevent unrealistic large absorption peaks being created based on
!       three low-absorption points)
!   If any of these tests fail, a normal quadratic fit is applied directly to
!   the interpolation points. The normal quadratic is *ALWAYS* applied if the
!   QAD option is enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE FINCOM_DAT ! Finemesh data
    USE GASCOM_DAT ! Molecule and isotope data
    USE STTCOM_DAT ! Widemesh statistics
    USE WIDCOM_DAT ! Widemesh data
    USE FLGCOM_DAT, ONLY: QADFLG, WIDFLG ! Option flags
!                  
  IMPLICIT NONE
!
! ARGUMENTS      
    INTEGER(I4), INTENT(IN) :: IWID ! Index of current wide mesh interval [IPW]
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: ABSMIN = 1.0E-30 ! Min value for inverse quadratic
    REAL(R4), PARAMETER :: CRVMIN = 1.0E-3  ! Min curvature for inverse quad.fit
!
! LOCAL VARIABLES
    LOGICAL     :: USEIQD      ! T=apply inverse quadratic fit algorithm
    INTEGER(I4) :: ICLC        ! Index of calculated path
    INTEGER(I4) :: IFIN        ! Counter for fine-mesh intervals
    INTEGER(I4) :: IGAS        ! Index of absorber in GAS
    INTEGER(I4) :: ILBL        ! Counter for LBL calculated paths
    REAL(R4)    :: A, B        ! Differences between Yn values
    REAL(R4)    :: AI,BI       ! Differences between YIn values
    REAL(R4)    :: F           ! Value of quadratic function 
    REAL(R4)    :: U           ! Scaled Value (0:1) of Fine-Mesh points
    REAL(R4)    :: UPEAK       ! Position of quadratic peak in scaled grid units
    REAL(R4)    :: Y0,Y1,Y2    ! Values of quadratic functions at Half-Mesh pts
    REAL(R4)    :: YI0,YI1,YI2 ! Values of inv.quadratic fns at Half-Mesh pts
    REAL(R8)    :: WNOOFF      ! Wavenumber [cm-1] at start of Widemesh Intvl
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  WNOOFF = WNOWID(IWID-1) 
  DO ILBL = 1, NLBL
    ICLC = IDXLBL(ILBL) 
    IGAS = CLC(ICLC)%IGS
    Y0 = ABSWID(1,IWID,ILBL)
    Y1 = ABSWID(2,IWID,ILBL)
    Y2 = ABSWID(3,IWID,ILBL)
    B  = Y2 - Y0
    A  = 2.0 * ( Y0 + Y2 - 2.0*Y1 )
!
    USEIQD = .NOT. QADFLG .AND. MIN ( Y0, Y1, Y2 ) .GT. ABSMIN 
    IF ( USEIQD ) THEN                   ! test conditions for using IQD fit
      YI0 = 1.0 / Y0
      YI1 = 1.0 / Y1
      YI2 = 1.0 / Y2
      BI  = YI2 - YI0
      AI  = 2.0 * ( YI0 + YI2 - 2.0*YI1 )
      USEIQD = ( ABS ( AI ) .GT. CRVMIN * ABS ( BI ) )
      IF ( USEIQD ) THEN
        UPEAK = 0.5 * ( AI - BI ) / AI
        USEIQD = ( UPEAK .LE. 0.0 .OR. UPEAK .GE. 1.0 )
      END IF
    END IF
!
    IF ( USEIQD ) THEN
      DO IFIN = 1, NFIN                     ! Loop over fine mesh frequency
        U = SNGL ( ( WNOFIN(IFIN) - WNOOFF ) / DELWID )
        F = YI0 + U * ( BI + AI * (U - 1.0) ) ! guaranteed .GE. ABSMIN
        ABSFIN(IFIN,ICLC) = 1.0 / F
      END DO
    ELSE                                    ! Normal quadratic fit
      DO IFIN = 1, NFIN                     ! Loop over fine mesh frequency
        U = SNGL ( ( WNOFIN(IFIN) - WNOOFF ) / DELWID )
        ABSFIN(IFIN,ICLC) = MAX ( Y0 + U * ( B + A*(U-1.0) ), 0.0 )
      END DO
      IF ( WIDFLG ) STT(IWID,IGAS)%NQD = STT(IWID,IGAS)%NQD + 1
    END IF
!
    IF ( GAS(IGAS)%NTE ) THEN                      ! always normal quadratic fit
      Y0 = CNTWID(1,IWID,ILBL)
      Y1 = CNTWID(2,IWID,ILBL)
      Y2 = CNTWID(3,IWID,ILBL)
      B  = Y2 - Y0
      A  = 2.0 * ( Y0 + Y2 - 2.0*Y1 )
      DO IFIN = 1, NFIN
        U = SNGL ( ( WNOFIN(IFIN) - WNOOFF ) / DELWID )
        CNTFIN(IFIN,ICLC) = MAX ( Y0 + U * ( B + A*(U-1.0) ), 0.0 )
      END DO
    END IF
  END DO
!   
END SUBROUTINE SPCWNG
END MODULE SPCWNG_SUB
