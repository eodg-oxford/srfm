MODULE REALUT_SUB
CONTAINS
SUBROUTINE REALUT ( ILUT, WNO, ZERABS, TAB, FAIL, ERRMSG )
!
! VERSION
!   31MAY24 AD Checked.
!   28MAR19 AD Bug#19: Check for WNO outside TAB data range. Checked.
!   01MAY17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   Interpolate data from LUT file for given wavenumber
!   Called by SPCLUT.
!   ZERABS is set TRUE (and TABWNO undefined) if all abs coeffs are zero, 
!   either due to interpolating between two records containing only zeros, 
!   or due to wavenumber being outside spectral range of tabulation.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LUTCOM_DAT ! TAB LUT data
    USE PHYCON_DAT, ONLY: GHZ2CM ! GHz to cm-1 conv.fac
!                  
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: ILUT   ! Index of LUT data
    REAL(R8),      INTENT(IN)  :: WNO    ! Wavenumber [cm-1] to be interpolated
    LOGICAL,       INTENT(OUT) :: ZERABS ! T=zero absorption
    REAL(R4),      INTENT(OUT) :: TAB(:) ! Interpolated ln(Abs.Coeff)
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: LNKZER = -80.0 ! Value of ln(k) treated as k=0
!                                           ln(k)=-80 equiv to k=1.8e-35
! LOCAL VARIABLES
    INTEGER(I4) :: IBUF   ! Counter over buffers (1:2)
    INTEGER(I4) :: IOS    ! Saved value of IOSTAT for error message
    REAL(R4)    :: DV     ! Spectral interpolation weight
    REAL(R8)    :: WNOUPP ! Upper wavenumber [cm-1] of two buffers
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
  IOS = 0        ! Initialise in case no data read on this call
!
! Possible that LUT is used outside range of tabulated data if spectral 
! convolution is applied, in which case assume zero absorption
  IF ( WNO .LT. LUT(ILUT)%WNL .OR. &
       WNO .GT. LUT(ILUT)%WNU      ) THEN
    ZERABS = .TRUE.
    RETURN
  END IF

  IBUF = LUT(ILUT)%IBU
  WNOUPP = LUT(ILUT)%WNO(IBUF)
  DO WHILE ( WNOUPP .LT. WNO ) 
    IBUF = 3 - IBUF
    IF ( LUT(ILUT)%BIN ) THEN
      READ ( LUT(ILUT)%LUN, IOSTAT=IOS, ERR=900 ) &
        WNOUPP, LUT(ILUT)%TAB(:,IBUF) 
    ELSE 
      READ ( LUT(ILUT)%LUN, *, IOSTAT=IOS, ERR=900 ) &
        WNOUPP, LUT(ILUT)%TAB(:,IBUF) 
    END IF
    IF ( LUT(ILUT)%GHZ ) WNOUPP = WNOUPP * GHZ2CM
    LUT(ILUT)%WNO(IBUF) = WNOUPP
    LUT(ILUT)%IBU = IBUF
  END DO
!
! Check if buffers are full of zero absorption values
  ZERABS = MAXVAL ( LUT(ILUT)%TAB ) .LT. LNKZER 
  IF ( ZERABS ) RETURN
!
! Check if WNO bracketed by loaded buffers
  DV = SNGL ( ( WNO - LUT(ILUT)%WNO(1) ) / &
              ( LUT(ILUT)%WNO(2) - LUT(ILUT)%WNO(1) ) )  
  ZERABS = ( DV .LT. 0.0 .OR. DV .GT. 1.0 )
  IF ( ZERABS ) RETURN
!
  TAB = (1.0-DV) * LUT(ILUT)%TAB(:,1) + DV * LUT(ILUT)%TAB(:,2)
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-REALUT: Input failure on TAB file. IOSTAT=', IOS
!
END SUBROUTINE REALUT
END MODULE REALUT_SUB

