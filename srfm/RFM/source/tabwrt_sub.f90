MODULE TABWRT_SUB
CONTAINS
SUBROUTINE TABWRT ( FAIL, ERRMSG )
!
! VERSION
!   25JUN24 AD Bug#46. Use PTH%AMT to scale path amount.
!   04JUL23 AD Checked.
!   13OCT17 AD F90 conversion of wrttab.for. Checked.
!
! DESCRIPTION
!   Write RFM tabulated absorption coefficients
!   Called by SPCTAB for each widemesh interval if TAB flag selected.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FINCOM_DAT ! Finemesh data
    USE PTHCOM_DAT ! Path segment data
    USE TABCOM_DAT ! Axes for tabulated absorption coefficients
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error occurs
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: ZERLNK = -99.0 ! ln(k) value used for k < TINY(0.0)
!
! LOCAL VARIABLES
    INTEGER(I4) :: IFIN ! Fine mesh grid-point counter
    INTEGER(I4) :: IOFF ! Path offset index for different gases
    INTEGER(I4) :: IOS  ! Value of IOSTAT saved for error messages.
    INTEGER(I4) :: ITAB ! Counter for .tab files/absorbers
    INTEGER(I4) :: IX   ! Index for p,T,q elements
    INTEGER(I4) :: LUN  ! Logical Unit Number.
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! ifort compiler gives segmentation fault if WHERE applied to complete ABSFIN
! array for large array, so carry this out for each wavenumber point instead
  DO IFIN = 1, NFIN
    ABSFIN(IFIN,:) = ABSFIN(IFIN,:) * PTH(:)%AMT
    WHERE ( ABSFIN(IFIN,:) .GE. TINY(1.0) ) 
      ABSFIN(IFIN,:) = LOG ( ABSFIN(IFIN,:) )
    ELSEWHERE
      ABSFIN(IFIN,:) = ZERLNK
    ENDWHERE
  END DO
!          
  LUN = LUNTAB
  IOFF = 0
  DO ITAB = 1, NTAB
    IF ( BINTAB ) THEN
      DO IFIN = 1, NFIN
        WRITE ( LUN, IOSTAT=IOS, ERR=900 ) &
          WNOFIN(IFIN)*GHZTAB, ( ABSFIN(IFIN,IX), IX = IOFF+1, IOFF+NATAB ) 
      END DO
    ELSE
      DO IFIN = 1, NFIN
        WRITE ( LUN, *, IOSTAT=IOS, ERR=900 ) &
          WNOFIN(IFIN)*GHZTAB, ( ABSFIN(IFIN,IX), IX = IOFF+1, IOFF+NATAB )
      END DO
    END IF
    LUN  = LUN + 1
    IOFF = IOFF + NATAB
  END DO
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) & 
    'F-TABWRT: I/O failure on output file. IOSTAT=', IOS
!
END SUBROUTINE TABWRT
END MODULE TABWRT_SUB

