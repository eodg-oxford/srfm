MODULE DIMPLV_SUB
CONTAINS
SUBROUTINE DIMPLV 
!
! VERSION
!   07OCT24 AD Checked.
!   01MAY17 AD F90 conversion of part of dimpre.for. Tested.
!
! DESCRIPTION
!   Set TAB p-axis to internal PRE ATM profile
!   Called by DRVDIM 
!   Called if TAB flag and PLV present in *DIM section of driver file.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT, ONLY: NATM, PREATM  ! Atmospheric pressure profile
    USE TABCOM_DAT, ONLY: NPTAB, PAXTAB ! p-axis for tabulated abs. coeffs.
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IATM ! Profile level counter
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  CALL WRTLOG ( 'I-DIMPLV: Setting TAB p-axis to pressure levels ' // &
                'of internal .atm profile' )
  NPTAB = NATM           ! NATM level values
  ALLOCATE ( PAXTAB(NPTAB) ) 
  PAXTAB = PREATM
!
! Print values to log file
  CALL WRTLOG ( 'p-axis values: ', .TRUE. )
  DO IATM = 1, NPTAB
    CALL WRTLOG ( C9REAL ( PAXTAB(IATM) ), .TRUE. )
  END DO
  CALL WRTLOG ( '', .FALSE. ) 
!
END SUBROUTINE DIMPLV
END MODULE DIMPLV_SUB
