MODULE TANMTX_SUB
CONTAINS
SUBROUTINE TANMTX
!
! VERSION
!   09SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Set matrix output levels
!   Called by DRVTAN if MTX flag is set
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LEVCOM_DAT ! Intermediate output levels 
    USE TANCOM_DAT ! Tangent path data
!
  IMPLICIT NONE 
!
! LOCAL VARIABLES
    INTEGER(I4) :: ILEV ! Level counter for matrix elements
    INTEGER(I4) :: ITAN ! Counter for output levels (1st matrix dimension)
    INTEGER(I4) :: JTAN ! Secondary level counter (2nd matrix dimension)
!
! EXECUTABLE CODE -------------------------------------------------------------
!    
  NLEV = NTAN
  MTAN = NTAN + NTAN**2
!
  ALLOCATE ( LEV(NLEV) ) 
  ALLOCATE ( ITNLEV(NTAN,NTAN) ) 
!
  JTAN = NTAN
  DO ILEV = 1, NLEV
    LEV(ILEV)%IAT = TAN(ILEV)%IAT
    LEV(ILEV)%HGT = TAN(ILEV)%HGT
    LEV(ILEV)%IDR = 0 
    DO ITAN = 1, NTAN 
      JTAN = JTAN + 1
      ITNLEV(ITAN,ILEV) = JTAN
    END DO
  END DO
!
END SUBROUTINE TANMTX
END MODULE TANMTX_SUB

