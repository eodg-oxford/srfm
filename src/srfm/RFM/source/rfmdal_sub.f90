MODULE RFMDAL_SUB
CONTAINS
SUBROUTINE RFMDAL 
!
! VERSION
!   21MAR23 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Deallocate program level pointers
!   Called once by RFM
!
! VARIABLE KINDS
    USE KIND_DAT   ! Variable kinds
!
! GLOBAL DATA
    USE CIACOM_DAT ! Collision-induced absorption data
    USE GASCOM_DAT ! Molecule and isotope data
    USE ILSCOM_DAT ! Instrument Lineshape data
    USE XSCCOM_DAT ! Tabulated Cross-Section data
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: I ! general counter 
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO I = 1, NCIA
    DEALLOCATE ( CIA(I)%ABS, CIA(I)%WNO )
  END DO  
!
  DO I = 1, NGAS
    IF ( GAS(I)%NIS .GT. 0 ) DEALLOCATE ( GAS(I)%ISO, GAS(I)%WGT ) 
  END DO
!
  DO I = 1, NILS
    DEALLOCATE ( ILS(I)%FNC )
  END DO
!
  DO I = 1, NXSC
    IF ( XSC(I)%NXT .GT. 0 ) DEALLOCATE ( XSC(I)%IOF, XSC(I)%ITRI, &
      XSC(I)%NPT, XSC(I)%PRE, XSC(I)%TEM, XSC(I)%ABS, XSC(I)%DWN, XSC(I)%WN1 )
  END DO
!
END SUBROUTINE RFMDAL
END MODULE RFMDAL_SUB
