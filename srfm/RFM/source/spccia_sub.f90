MODULE SPCCIA_SUB
CONTAINS
SUBROUTINE SPCCIA
!
! VERSION
!   31AUG24 AD Checked.
!   01MAY17 AD F90 conversion of part of rfmcia.for. Checked.
! 
! DESCRIPTION    
!   Calculate the Collision-induced absorption 
!   Called by RFMSPC for each wide mesh interval.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CIPCOM_DAT ! CIA path data
    USE FINCOM_DAT ! Finemesh data
!
! SUBROUTINES
    USE CIAINT_FNC ! Interpolate CIA spectrum
!                  
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: ICLC ! Index of calculated path
    INTEGER(I4) :: ICIP ! Counter for CIA paths
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO ICIP = 1, NCIP
    ICLC = CIP(ICIP)%ICL
    ABSFIN(:,ICLC) = ABSFIN(:,ICLC) +  CIP(ICIP)%AM2 * &
       CIAINT ( CIP(ICIP)%ID1, CIP(ICIP)%ID2, CIP(ICIP)%TEM, NFIN, WNOFIN ) 
  END DO
!
END SUBROUTINE SPCCIA
END MODULE SPCCIA_SUB
