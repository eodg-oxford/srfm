MODULE SPCJAC_SUB
CONTAINS
SUBROUTINE SPCJAC 
!
! VERSION
!   25NOV23 AD Checked.
!   30MAY18 AD Bug#7 Only subtract for TAN%CLC is TRUE. Checked.
!   01MAY17 AD F90 conversion of rfmjac.for. Checked.
!
! DESCRIPTION
!   Calculate Jacobians
!   Called by RFMSPC for each spectral range if JAC flag enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FULCOM_DAT ! Full grid data 
    USE JACCOM_DAT ! Jacobian data
    USE TANCOM_DAT ! Tangent path data
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IJAC   ! Counter for Jacobians
    INTEGER(I4) :: ITAN   ! Counter for output tangent paths
    INTEGER(I4) :: JTAN   ! Index of tan.pth containing FOV contribution
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO ITAN = 1, NTNJAC                 ! Loop over all unperturbed paths
    IF ( .NOT. TAN(ITAN)%CLC ) CYCLE
    DO IJAC = 1, NJAC 
      JTAN = ITNJAC(ITAN,IJAC) 
      IF ( JTAN .EQ. 0 ) CYCLE
      RADFUL(:,JTAN) = RADFUL(:,JTAN) - RADFUL(:,ITAN) 
      OPTFUL(:,JTAN) = OPTFUL(:,JTAN) - OPTFUL(:,ITAN) 
      TRAFUL(:,JTAN) = TRAFUL(:,JTAN) - TRAFUL(:,ITAN) 
    END DO
  END DO
!
END SUBROUTINE SPCJAC
END MODULE SPCJAC_SUB

