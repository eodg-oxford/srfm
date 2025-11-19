MODULE SPCFOV_SUB
CONTAINS
SUBROUTINE SPCFOV 
!
! VERSION
!   19FEB20 AD Bug#25: Set TRA=0 before adding convolutions. Checked
!   30MAY18 AD Bug#7: rewritten and simplified
!   02MAY18 AD Bug#1.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Convolve spectra with FOV
!   Called by RFMSPC for each spectral range if FOV flag enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FOVCOM_DAT ! Field of View data
    USE FULCOM_DAT ! Full grid data    
    USE JACCOM_DAT, ONLY: NJAC, ITNJAC  ! No.Jacobians, Jacobians tan.paths
    USE TANCOM_DAT, ONLY: NTAN    ! No. Tot, nominal tangent paths
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IFOV   ! Counter for FOV tabulation points
    INTEGER(I4) :: IJAC   ! Counter for Jacobian elements
    INTEGER(I4) :: ITAN   ! Counter for output tangent paths
    INTEGER(I4) :: JTAN   ! Index of tan.pth containing FOV contribution
    INTEGER(I4) :: KTAN   ! Counter for final tangent ray path
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO ITAN = 1, NTAN
    TRAFUL(:,ITAN) = 0.0D0
    DO IFOV = 1, NFOV
      JTAN = ITNFOV(ITAN,IFOV) 
      RADFUL(:,ITAN) = RADFUL(:,ITAN) + FOV(IFOV)%FNC * RADFUL(:,JTAN) 
      TRAFUL(:,ITAN) = TRAFUL(:,ITAN) + FOV(IFOV)%FNC * TRAFUL(:,JTAN) 
    END DO
  END DO
!
  DO IJAC = 1, NJAC
    DO ITAN = 1, NTAN
      KTAN = ITNJAC(ITAN,IJAC)
      IF ( KTAN .EQ. 0 ) CYCLE     ! zero Jacobian after FOV-conv
      TRAFUL(:,KTAN) = 0.0D0
      DO IFOV = 1, NFOV
        JTAN = ITNFOV(ITAN,IFOV)   ! Unperturbed FOV path
        JTAN = ITNJAC(JTAN,IJAC)   ! Jacobian path
        IF ( JTAN .EQ. 0 ) CYCLE   ! zero Jacobian from this JAC,FOV path
        RADFUL(:,KTAN) = RADFUL(:,KTAN) + FOV(IFOV)%FNC * RADFUL(:,JTAN) 
        TRAFUL(:,KTAN) = TRAFUL(:,KTAN) + FOV(IFOV)%FNC * TRAFUL(:,JTAN) 
      END DO
    END DO
  END DO
!
END SUBROUTINE SPCFOV
END MODULE SPCFOV_SUB

