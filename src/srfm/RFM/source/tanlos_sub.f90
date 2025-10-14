MODULE TANLOS_SUB
CONTAINS
SUBROUTINE TANLOS 
!
! VERSION
!   31MAR24 AD Checked.
!   04FEB19 AD Original. Checked.
!
! DESCRIPTION
!   Set up LOS Jacobians
!   Called once by DRVTAN if LOS flag enabled.
!   Just adds 'Jacobian' element for 'los' perturbation
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE JACCOM_DAT ! Jacobian data
    USE FLGCOM_DAT, ONLY: JACFLG ! T = calculate Jacobian spectra
    USE TANCOM_DAT, ONLY: NTAN   ! No. of nominal tangent paths
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    TYPE(JACTYP), ALLOCATABLE :: JACSAV(:) ! Saved JAC during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( ALLOCATED ( JAC ) ) CALL MOVE_ALLOC ( JAC, JACSAV ) 
  NJAC = NJAC + 1
  ALLOCATE ( JAC(NJAC) ) 
  IF ( ALLOCATED ( JACSAV ) ) JAC(1:NJAC-1) = JACSAV
  JAC(NJAC)%COD = 'los'
  JAC(NJAC)%COL = .TRUE.
!
! Ensure ITNJAC is defined if no further Jacobians added
  IF ( .NOT. JACFLG ) THEN 
    ALLOCATE ( ITNJAC(NTAN,NJAC) ) 
    ITNJAC = 0
  END IF
!
END SUBROUTINE TANLOS
END MODULE TANLOS_SUB
