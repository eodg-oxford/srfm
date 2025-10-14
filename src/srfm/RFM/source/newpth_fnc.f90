MODULE NEWPTH_FNC
CONTAINS
LOGICAL FUNCTION NEWPTH ( PTB, ORG, JDX )
!
! VERSION
!   03DEC21 AD Only check criterion which matches Jacobian type. Checked.
!   06OCT21 AD Original.
!
! DESCRIPTION
!   T=significantly different Jacobian path segment
!   Called by JACPTH
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE PTBCON_DAT ! Jacobian perturbation sizes
    USE JACCOM_DAT, ONLY: JDXTEM, JDXPRE  ! JDX values for Tem, Pre Jacob.
    USE PTHCOM_DAT, ONLY: PTHTYP ! Path data type
!
  IMPLICIT NONE
!
! ARGUMENTS
    TYPE(PTHTYP), INTENT(IN) :: PTB ! Perturbed path
    TYPE(PTHTYP), INTENT(IN) :: ORG ! Unperturbed path
    INTEGER(I4),  INTENT(IN) :: JDX ! Index of target species
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: SIGVMR = 0.01  ! Fraction of ptb requiring new path
    REAL(R4), PARAMETER :: SIGTEM = 0.01  ! Fraction of ptb requiring new path
    REAL(R4), PARAMETER :: SIGPRE = 0.01  ! Fraction of ptb requiring new path
!
! EXECUTABLE CODE ------------------------------------------------------------
!
  SELECT CASE ( JDX ) 
  CASE ( JDXTEM )  
    NEWPTH = ABS ( PTB%TEM - ORG%TEM ) .GT. SIGTEM * PTBTEM 
  CASE ( JDXPRE ) 
    NEWPTH = ABS ( PTB%PRE/ORG%PRE - 1.0 ) .GT. SIGPRE * PTBPRE
  CASE DEFAULT 
    NEWPTH = ABS ( PTB%AMT/ORG%AMT - 1.0 ) .GT. SIGVMR * PTBVMR 
  END SELECT
!
END FUNCTION NEWPTH
END MODULE NEWPTH_FNC
