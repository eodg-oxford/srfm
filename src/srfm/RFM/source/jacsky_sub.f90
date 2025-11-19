MODULE JACSKY_SUB
CONTAINS
SUBROUTINE JACSKY
!
! VERSION
!   27APR24 AD Checked.
!   05MAR19 AD Original. Checked.
!
! DESCRIPTION
!   Allocate diffuse sky tan paths for Jacobians
!   Called by JACPTH
!   Generally separate Jacobian tan paths are assigned for diffuse sky rads
!   and normal tan paths. Four options depending on established sensitivity
!        Nom  Sky
!    1:  no    no  - eg surface Jacobian for limb path
!    2:  no   yes  - eg atm Jacobian above observer viewing surface
!    3:  yes   no  - eg surface Jacobian for surface intersect.path
!    4:  yes  yes  - eg general atm Jacobian
!   case 1: ITNJAC = 0 for both nom and sky paths - no change required
!   case 2: Need to define additional nom tan path to allow for sky Jac
!   case 3: JTAN=ITNJAC(ITAN,IJAC) will have copied %SKY from orig, so OK
!   case 4: Need to assign TAN(JTAN)%ISK to ptb sky for same Jacobian
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE JACCOM_DAT ! Jacobian data
    USE TANCOM_DAT ! Tangent path data
!
! SUBROUTINES
    USE ADDTAN_SUB ! Add tangent ray path
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IJAC   ! Jacobian element counter
    INTEGER(I4) :: ISKTAN ! Original (unperturbed) Diff.Sky path
    INTEGER(I4) :: ITAN   ! Unperturbed tangent ray indices
    INTEGER(I4) :: JSKTAN ! Diff.Sky path for Jacobian
    INTEGER(I4) :: JTAN   ! Tan path used for Jacobian
!
! EXECUTABLE CODE ------------------------------------------------------------
!
  DO ITAN = 1, NTNJAC
!
! Ignore diff.sky paths or paths which do not depend on diff.sky calc.
    IF ( TAN(ITAN)%SKY .OR. TAN(ITAN)%ISK .EQ. 0 ) CYCLE
!
    ISKTAN = TAN(ITAN)%ISK  ! index of original diff.sky path
    DO IJAC = 1, NJAC
      JSKTAN = ITNJAC(ISKTAN,IJAC) ! Diffuse sky Jacobian path
      IF ( JSKTAN .EQ. 0 ) CYCLE   ! Cases 1&3: Diff.sky insensitive to Jacobian
      JTAN = ITNJAC(ITAN,IJAC)
      IF ( JTAN .EQ. 0 ) THEN      ! Case 2
        CALL ADDTAN ( ITAN, .TRUE. ) 
        ITNJAC(ITAN,IJAC) = MTAN
        TAN(MTAN)%ISK = JSKTAN
        TAN(MTAN)%JDX = JAC(IJAC)%JDX
      ELSE                         ! Case 4
        TAN(JTAN)%ISK = JSKTAN
      END IF
    END DO
  END DO
!
END SUBROUTINE JACSKY
END MODULE JACSKY_SUB

