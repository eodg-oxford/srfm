MODULE JACFOV_SUB
CONTAINS
SUBROUTINE JACFOV
!
! VERSION
!   24NOV23 AD Checked.
!   30MAY18 AD Bug#7 rewritten and simplified. Checked.
!   02MAY18 AD original.
!
! DESCRIPTION
!   Assign tan paths for Jacobians after FOV convolution
!   Called by JACPTH.
!   Ensure that ITNJAC(1:NTAN,1:NJAC) is non-zero for any convolved paths
!   which have finite Jacobian. 
!   This adds further elements to TAN(:), increasing MTAN, although the new 
!   tangent paths are flagged for no calculation.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FOVCOM_DAT, ONLY: ITNFOV, NFOV ! Field of View data
    USE JACCOM_DAT, ONLY: ITNJAC, NJAC ! Jacobian data
    USE TANCOM_DAT, ONLY: MTAN, NTAN   ! Tot, Nom. tangent paths
!
! SUBROUTINES
    USE ADDTAN_SUB ! Add tangent ray path
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IFOV   ! FOV counter
    INTEGER(I4)   :: IJAC   ! Jacobian element counter
    INTEGER(I4)   :: ITAN   ! Unperturbed tangent ray indices
    INTEGER(I4)   :: JTAN   ! Tan path used for Jacobian
!
! EXECUTABLE CODE ------------------------------------------------------------
!
  DO ITAN = 1, NTAN
    DO IJAC = 1, NJAC
      DO IFOV = 1, NFOV
        JTAN = ITNFOV(ITAN,IFOV)
        IF ( ITNJAC(JTAN,IJAC) .NE. 0 ) THEN   ! need to assign a tan path
          CALL ADDTAN ( ITAN, .FALSE. )
          ITNJAC(ITAN,IJAC) = MTAN
          EXIT
        END IF
      END DO 
    END DO
  END DO
!
END SUBROUTINE JACFOV
END MODULE JACFOV_SUB

