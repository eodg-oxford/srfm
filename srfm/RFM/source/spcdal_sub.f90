MODULE SPCDAL_SUB
CONTAINS
SUBROUTINE SPCDAL 
!
! VERSION
!   04JUL24 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Deallocate spectral range level pointers
!   Called by RFMSPC at the completion of each spectral range.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LUTCOM_DAT ! TAB LUT data
    USE SVDCOM_DAT ! SVD-compressed LUT data
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: I ! General counter 
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( ALLOCATED ( LUT ) ) THEN
    DO I = 1, NLUT
      CLOSE ( LUT(I)%LUN ) 
      DEALLOCATE ( LUT(I)%ICL, LUT(I)%IDX, LUT(I)%WGT, LUT(I)%TAB )
    END DO
    DEALLOCATE ( LUT )
  END IF
!
  IF ( ALLOCATED ( SVD ) ) THEN
    DO I = 1, NSVD
      DEALLOCATE ( SVD(I)%U, SVD(I)%K, SVD(I)%ICL, SVD(I)%IDX, SVD(I)%WGT ) 
    END DO
    DEALLOCATE ( SVD ) 
  END IF
!
END SUBROUTINE SPCDAL
END MODULE SPCDAL_SUB
