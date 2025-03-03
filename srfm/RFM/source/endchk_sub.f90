MODULE ENDCHK_SUB
CONTAINS
SUBROUTINE ENDCHK ( LUNDRV, KEY, FAIL, ERRMSG )
!
! VERSION
!   19DEC23 AD Checked.
!   01APR17 AD F90 version. Tested

! DESCRIPTION
!   Check end of Driver Table section has been reached
!   General purpose module.
!   Returns with FAIL=TRUE if any fields found before end of current section.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE NXTFLD_SUB ! Load next field from section of driver file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for driver table
    CHARACTER(4),  INTENT(IN)  :: KEY    ! Driver Table section
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: LENGTH ! No. of characters in FIELD
    CHARACTER(30) :: FIELD  ! Extra field read from driver table
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
  IF ( LENGTH .NE. 0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-ENDCHK: Unexpected extra field in ' // KEY // &
             ' section: ' // FIELD
  END IF
!
END SUBROUTINE ENDCHK
END MODULE ENDCHK_SUB

