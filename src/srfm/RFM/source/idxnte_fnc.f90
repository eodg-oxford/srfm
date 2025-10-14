MODULE IDXNTE_FNC
CONTAINS
INTEGER(I4) PURE FUNCTION IDXNTE ( IDXVIB )
!
! VERSION
!   14SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Index in NTECOM of encoded Vib.index
!   General purpose module.
!   Returns 0 if CODE not assigned as a retrieval parameter
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL VARIABLES
    USE NTECOM_DAT ! Non-LTE data
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IDXVIB ! Code of retrieved parameter
!
! LOCAL VARIABLES
    INTEGER(I4) :: INTE ! Counter
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  DO INTE = 1, NNTE
    IF ( NTE(INTE)%IDX .EQ. IDXVIB ) THEN
      IDXNTE = INTE
      RETURN
    END IF
  END DO
!
  IDXNTE = 0
!
END FUNCTION IDXNTE
END MODULE IDXNTE_FNC


