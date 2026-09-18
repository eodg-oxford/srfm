MODULE IDXVIB_FNC
CONTAINS
INTEGER(I4) FUNCTION IDXVIB ( IDM, VIBSTR )
!
! VERSION
!   24APR26 AD Original.
!
! DESCRIPTION
!   Index in VIBCOM of vibrational level index
!   Called by RECPAR for each line if Vib level data loaded
!   Returns 0 if no index assigned to HITRAN vib level ID string
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL VARIABLES
    USE VIBCOM_DAT ! Vibrational level indices
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN) :: IDM    ! HITRAN molecule index
    CHARACTER(15), INTENT(IN) :: VIBSTR ! HITRAN vib.lev string
!
! LOCAL VARIABLES
    INTEGER(I4) :: ILEV ! Counter for levels within molec.class
    INTEGER(I4) :: IMOL ! Counter for stored molecules
    INTEGER(I4) :: IVIB ! Index of molec class in VIB array
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  DO IMOL = 1, NMOL
    IF ( IDXMOL(IMOL) .EQ. IDM ) THEN
      IVIB = IVBMOL(IMOL)
      DO ILEV = 1, VIB(IVIB)%NLV
        IF ( VIB(IVIB)%STR(ILEV) .EQ. VIBSTR ) THEN
          IDXVIB = VIB(IVIB)%IDX(ILEV)
          RETURN
        END IF
      END DO
      IDXVIB = 0
      RETURN
    END IF
  END DO
!
! All IDM values should be listed in IDXMOL
  STOP 'F-IDXVIB: Logical error'
!
END FUNCTION IDXVIB
END MODULE IDXVIB_FNC


