MODULE ADDQFN_SUB
CONTAINS
SUBROUTINE ADDQFN ( IDXMOL, IDXISO, IQFN )
!
! VERSION
!   17SEP24 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Add vibrational partition function
!   Called by NTEFIL
!   Returns index of existing/new vib partition function
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE NTECOM_DAT ! Non-LTE data
    USE QFNCOM_DAT ! Non-LTE Vib.Partition Fn data
    USE ATMCOM_DAT, ONLY: NATM, QFNATM ! Vib.Partition Fn profiles
!
! SUBROUTINES
    USE IDXQFN_FNC ! Index in QFNCOM of molec,isotope
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN)  :: IDXMOL ! HITRAN/RFM index of molecule
    INTEGER(I4), INTENT(IN)  :: IDXISO ! HITRAN isotope#
    INTEGER(I4), INTENT(OUT) :: IQFN   ! Index of new/existing Vib.Part.Fn
!
! LOCAL VARIABLES
    INTEGER(I4)               :: INTE        ! Counter for non-LTE levels
    REAL(R4),     ALLOCATABLE :: R2DSAV(:,:) ! Saved 2D real array
    TYPE(QFNTYP), ALLOCATABLE :: QFNSAV(:)   ! Saved QFN
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IQFN = IDXQFN ( IDXMOL, IDXISO ) 
  IF ( IQFN .NE. 0 ) RETURN          ! already loaded
!
  IF ( NQFN .GT. 0 ) THEN
    CALL MOVE_ALLOC ( QFN, QFNSAV )
    CALL MOVE_ALLOC ( QFNATM, R2DSAV )
  END IF
!
  NQFN = NQFN + 1
  ALLOCATE ( QFN(NQFN) )
  ALLOCATE ( QFNATM(NATM,NQFN) )
!
  IF ( NQFN .GT. 1 ) THEN
    QFN(1:NQFN-1) = QFNSAV
    QFNATM(:,1:NQFN-1) = R2DSAV
  END IF
!
  QFN(NQFN)%IDM = IDXMOL
  QFN(NQFN)%IDI = IDXISO
!
! Assign this Vib.Partition Fn to any existing non-LTE levels
  DO INTE = 1, NNTE
    IF ( NTE(INTE)%IDM .EQ. IDXMOL .AND. NTE(INTE)%IDI .EQ. IDXISO ) &
      NTE(INTE)%IQF = NQFN
  END DO
!
  IQFN = NQFN
!
END SUBROUTINE ADDQFN
END MODULE ADDQFN_SUB
