MODULE ADDGFL_SUB
CONTAINS
SUBROUTINE ADDGFL ( NAMGRD, TYP, WMN, WMX, IGFL ) 
!
! VERSION
!   26OCT23 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Add irreg. grid file to GFLCOM
!   Called by GRDFIL, SPCFIL.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GFLCOM_DAT ! Irregular grid files
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*), INTENT(IN)  :: NAMGRD ! Name of .grd file
    CHARACTER(3), INTENT(IN)  :: TYP    ! Type of file: 'BIN', 'ASC' or 'GRD'
    REAL(R8),     INTENT(IN)  :: WMN    ! Lower Wno [cm-1] of grid
    REAL(R8),     INTENT(IN)  :: WMX    ! Upper Wno [cm-1] of grid
    INTEGER(I4),  INTENT(OUT) :: IGFL   ! Index assigned to file in GFLCOM 
!
! LOCAL VARIABLES
    TYPE(GFLTYP), ALLOCATABLE :: GFLSAV(:) ! Saved GFL
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( ALLOCATED ( GFL ) ) CALL MOVE_ALLOC ( GFL, GFLSAV ) 
  NGFL = NGFL + 1
  ALLOCATE ( GFL(NGFL) ) 
  IF ( ALLOCATED ( GFLSAV ) ) GFL(1:NGFL-1) = GFLSAV
!
  GFL(NGFL)%NAM = NAMGRD
  GFL(NGFL)%BIN = TYP .EQ. 'BIN'
  GFL(NGFL)%GRD = TYP .EQ. 'GRD'
  GFL(NGFL)%WMN = WMN
  GFL(NGFL)%WMX = WMX
!
  IGFL = NGFL
!
END SUBROUTINE ADDGFL
END MODULE ADDGFL_SUB
