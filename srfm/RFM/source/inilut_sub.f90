MODULE INILUT_SUB
CONTAINS
SUBROUTINE INILUT ( ISPC, FAIL, ERRMSG )
!
! VERSION
!   08JUN24 AD Checked.
!   20DEC17 AD F90 conversion of spclut.for. Checked.
!
! DESCRIPTION
!   Initialise LUT data for each new spectral range
!   Called by SPCINI at the start of each spectral range.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE LFLCOM_DAT ! Look-Up Table files
    USE LUTCOM_DAT ! TAB LUT data
!
! SUBROUTINES
    USE LUTTAB_SUB ! Initialise LUT-TAB file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: ISPC   ! Current spectral range
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IGAS   ! Counter, index for absorbers
    INTEGER(I4) :: ILFL   ! Index of LUT files
    INTEGER(I4) :: ILUT   ! Counter for TAB-LUT files
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
!
! First calculate how many LUT files are to be used for spectral range
  NLUT = COUNT ( IDXLFL(:,ISPC) .GT. 0 )
  IF ( NLUT .EQ. 0 ) RETURN ! No LUTs for this spec.range
!
  ALLOCATE ( LUT(NLUT) ) 
!
  ILUT = 0
  DO IGAS = 1, NGAS
    ILFL = IDXLFL(IGAS,ISPC) 
    IF ( ILFL .EQ. 0 ) CYCLE
    GAS(IGAS)%HIT = .FALSE.    ! set flag=F since replaced by LUT
    GAS(IGAS)%XSC = .FALSE.
    ILUT = ILUT + 1
    CALL LUTTAB ( LFL(ILFL)%LUN, LFL(ILFL)%NAM, LFL(ILFL)%BIN, &
                  ILUT, LFL(ILFL)%NDP, LFL(ILFL)%NDT, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
  END DO
!
! Normal exit
!
END SUBROUTINE INILUT
END MODULE INILUT_SUB
