MODULE INISVD_SUB
CONTAINS
SUBROUTINE INISVD ( ISPC, FAIL, ERRMSG )
!
! VERSION
!   07JUN24 AD Checked.
!   20DEC17 AD Adapted from INILUT. Checked.
!
! DESCRIPTION
!   Initialise SVD data for each new spectral range
!   Called by SPCINI at the start of each spectral range.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE SFLCOM_DAT ! SVD-compressed LUT files
    USE SVDCOM_DAT ! SVD-compressed LUT data
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE REASVD_SUB ! Load SVD-compressed LUT data from file
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
    INTEGER(I4) :: ISFL   ! Index of SVD-LUT files
    INTEGER(I4) :: ISVD   ! Counter for SVD-LUT files
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
!
! First calculate how many SVD and LUT files are to be used for spectral range
  NSVD = COUNT ( IDXSFL(:,ISPC) .GT. 0 )
  IF ( NSVD .EQ. 0 ) RETURN ! No SVD-LUTs for this spec.range
!
  ALLOCATE ( SVD(NSVD) ) 
!
  ISVD = 0
  DO IGAS = 1, NGAS
    ISFL = IDXSFL(IGAS,ISPC) 
    IF ( ISFL .EQ. 0 ) CYCLE
    GAS(IGAS)%HIT = .FALSE.    ! set flag=F since replaced by LUT
    GAS(IGAS)%XSC = .FALSE.
    ISVD = ISVD + 1
    CALL REASVD ( LUNTMP, SFL(ISFL)%NAM, ISVD, SFL(ISFL)%NL, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
  END DO
!
END SUBROUTINE INISVD
END MODULE INISVD_SUB
