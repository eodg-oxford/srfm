MODULE INILBL_SUB
CONTAINS
SUBROUTINE INILBL
!
! VERSION
!   20AUG24 AD Checked.
!   11AUG23 AD Allow for multiple HITRAN files
!   22MAY22 AD Change HFLCOM variable name USEIDG to USEIDM
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION  
!   Initialise line-by-line calc segments
!   Called by SPCINI at start of each spectral range
!   Only needs to be called once unless LUTs are used, which are an alternative
!   source of absorption data to line-by-line calculations.
!
! VARIABLE KINDS
    USE KIND_DAT  
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE GASCOM_DAT ! Molecule and isotope data
    USE HFLCOM_DAT ! HITRAN file data
    USE WIDCOM_DAT, ONLY: NLBL, IDXLBL ! line-by-line calc paths
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: ICLC   ! Index of calculated path segment
    INTEGER(I4) :: IDXMOL ! HITRAN/RFM index of (line) molecule
    INTEGER(I4) :: IGAS   ! Counter for absorbers
    INTEGER(I4) :: IHFL   ! Index of HITRAN data file
    INTEGER(I4) :: ILBL   ! Index of line-by-line calc.path segment
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Clear all use-molecule flags for HITRAN files (info saved in IFLHFL)
  DO IHFL = 1, NHFL
    HFL(IHFL)%IDM(:) = .FALSE.
  END DO
!
  CLC%LBL = .FALSE.
  DO IGAS = 1, NGAS
    IF ( GAS(IGAS)%HIT ) THEN
      IDXMOL = GAS(IGAS)%IDM
      IHFL = IFLIDM(IDXMOL)
      HFL(IHFL)%IDM(IDXMOL) = .TRUE.
      WHERE ( CLC%IGS .EQ. IGAS ) CLC%LBL = .TRUE.
    ELSE IF ( GAS(IGAS)%CTM ) THEN
! Also need to flag continuum-only gases for LBL calc since continuum is
! interpolated to widemesh and them subsequently interpolated to finemesh in
! the same way as line wing contributions
      WHERE ( CLC%IGS .EQ. IGAS ) CLC%LBL = .TRUE.
    END IF
  END DO
!
  NLBL = COUNT ( CLC%LBL ) 
  IF ( ALLOCATED ( IDXLBL ) ) DEALLOCATE ( IDXLBL ) 
  ALLOCATE ( IDXLBL(NLBL) ) 
  ILBL = 0
  DO ICLC = 1, NCLC
    IF ( CLC(ICLC)%LBL ) THEN
      ILBL = ILBL + 1
      IDXLBL(ILBL) = ICLC
    END IF
  END DO
!  
END SUBROUTINE INILBL
END MODULE INILBL_SUB
