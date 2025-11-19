MODULE CIAMOL_SUB
CONTAINS
SUBROUTINE CIAMOL ( MOLSTR, IDX1, IDX2, FAIL, ERRMSG ) 
!
! VERSION
!   30JUL24 AD Checked.
!   24JUN19 AD Change to use IDXMOL to determine indices. Checked.
!   01MAY17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   Identify molecules from .cia file header record
!   Called by CIAFIL.
!   The .cia header record contains a C*20 right-justified string containing
!   the pair of molecules, eg: '     O2-N2'
!   This will return RFM indices of both molecules, or 0 if not a recognised
!   molecule or not a molecule name added in *GAS section (eg 'ar') 
!   The indices will be in order IDX1 .LE. IDX2, but with 'air' always as IDX2
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE IDXCON_DAT, ONLY: IDXAIR ! RFM 'molecular' index for air
!
! SUBROUTINES
    USE MOLIDX_SUB ! Give molecule name for HITRAN/RFM index, or vice-versa
!
  IMPLICIT NONE
!
! ARGUMENTS      
    CHARACTER(*),  INTENT(IN)  :: MOLSTR ! String containing molecules
    INTEGER(I4),   INTENT(OUT) :: IDX1   ! HITRAN/RFM index of first molecule
    INTEGER(I4),   INTENT(OUT) :: IDX2   ! HITRAN/RFM index of second molecule
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)            :: I      ! Location of '-' in MOLSTR
    INTEGER(I4)            :: IDUMMY ! Temporary storage for swapping 
    CHARACTER(LEN(MOLSTR)) :: MOL    ! Molecule name extracted from MOLSTR
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Locate '-' character in MOLSTR
  I = INDEX ( MOLSTR, '-' ) 
  IF ( I .LE. 1 .OR. I .EQ. LEN ( MOLSTR ) ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-CIAMOL: Unable to read molecules from .cia header:' // MOLSTR
    RETURN
  END IF
!
  IDX1 = 0
  MOL = MOLSTR(1:I-1)
  CALL MOLIDX ( IDX1, MOL ) 
!
  IDX2 = 0
  MOL = MOLSTR(I+1:) 
  CALL MOLIDX ( IDX2, MOL ) 
!
! Ensure IDX1 .le. IDX2
  IF ( IDX2 .LT. IDX1 .OR. IDX1 .EQ. IDXAIR ) THEN
    IDUMMY = IDX1
    IDX1 = IDX2
    IDX2 = IDUMMY
  END IF
!
  FAIL = .FALSE.
! 
END SUBROUTINE CIAMOL
END MODULE CIAMOL_SUB
