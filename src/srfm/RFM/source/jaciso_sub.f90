MODULE JACISO_SUB
CONTAINS
SUBROUTINE JACISO ( IDXMOL, IDXISO, FAIL, ERRMSG )
!
! VERSION
!   20MAR23 AD Checked.
!   07NOV17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Check valid parameters for isotopic Jacobian
!   Called by JACTGT each target where isotope is specified.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LFLCOM_DAT ! Look-Up Table files
!
! SUBROUTINES
    USE ADDGAS_SUB ! Add new molecule/isotope to list of absorbers
    USE ADDVMR_SUB ! Add extra VMR profile to ATMCOM
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
!
  IMPLICIT NONE
!
! ARGUMENTS 
    INTEGER(I4),   INTENT(IN)  :: IDXMOL ! HITRAN/RFM ID for molecule
    INTEGER(I4),   INTENT(IN)  :: IDXISO ! Isotope# (1=most abundant, etc)
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL     :: NEWGAS ! T=new isotope, F=isotope already separated
    INTEGER(I4) :: IGAS   ! Default index of molecule in GASCOM
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
  CALL ADDGAS ( IDXMOL, IDXISO, NEWGAS ) ! Check if already exists, else add
!
  IF ( .NOT. NEWGAS ) RETURN  ! Separate gas,profile already set up
!
  IGAS = IDXGAS ( IDXMOL )  
!
! Cannot introduce new isotopes if LUTs have already been loaded
  IF ( NLFL .GT. 0 ) THEN  ! Some LUTs already loaded
    IF ( ANY ( LFL%IGS .EQ. IGAS ) ) THEN
      ERRMSG = 'F-JACISO: *JAC section should precede *LUT' & 
               // ' section if defining isotopic Jac.'
      FAIL = .TRUE.
      RETURN
    END IF
  END IF
!            
  CALL ADDVMR ( IGAS ) 
!
END SUBROUTINE JACISO
END MODULE JACISO_SUB

