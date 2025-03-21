MODULE REAHIT_SUB
CONTAINS
SUBROUTINE REAHIT ( EOF, FAIL, ERRMSG ) 
!
! VERSION
!   23AUG24 AD Checked.
!   11AUG23 AD Use single HITREC to read data
!   23MAY23 AD Different read subroutines depending on HITRAN data type
!   29JAN20 AD Use HITREC to read from all types of HITRAN input files. Checked
!   13JUN17 AD Add REAPAR
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Read record from HITRAN line data file
!   Called by REACYC, SPCWID.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE HITCOM_DAT ! HITRAN line data
    USE REJCOM_DAT ! Minimum line strength limits
!
! SUBROUTINES
    USE HITREC_SUB ! Read record from HITRAN line data file
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
    USE SETNTE_SUB ! Set IUSNTE,ILSNTE in HITCOM
    USE USEQAL_FNC ! Set TRUE if listed in line molecule qualifiers
    USE VALISO_FNC ! Check recognised isotope
!    
  IMPLICIT NONE
!
! ARGUMENTS      
    LOGICAL,       INTENT(OUT) :: EOF    ! Set TRUE if end-of-file reached
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)  :: IDXMOL ! HITRAN index of molecule
    INTEGER(I4)  :: IGAS   ! Index of molecule in GASCOM
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
!
! Find next useful record in file. 
  DO 
    CALL HITREC ( HIT, EOF, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
    IF ( EOF ) RETURN
!
    IDXMOL = HIT%IDM
! Read routines should only return required molecules
    IF ( IGSMOL(IDXMOL) .EQ. 0 ) STOP 'F-REAHIT: Logical error'
    IF ( .NOT. VALISO ( IGSMOL(IDXMOL), HIT%IDI ) ) CYCLE
!
    IGAS = IDXGAS ( HIT%IDM, HIT%IDI )
    IF ( USEREJ ) THEN
      IF ( HIT%STR .LT. STRREJ(IGAS) ) CYCLE
    END IF
    IF ( GAS(IGAS)%QAL ) THEN
      IF ( .NOT. USEQAL ( HIT%IDM, HIT%IDI, HIT%ILS, HIT%IUS ) ) CYCLE
    END IF
    HIT%IGS = IGAS
    HIT%WGT = GAS(IGAS)%WGT(HIT%IDI)
    IF ( GAS(IGAS)%NTE ) CALL SETNTE
!
    EXIT    ! next required line loaded
  END DO
!
END SUBROUTINE REAHIT
END MODULE REAHIT_SUB
