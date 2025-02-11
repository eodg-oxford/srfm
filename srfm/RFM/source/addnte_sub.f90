MODULE ADDNTE_SUB
CONTAINS
SUBROUTINE ADDNTE ( IDXVIB, ENERGY, NEWNTE )
!
! VERSION
!   23OCT23 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Add vibrational temperature to NTE data
!   Called by ATMPRF and NTEFIL
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE NTECOM_DAT ! Non-LTE data
    USE ATMCOM_DAT, ONLY: NTEVMR ! T=use non-LTE for molecule
!
! SUBROUTINES
    USE ADDVIB_SUB ! Add extra Vib Temp profile to ATMCOM
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
    USE IDXNTE_FNC ! Index in NTECOM of encoded Vib.index
    USE IDXQFN_FNC ! Index in QFNCOM of molec,isotope
    USE MOLIDX_SUB ! Give molecule name for HITRAN/RFM index, or vice-versa
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),        INTENT(IN)  :: IDXVIB ! Encoded Vib Temp data
    REAL(R4), OPTIONAL, INTENT(IN)  :: ENERGY ! Transition energy [cm-1]
    LOGICAL,  OPTIONAL, INTENT(OUT) :: NEWNTE ! T=new NTE data, F=already listed
!
! LOCAL CONSTANTS
! local list of energy levels for required VT profiles if not supplied as ENERGY
! Get info from .nte files if more are required
    INTEGER(I4), PARAMETER :: NVIB = 2  ! No. of Vib.Temps 
    INTEGER(I4), PARAMETER :: LSTVIB(NVIB) = (/ 5000002, 5001002 /) ! CO(0:1)(2)
    REAL(R4),    PARAMETER :: ENGVIB(NVIB) = (/ 2143.2711, 2143.2711 /)
!
! LOCAL VARIABLES
    INTEGER(I4)       :: IDXISO ! RFM/HITRAN Isotopic ID
    INTEGER(I4)       :: IDXLEV ! Vibrational level# (=global quantum no.)
    INTEGER(I4)       :: IDXMOL ! RFM/HITRAN ID for VT molecule 
    INTEGER(I4)       :: IGAS   ! Index of mol/iso in GASCOM
    INTEGER(I4)       :: INTE   ! Index in NTE 
    INTEGER(I4)       :: IVIB   ! Counter for locally stored Energy levels info
    CHARACTER(LENGAS) :: MOL    ! Molecule name 
    TYPE(NTETYP), ALLOCATABLE :: NTESAV(:) ! Copy of NTE during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  INTE = IDXNTE ( IDXVIB ) 
  IF ( INTE .GT. 0 ) THEN                             ! already listed
    IF ( PRESENT ( ENERGY ) ) NTE(INTE)%ENG = ENERGY  ! just update Energy value
    IF ( PRESENT ( NEWNTE ) ) NEWNTE = .FALSE.
    RETURN
  END IF
!
  IDXMOL = IDXVIB/1000000     ! Extract Molecule ID
  IDXISO = MOD ( IDXVIB, 1000000 ) / 1000
  IDXLEV = MOD ( IDXVIB, 1000 )
!
  IGAS = IDXGAS ( IDXMOL, IDXISO )
  GAS(IGAS)%NTE = .TRUE.
  NTEVMR(IGAS)  = .TRUE.
!
  IF ( ALLOCATED ( NTE ) ) CALL MOVE_ALLOC ( NTE, NTESAV )
  NNTE = NNTE + 1
  ALLOCATE ( NTE(NNTE) )
  IF ( ALLOCATED ( NTESAV ) ) NTE(1:NNTE-1) = NTESAV
!
  NTE(NNTE)%IDX = IDXVIB
  NTE(NNTE)%IDM = IDXMOL
  NTE(NNTE)%IDI = IDXISO
  NTE(NNTE)%IDL = IDXLEV
  NTE(NNTE)%IQF = IDXQFN ( IDXMOL, IDXISO ) 
  CALL MOLIDX ( IDXMOL, MOL )
  NTE(NNTE)%COD = TRIM ( MOL ) // '(' // TRIM ( C11INT(IDXISO) ) &
                  // ')(' // TRIM ( C11INT(IDXLEV) ) // ')'
!
! Check that the energy level for this VT is stored here 
  IF ( PRESENT ( ENERGY ) ) THEN
    NTE(NNTE)%ENG = ENERGY
  ELSE
    NTE(NNTE)%ENG = 0.0
    DO IVIB = 1, NVIB
      IF ( IDXVIB .EQ. LSTVIB(IVIB) ) THEN
        NTE(NNTE)%ENG = ENGVIB(IVIB)
        EXIT
      ENDIF
    END DO
  END IF
!
  CALL ADDVIB  ! Create extra Vib.Temp profile
!
  IF ( PRESENT ( NEWNTE ) ) NEWNTE = .TRUE.
!
END SUBROUTINE ADDNTE
END MODULE ADDNTE_SUB
