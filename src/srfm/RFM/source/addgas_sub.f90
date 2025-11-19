MODULE ADDGAS_SUB
CONTAINS
SUBROUTINE ADDGAS ( IDXMOL, IDXISO, NEWGAS )
!
! VERSION
!   29SEP24 AD Checked.
!   18APR22 AD Set IAIGAS. Checked.
!   24JUN19 AD Add %CIA. Checked.
!   07NOV17 AD F90 Original. Based on gaschk.for. Checked.
!
! DESCRIPTION
!   Add new molecule/isotope to list of absorbers
!   General purpose module.
!   If already loaded, just sets NEWGAS=FALSE
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE GASCOM_DAT ! Molecule and isotope data
    USE IDXCON_DAT ! RFM/HITRAN indices for specific molecules
    USE SHPCON_DAT ! Line-shape indices
!
! SUBROUTINES
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE ISOLST_SUB ! List isotopomer weights of specific molecule
    USE MOLIDX_SUB ! Give molecule name for HITRAN/RFM index, or vice-versa
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),           INTENT(IN)  :: IDXMOL ! HITRAN/RFM index of molecule
    INTEGER(I4), OPTIONAL, INTENT(IN)  :: IDXISO ! HITRAN isotope#
    LOGICAL,     OPTIONAL, INTENT(OUT) :: NEWGAS ! T=new gas, F=already listed
!
! LOCAL VARIABLES
    INTEGER(I4) :: IGAS ! Index in GAS of previously separated isotopes
    INTEGER(I4) :: IISO ! Counter for isotopes
    INTEGER(I4) :: IMOL ! Local IDXMOL since MOLIDX argument is INOUT 
    INTEGER(I4) :: LGAS ! Initial value of NGAS
    INTEGER(I4) :: MGAS ! Index in GAS of default molecule isotope 
    INTEGER(I4) :: NISO ! No. isotopomers for line molecule
    REAL(R4)    :: WGTISO(MAXISO) = 0      ! List of isotope weights 
    TYPE(GASTYP), ALLOCATABLE :: GASSAV(:) ! temporary store of GAS
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( IDXMOL .LE. 0 .OR. IDXMOL .GT. MAXMOL ) STOP 'F-ADDGAS: Logical error#1'
  LGAS = NGAS
!
  IF ( IGSMOL(IDXMOL) .EQ. 0 ) THEN  ! Molecule not yet used
!
    IF ( ALLOCATED ( GAS ) ) CALL MOVE_ALLOC ( GAS, GASSAV )
    NGAS = NGAS + 1
    ALLOCATE ( GAS(NGAS) )
    IF ( ALLOCATED ( GASSAV ) ) GAS(1:NGAS-1) = GASSAV
!
    IGSMOL(IDXMOL) = NGAS
    GAS(NGAS)%IDM = IDXMOL
    IF ( IDXMOL .EQ. IDXAER ) IAXGAS = NGAS
    IF ( IDXMOL .EQ. IDXAIR ) IAIGAS = NGAS
    GAS(NGAS)%IDI = 0                ! default isotope profile
    GAS(NGAS)%HIT = .FALSE.
    GAS(NGAS)%NTE = .FALSE.
    GAS(NGAS)%QAL = .FALSE.
    GAS(NGAS)%XSC = .FALSE.
    GAS(NGAS)%CIA = .FALSE.
    IMOL = IDXMOL
    CALL MOLIDX ( IMOL, GAS(NGAS)%COD )
    SELECT CASE ( IDXMOL ) 
    CASE ( 1:MAXHLN ) 
      GAS(NGAS)%SHP = SHPVOI
      CALL ISOLST ( IDXMOL, NISO, WGTISO )  
      GAS(NGAS)%NIS = NISO
      ALLOCATE ( GAS(NGAS)%ISO(0:NISO), GAS(NGAS)%WGT(1:NISO) ) 
      GAS(NGAS)%ISO(0:NISO) = NGAS
      GAS(NGAS)%WGT(1:NISO) = WGTISO(1:NISO)
    CASE ( IDXAIR ) 
      GAS(NGAS)%SHP = SHPAIR
      GAS(NGAS)%NIS = 0
    CASE ( IDXAER:MAXMOL ) 
      GAS(NGAS)%SHP = SHPXSC
      GAS(NGAS)%NIS = 0
    CASE DEFAULT 
      STOP 'F-ADDGAS: unrecognised molecule#'
    END SELECT
    GAS(NGAS)%CTM = CTMFLG .AND. &
                  ( IDXMOL .EQ. IDXH2O .OR. IDXMOL .EQ. IDXCO2 .OR. &
                    IDXMOL .EQ. IDXO2  .OR. IDXMOL .EQ. IDXN2 ) 
! For H2O, default is to subtract 25cm-1 baseline when used in combination
! with continuum, or when generating LUTs
    IF ( IDXMOL .EQ. IDXH2O ) SUBH2O = CTMFLG .OR. TABFLG
  END IF
!
  IF ( PRESENT ( IDXISO ) ) THEN
    IF ( IDXISO .LT. 0 .OR. IDXISO .GT. MAXISO ) STOP 'F-ADDGAS: Logical error#2'
    IF ( IDXISO .GT. 0 ) THEN
! 
      MGAS = IGSMOL(IDXMOL)     ! index of default molecule data
      NISO = GAS(MGAS)%NIS 
      IF ( IDXISO .GT. NISO ) STOP 'F-ADDGAS: Logical Error#3'
      IF ( GAS(MGAS)%ISO(IDXISO) .EQ. MGAS ) THEN ! not yet separated
        CALL MOVE_ALLOC ( GAS, GASSAV ) 
        NGAS = NGAS + 1
        ALLOCATE ( GAS(NGAS) )
        GAS(1:NGAS-1) = GASSAV
        GAS(MGAS)%ISO(IDXISO) = NGAS  ! modify iso index in def.molec data
        GAS(NGAS) = GAS(MGAS)         ! copy default molecule data
        ALLOCATE ( GAS(NGAS)%ISO(0:NISO), GAS(NGAS)%WGT(NISO) )
        GAS(NGAS)%WGT = GAS(MGAS)%WGT
        GAS(NGAS)%IDI = IDXISO        ! note current isotope
        GAS(NGAS)%COD = GAS(MGAS)%COD // '(' // TRIM ( C11INT(IDXISO) ) // ')'
        DO IISO = 1, NISO         ! update iso lists any previously sep. isotopes
          IGAS = GAS(MGAS)%ISO(IISO)
          GAS(IGAS)%ISO = GAS(MGAS)%ISO
        END DO
        ISOMOL(IDXMOL) = .TRUE.
      END IF
    END IF
  END IF
!
  IF ( PRESENT ( NEWGAS ) ) NEWGAS = NGAS .GT. LGAS
!
END SUBROUTINE ADDGAS
END MODULE ADDGAS_SUB
