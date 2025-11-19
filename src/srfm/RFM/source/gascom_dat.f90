MODULE GASCOM_DAT
!
! VERSION
!   03OCT24 AD Checked.
!   18APR22 AD Add IAIGAS. Checked.
!   24JUN19 AD Add %CIA. Checked.
!   01JUN17 AD F90 version. Checked.
!
! DESCRIPTION
!   Molecule and isotope data.
!   Loaded by ADDGAS and cleared by GASCOM_RESET to support repeated in-memory
!   runs. GAS%QAL set TRUE by GASQAL. Pointers GAS%ISO, GAS%WGT are deallocated
!   by RFMDAL and the reset helper.
!
!   MAXMOL should be as least as large as the highest index assigned in 
!   molidx_sub.f90 plus extra to allow for any user-defined x/s molecules.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
  SAVE
  PUBLIC :: GASCOM_RESET
!
! GLOBAL CONSTANTS
    INTEGER(I4), PARAMETER :: LENGAS = 7    ! Maximum length of molecule name
    INTEGER(I4), PARAMETER :: MAXHLN = 98   ! Maximum index for a line molecule
    INTEGER(I4), PARAMETER :: MAXISO = 11   ! Maximum no.isotopologues
    INTEGER(I4), PARAMETER :: MAXMOL = 200  ! Maximum recognised molec.ID
!
  TYPE :: GASTYP
    LOGICAL      :: CIA    ! T=use collision-induced absorption for this gas
    LOGICAL      :: CTM    ! T=use continuum for this gas
    LOGICAL      :: HIT    ! T= line data found in HITRAN file
    LOGICAL      :: NTE    ! T=use non-LTE for molecule
    LOGICAL      :: QAL    ! T=use band/isotope selection qualifiers 
    LOGICAL      :: XSC    ! T=.xsc file identified for molecule
    INTEGER(I4)  :: IDI    ! HITRAN isotope index (or 0)
    INTEGER(I4)  :: IDM    ! HITRAN index for each species 
    INTEGER(I4)  :: NIS    ! No.isotopes for species
    INTEGER(I4)  :: SHP    ! Lineshape to be used for gas
    CHARACTER(LENGAS)    :: COD    ! Character codes for species
    INTEGER(I4), POINTER :: ISO(:) ! [0:NISO] IDXGAS value for each isotope
    REAL(R4),    POINTER :: WGT(:) ! [NISO] Molecular wts [Atomic units] 
  END TYPE GASTYP
!
! GLOBAL VARIABLES
    TYPE(GASTYP), ALLOCATABLE :: GAS(:)
!
    LOGICAL     :: SUBH2O = .FALSE.         ! T=subtract H2O abs at 25cm-1
    LOGICAL     :: ISOMOL(MAXMOL) = .FALSE. ! T=molecule split by isotopes
    INTEGER(I4) :: IAIGAS = 0               ! Index of air
    INTEGER(I4) :: IAXGAS = 0               ! Index of aerosol 
    INTEGER(I4) :: IGSMOL(MAXMOL) = 0       ! IGAS for HITRAN index (or 0).
    INTEGER(I4) :: NGAS = 0                 ! No.species to be used
!
CONTAINS

  SUBROUTINE GASCOM_RESET()
    INTEGER(I4) :: I

    IF ( ALLOCATED ( GAS ) ) THEN
      DO I = 1, SIZE ( GAS )
        IF ( ASSOCIATED ( GAS(I)%ISO ) ) DEALLOCATE ( GAS(I)%ISO )
        IF ( ASSOCIATED ( GAS(I)%WGT ) ) DEALLOCATE ( GAS(I)%WGT )
      END DO
      DEALLOCATE ( GAS )
    END IF

    SUBH2O = .FALSE.
    ISOMOL = .FALSE.
    IAIGAS = 0
    IAXGAS = 0
    IGSMOL = 0
    NGAS   = 0
  END SUBROUTINE GASCOM_RESET

END MODULE GASCOM_DAT
