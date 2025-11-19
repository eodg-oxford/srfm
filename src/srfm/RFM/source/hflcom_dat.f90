MODULE HFLCOM_DAT
!
! VERSION
!   10AUG24 AD Checked.
!   11AUG23 AD Rewritten again
!   23MAY23 AD Extensive modifications
!   29JAN20 AD Add BASHFL, BINHFL, NFPHFL, LSTFWD, LSTLIN. Checked.
!   13JUN17 AD Add PARHFL
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   HITRAN line file general data.
!   Initially loaded by HITOPN, updated by various routines, and released via
!   HFLCOM_RESET to close any open LUNs when re-running within one process.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE HDBCOM_DAT ! HDB data structure
    USE HITCOM_DAT, ONLY: HITTYP ! HITRAN line data structure
!
  IMPLICIT NONE
  SAVE
  PUBLIC :: HFLCOM_RESET
!
! GLOBAL CONSTANTS
    INTEGER(I4), PARAMETER :: MAXIDM = 98 ! Max HITRAN index of reqd molecules
!
  TYPE :: HFLTYP
    INTEGER(I4)  :: LUN         ! Logical Unit Number of file
    CHARACTER(3) :: TYP         ! Type of HITRAN file, 'BIN','PAR' or 'HDB'
    TYPE(HITTYP) :: HIT         ! HITRAN line data structure
    TYPE(HDBTYP) :: HDB         ! HITRAN database file auxiliary data
    LOGICAL      :: IDM(MAXIDM) ! T=use molec# from file
    INTEGER(I4)  :: IFP(MAXIDM) ! Forward pointers for binary files
  END TYPE HFLTYP

! GLOBAL VARIABLES
    INTEGER(I4)  :: NHFL = 0            ! No. different HITRAN files being used
    REAL(R8)     :: WNOMAX = 0.0D0      ! Max. wavenumber to read to
    INTEGER(I4)  :: IFLIDM(MAXIDM)      ! Index of file used for each molec
    TYPE(HFLTYP), ALLOCATABLE :: HFL(:) ! Structure for each HITRAN file
!
CONTAINS

  SUBROUTINE HFLCOM_RESET()
    INTEGER(I4) :: I
    LOGICAL :: IS_OPEN

    IF ( ALLOCATED ( HFL ) ) THEN
      DO I = 1, SIZE ( HFL )
        INQUIRE ( UNIT=HFL(I)%LUN, OPENED=IS_OPEN )
        IF ( IS_OPEN ) CLOSE ( UNIT=HFL(I)%LUN )
      END DO
      DEALLOCATE ( HFL )
    END IF
    NHFL   = 0
    WNOMAX = 0.0_R8
    IFLIDM = 0
  END SUBROUTINE HFLCOM_RESET

END MODULE HFLCOM_DAT
