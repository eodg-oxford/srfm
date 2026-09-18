MODULE REACIA_SUB
CONTAINS
SUBROUTINE REACIA ( LUNCIA, ID1, ID2, NPT, TEM, WNL, WNU, FAIL, ERRMSG )
!
! VERSION
!   16AUG26 AD Checked.
!   18AUG25 AD Simplified: remove CIA%IGG, IGGCIA, IDDCIA
!   30AUG24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   Read data from .cia file
!   Called by CIAFIL.
!   A HITRAN .cia file is subdivided into tabulations of wavenumber and 
!   collisions-induced absorption for a particular pair of molecules,
!   spectral range and temperature.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CIACOM_DAT ! Collision-induced absorption data
    USE PHYCON_DAT, ONLY: AVOG ! Avogadro's number (~6e26) [/kmole]
!
! SUBROUTINES
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
!
  IMPLICIT NONE
!
! ARGUMENTS      
    INTEGER(I4),   INTENT(IN)  :: LUNCIA ! LUN for reading in data
    INTEGER(I4),   INTENT(IN)  :: ID1    ! HITRAN/RFM index of first molecule
    INTEGER(I4),   INTENT(IN)  :: ID2    ! HITRAN/RFM index of second molecule
    INTEGER(I4),   INTENT(IN)  :: NPT    ! No. array elements WNO,CIA filled
    REAL(R4),      INTENT(IN)  :: TEM    ! Temperature [K] of tabulation
    REAL(R8),      INTENT(IN)  :: WNL    ! Lower Wno limit [cm-1] of data
    REAL(R8),      INTENT(IN)  :: WNU    ! Upper Wno limit [cm-1] of data
    LOGICAL,       INTENT(OUT) :: FAIL   ! T=a fatal error has been detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    REAL(R8), PARAMETER :: DAVOG2 = DBLE ( AVOG )**2 ! (Avog.no)^2
!
! LOCAL VARIABLES
    INTEGER(I4) :: IOS  ! Saved value of IOSTAT for error messages
    INTEGER(I4) :: IPT  ! Counter for points within tabulation
    REAL(R8)    :: DCIA ! Abs.coefficient [cm^5/molec^2] read from .cia file
    TYPE(CIATYP), ALLOCATABLE :: CIASAV(:)   ! Saved CIA during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( ALLOCATED ( CIA ) ) CALL MOVE_ALLOC ( CIA, CIASAV )
  NCIA = NCIA + 1
  ALLOCATE ( CIA(NCIA) )
  IF ( ALLOCATED ( CIASAV ) ) CIA(1:NCIA-1) = CIASAV
!
  CIA(NCIA)%ID1 = ID1
  CIA(NCIA)%ID2 = ID2
  CIA(NCIA)%NPT = NPT
  CIA(NCIA)%TEM = TEM
  CIA(NCIA)%WNL = WNL
  CIA(NCIA)%WNU = WNU
  CIA(NCIA)%IG1 = IDXGAS ( ID1 ) 
  CIA(NCIA)%IG2 = IDXGAS ( ID2 ) 
  ALLOCATE ( CIA(NCIA)%WNO(NPT), CIA(NCIA)%ABS(NPT) )
  DO IPT = 1, NPT
    READ ( LUNCIA, *, IOSTAT=IOS, ERR=900 ) CIA(NCIA)%WNO(IPT), DCIA
    CIA(NCIA)%ABS(IPT)  = MAX ( 0.0, SNGL ( DCIA * DAVOG2 ) ) 
  END DO
!
900 CONTINUE
  FAIL = IOS .NE. 0
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-REACIA: I/O error reading CIA file, IOSTAT=', IOS
!
END SUBROUTINE REACIA
END MODULE REACIA_SUB
