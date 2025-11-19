MODULE HITCOM_DAT
!
! VERSION
!   13AUG24 AD Checked.
!   11AUG23 AD Remove SET_* flags, add %IST 
!   31MAY23 AD Remove PARHIT,BASHIT,IFWDPT. Add SET_* flags.
!   29JAN20 AD Removed redundant fields, add/rename others 
!   24FEB17 AD F90 version. Checked.
!
! DESCRIPTION
!   HITRAN line data.
!   Data Type representing structure of HITRAN record. Invoke HITCOM_RESET to
!   release cached buffers when the model is reused in the same session.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
  PUBLIC :: HITCOM_RESET
!
  TYPE :: HITTYP
    INTEGER(I4)  :: IDM ! HITRAN Gas ID
    INTEGER(I4)  :: IDI ! Isotope Number (1=most abundant,2=2nd,3 etc)
    INTEGER(I4)  :: IGS ! Index of molec,iso in GAS
    INTEGER(I4)  :: IUS ! Upper state global quanta index.
    INTEGER(I4)  :: ILS ! Lower state global quanta index
    INTEGER(I4)  :: IUV ! Index of Vib.Temp profile affecting upper level
    INTEGER(I4)  :: ILV ! Index of Vib.Temp profile affecting lower level
    INTEGER(I4)  :: IST ! Index of line parameter set
    REAL(R4)     :: STR ! Line strength  [cm-1./(kg.moles.cm-2)]@296K
    REAL(R4)     :: ELS ! Lower-state energy [cm-1]
    REAL(R4)     :: HWA ! Air-broad halfwidth  (HWHM) [cm-1/atm] @ 296K
    REAL(R4)     :: HWS ! Self-broad halfwidth (HWHM) [cm-1/atm] @ 296K.
    REAL(R4)     :: TCA ! Coeff.of temp.dependence of air-broadened HW
    REAL(R4)     :: TCS ! Coeff.of temp.dependence of self-broadened HW
    REAL(R4)     :: PSA ! Transition shift due to air pressure
    REAL(R4)     :: PSS ! Transition shift due to self-pressure
    REAL(R4)     :: LMA ! Line coupling - air-broadened
    REAL(R4)     :: LMS ! Line coupling - self-broadened
    REAL(R4)     :: WGT ! Molecular weight [atomic units]
    REAL(R8)     :: WNO ! Cyclic buf of line wnos. [cm-1]
    CHARACTER(9) :: BLQ ! Lower State local quanta
    CHARACTER(9) :: ULQ ! Upper State local quanta
  END TYPE HITTYP
!
! GLOBAL VARIABLES
    TYPE(HITTYP) :: HIT
    TYPE(HITTYP), ALLOCATABLE :: CYC(:)  ! Cyclic buffer
!
    INTEGER(I4) :: ICYC1  ! Index for lowest wavenumber line
    INTEGER(I4) :: NCYC   ! Current size of CYC array
    INTEGER(I4) :: NLIN   ! No.lines currently stored
!
CONTAINS

  SUBROUTINE HITCOM_RESET()
    HIT%IDM = 0
    HIT%IDI = 0
    HIT%IGS = 0
    HIT%IUS = 0
    HIT%ILS = 0
    HIT%IUV = 0
    HIT%ILV = 0
    HIT%IST = 0
    HIT%STR = 0.0_R4
    HIT%ELS = 0.0_R4
    HIT%HWA = 0.0_R4
    HIT%HWS = 0.0_R4
    HIT%TCA = 0.0_R4
    HIT%TCS = 0.0_R4
    HIT%PSA = 0.0_R4
    HIT%PSS = 0.0_R4
    HIT%LMA = 0.0_R4
    HIT%LMS = 0.0_R4
    HIT%WGT = 0.0_R4
    HIT%WNO = 0.0_R8
    HIT%BLQ = ''
    HIT%ULQ = ''

    IF ( ALLOCATED ( CYC ) ) THEN
      DEALLOCATE ( CYC )
    END IF

    ICYC1 = 0
    NCYC  = 0
    NLIN  = 0
  END SUBROUTINE HITCOM_RESET

END MODULE HITCOM_DAT
