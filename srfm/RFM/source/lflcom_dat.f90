MODULE LFLCOM_DAT
!
! VERSION
!   15JUN24 AD Checked.
!   19DEC17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Look-Up Table files
!   Initialised in LUTFIL.
!   IDXLFL allocated in DRVLUT.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: LFLTYP
    LOGICAL           :: BIN ! T=binary file, F=ASCII
    INTEGER(I4)       :: IGS ! Absorber index
    INTEGER(I4)       :: LUN ! LUN assigned to TAB LUT files
    INTEGER(I4)       :: NDP ! Undersampling factor for p-axis (1=default)
    INTEGER(I4)       :: NDT ! Undersampling factor for T-axis (1=default)
    CHARACTER(LENREC) :: NAM ! Name of LUT file
  END TYPE LFLTYP
!
! GLOBAL VARIABLES
    TYPE(LFLTYP), ALLOCATABLE :: LFL(:)
!
    INTEGER(I4) :: NLFL = 0  ! No. of LUT files to be used
    INTEGER(I4), ALLOCATABLE :: IDXLFL(:,:) ! [ngas,nspc] indices of files
!
END MODULE LFLCOM_DAT

