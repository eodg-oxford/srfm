MODULE SFLCOM_DAT
!
! VERSION
!   11JUN24 AD Checked.
!   19DEC17 AD Adapted from LFLCOM. Checked.
!
! DESCRIPTION
!   SVD-compressed LUT files
!   Initialised in SVDFIL.
!   IDXSFL allocated in DRVSVD.
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
  TYPE :: SFLTYP
    INTEGER(I4)       :: IGS ! Absorber index
    INTEGER(I4)       :: NL  ! No. S.Vs qualifier (0=not defined)
    CHARACTER(LENREC) :: NAM ! Name of SVD-LUT file
  END TYPE SFLTYP
!
! GLOBAL VARIABLES
    TYPE(SFLTYP), ALLOCATABLE :: SFL(:)
!
    INTEGER(I4)              :: NSFL = 0    ! No. of SVD-LUT files to be used
    INTEGER(I4), ALLOCATABLE :: IDXSFL(:,:) ! [ngas,nspc] indices of files
!
END MODULE SFLCOM_DAT

