MODULE HDBCOM_DAT
!
! VERSION
!   08AUG24 AD Checked
!   11AUG23 AD Change HDB to structure
!   30MAY23 AD Original.
!
! DESCRIPTION
!   HDB data structure
!   Required for reading HITRAN data base files
!   Variables are set in HDBCHK.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL CONSTANTS
    INTEGER(I4), PARAMETER :: MAXLEN = 16 ! max length of any field
!
  TYPE :: HDBTYP
    INTEGER(I4)   :: IST    ! Index of line parameter set
    INTEGER(I4)   :: LEN    ! Length of record in HITRAN database file
    INTEGER(I4)   :: NFD    ! No. useful fields in database file
!    INTEGER(I4)   :: NPR    ! No. useful parameters in database file
    CHARACTER(14) :: IDMFMT ! Format string for reading IDM from database file
    CHARACTER(14) :: WNOFMT ! Format string for reading WNO from database file
    INTEGER(I4), POINTER :: IPR(:)  ! IPAR associated with each field
    INTEGER(I4), POINTER :: IPT(:)  ! Start of field in database record
    INTEGER(I4), POINTER :: JPT(:)  ! End of field in database record
  END TYPE HDBTYP
!
END MODULE HDBCOM_DAT
