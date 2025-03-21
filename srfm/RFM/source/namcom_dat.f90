MODULE NAMCOM_DAT
!
! VERSION
!   24MAR23 AD Checked.
!   08NOV17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   RFM output filenames
!   The following are the default filename templates where the '*' is replaced
!   by spectral range, path etc data to distinguish each file.
!   These templates can be changed by inserting appropriate sections into the
!   RFM driver table (*ABS, *BBT etc).
!   The '.asc' component will be changed to '.bin' if binary output files are
!   selected (BIN flag).
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL CONSTANTS
    INTEGER(I4), PARAMETER :: LENNAM = 200 ! Max length of any output filename 
!
! GLOBAL VARIABLES
    CHARACTER(LENNAM), TARGET :: ABSNAM = 'abs_*.asc' ! ABS filename template
    CHARACTER(LENNAM), TARGET :: BBTNAM = 'bbt_*.asc' ! BBT filename template
    CHARACTER(LENNAM), TARGET :: COONAM = 'coo_*.asc' ! COO filename template
    CHARACTER(LENNAM), TARGET :: OPTNAM = 'opt_*.asc' ! OPT filename template
    CHARACTER(LENNAM), TARGET :: PRFNAM = 'prf*.asc'  ! PRF filename template
    CHARACTER(LENNAM), TARGET :: PTHNAM = 'pth_*.asc' ! PTH filename template
    CHARACTER(LENNAM), TARGET :: RADNAM = 'rad_*.asc' ! RAD filename template
    CHARACTER(LENNAM), TARGET :: RJTNAM = 'rjt_*.asc' ! RJT filename template
    CHARACTER(LENNAM), TARGET :: TABNAM = 'tab_*.asc' ! TAB filename template
    CHARACTER(LENNAM), TARGET :: TRANAM = 'tra_*.asc' ! TRA filename template
    CHARACTER(LENNAM), TARGET :: WIDNAM = 'wid_*.asc' ! WID filename template
    CHARACTER(LENNAM)         :: DIRNAM = ''          ! Output directory path
!
END MODULE NAMCOM_DAT
