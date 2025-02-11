MODULE SVDHDR_SUB
CONTAINS
SUBROUTINE SVDHDR ( RECORD, TAB, IDXMOL, IDXISO, FAIL, ERRMSG ) 
!
! VERSION
!   06NOV23 AD Checked.
!   01MAY17 AD F90 conversion of svdinf.for. Checked.
!
! DESCRIPTION
!   Read SVD-LUT header record
!   Called by SVDINF, REASVD for each SVD LUT file
!   Opens file, reads header, and leaves pointer ready to read matrix data.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE IDGNEW_FNC ! Convert old RFM index for .xsc data to new value
    USE UPCASE_FNC ! Convert text string to upper case
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(80), INTENT(IN)  :: RECORD ! Header record from SVD-LUT file
    CHARACTER(3),  INTENT(OUT) :: TAB    ! Tab. function ('LOG', 'LIN' or '4RT')
    INTEGER(I4),   INTENT(OUT) :: IDXMOL ! RFM/HITRAN ID of gas
    INTEGER(I4),   INTENT(OUT) :: IDXISO ! Isotop.ID, or -1 if not isotopic.
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IOS ! Saved value of IOSTAT for error messages
!
! EXECUTABLE CODE -------------------------------------------------------------
!
!  MWCODE = RECORD(1:8)
!     
  READ ( RECORD(10:11), *, IOSTAT=IOS, ERR=900 ) IDXMOL
! for historical reasons SVD LUTs use 30 for SF6
  if ( idxmol .eq. 30 ) idxmol = 64  
  IDXMOL = IDGNEW ( IDXMOL )
  IF ( RECORD(12:12) .EQ. '.' ) THEN
    READ ( RECORD(13:13), *, IOSTAT=IOS, ERR=900 ) IDXISO
    TAB = UPCASE ( RECORD(15:17) ) 
  ELSE
    IDXISO = -1                        ! Flag for no isotope separation
    TAB = UPCASE ( RECORD(13:15) ) 
  END IF
!
  IF ( TAB .NE. 'LOG' .AND. TAB .NE. 'LIN' .AND. TAB .NE. '4RT' ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-SVDHDR: SVD-LUT file contains unrecognised tab.fnc=' // TAB
    RETURN
  END IF
!
900 CONTINUE
  FAIL = IOS .NE. 0
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-SVDHDR: Read failure on SVD-LUT file header. IOSTAT=', IOS
!
END SUBROUTINE SVDHDR
END MODULE SVDHDR_SUB
