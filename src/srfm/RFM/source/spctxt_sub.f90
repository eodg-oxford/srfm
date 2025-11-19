MODULE SPCTXT_SUB
CONTAINS
SUBROUTINE SPCTXT ( RECORD, FAIL, ERRMSG )
!
! VERSION
!   24SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Read Spectral Range & Resln from text string
!   Called by DRVSPC
!   Assumes that Spectral range label has just been read and NSPC set to
!   index for new range data.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE SPCCOM_DAT ! Spectral range data
!
  IMPLICIT NONE 
!
! ARGUMENTS 
    CHARACTER(*),  INTENT(IN)  :: RECORD ! Record containing range & resln info
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER :: IOS ! Saved value of IOSTAT for error message
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  READ ( RECORD, *, IOSTAT=IOS ) SPC(NSPC)%WNL, SPC(NSPC)%WNU, SPC(NSPC)%WNR
  FAIL = IOS .NE. 0 
  IF ( FAIL ) THEN
    IF ( SPC(NSPC)%LAB .EQ. '' ) THEN
      WRITE ( ERRMSG, * )  &
        'F-SPCTXT: Error reading Rng/Res params. IOSTAT=', IOS
    ELSE
      WRITE ( ERRMSG, * )  'F-SPCTXT: Error reading Rng/Res params ' // &
        'for Range=' // TRIM(SPC(NSPC)%LAB) // '. IOSTAT=', IOS
    END IF
    RETURN
  END IF
!
  IF ( SPC(NSPC)%WNR .LE. 0.0D0 ) THEN
    FAIL = .TRUE.
    IF ( SPC(NSPC)%LAB .EQ. '' ) THEN 
      ERRMSG = 'F-SPCTXT: Spectral Resln must be .GE. 0'
    ELSE
      ERRMSG = 'F-SPCTXT: Spectral Resln for Range=' // TRIM(SPC(NSPC)%LAB) &
              // ' must be .GE. 0'
    END IF
    RETURN
  END IF
!
  SPC(NSPC)%NPT =  1 + NINT ( (SPC(NSPC)%WNU - SPC(NSPC)%WNL) / SPC(NSPC)%WNR )
!
END SUBROUTINE SPCTXT
END MODULE SPCTXT_SUB
