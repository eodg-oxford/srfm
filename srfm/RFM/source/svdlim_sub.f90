MODULE SVDLIM_SUB
CONTAINS
SUBROUTINE SVDLIM ( NLREQ, NL )
!
! VERSION
!   30JUN24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Limit No. Singular Vectors used in SVD-LUT
!   Called by REASVD if qualifier attached to filename in *LUT section
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN)    :: NLREQ ! Requirement from file qualifier
    INTEGER(I4), INTENT(INOUT) :: NL    ! No. of S.V.s to use
!
! LOCAL VARIABLES
    CHARACTER(11) :: REQSTR ! String for required number of SVs
    CHARACTER(11) :: FILSTR ! String for number of SVs in file
    CHARACTER(80) :: LOGMSG ! Message sent to log file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  REQSTR = C11INT ( ABS ( NLREQ ) )
  FILSTR = C11INT ( NL ) 
!
  IF ( NLREQ .GT. NL ) THEN
    LOGMSG = 'W-SVDLIM: Including all ' // TRIM ( FILSTR ) // &
             ' basis vectors in file (<No.Rqd=' // TRIM ( REQSTR ) // ')' 
!
  ELSE IF ( NLREQ .GT. 0 ) THEN
    LOGMSG = 'I-SVDLIM: Including ' // TRIM ( REQSTR ) // ' out of ' // &
             TRIM ( FILSTR ) // ' basis vectors from LUT file' 
    NL = NLREQ
!
  ELSE IF ( NLREQ .GE. -NL ) THEN   ! -NL < NLREQ < 0
    LOGMSG = 'I-SVDLIM: Excluding ' // TRIM ( REQSTR ) // ' out of ' // &
             TRIM ( FILSTR ) // ' basis vectors from LUT file' 
    NL = NL + NLREQ
!
  ELSE                                        ! NLREQ < -NL
    LOGMSG = 'W-SVDLIM: Excluding all ' // TRIM ( FILSTR ) // &
             ' basis vectors in file (<No.Rqd=' // TRIM ( REQSTR ) // ')' 
    NL = 0
  END IF
!
  CALL WRTLOG ( LOGMSG )
!
END SUBROUTINE SVDLIM
END MODULE SVDLIM_SUB
