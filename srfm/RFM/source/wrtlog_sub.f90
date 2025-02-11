MODULE WRTLOG_SUB
CONTAINS
SUBROUTINE WRTLOG ( MESSGE, LACCUM )
!
! VERSION
!   17MAY24 AD Checked.
!   22NOV17 AD Temporary fix avoiding allocatable character strings
!   01MAY17 AD F90 conversion of rfmlog.for. Tested.
!
! DESCRIPTION
!   Write text message to log file
!   General purpose module.
!   Limits output to 80 characters, uses '...' if continuation required
!   If argument LACCUM is present:
!     If TRUE, the log message is stored or appended, after an inserted space, 
!       to any currently stored rather than output
!     If FALSE, the accumulated message is written to the log file
!   Note: also possible to directly output messages while another message is 
!   being accumulated.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE RFMLUN_DAT, ONLY: LUNLOG ! LUN for log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),      INTENT(IN) :: MESSGE ! Message to be written
    LOGICAL, OPTIONAL, INTENT(IN) :: LACCUM ! T=accumulate, F=write 
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: NCH = 500 ! Max length of output message 
!
! LOCAL VARIABLES
    INTEGER(I4) :: IACC     ! End of current line of text
    INTEGER(I4) :: IBRK     ! Location to insert line break in text
    INTEGER(I4) :: IEND     ! Pointer to location of last character in MESSGE
    INTEGER(I4) :: IOS      ! Saved value of IOSTAT for error messages
    INTEGER(I4) :: ISPC     ! Location of space character
    INTEGER(I4) :: ISTA     ! Pointer to start of MESSGE to be written out next
    CHARACTER(NCH), SAVE :: ACCTXT = '' ! Accumulated text 
    CHARACTER(NCH), SAVE :: OUTTXT = '' ! Output text
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IEND = LEN_TRIM ( MESSGE ) 
!  OUTTXT = MESSGE(1:IEND)
  IF ( PRESENT ( LACCUM ) ) THEN
    IACC = LEN_TRIM ( ACCTXT ) 
    ACCTXT = ACCTXT(1:IACC) // ' ' // MESSGE(1:IEND)
    IF ( LACCUM ) RETURN             ! Normal exit with saved text
    OUTTXT = ACCTXT
    ACCTXT = ''
  ELSE
    OUTTXT = MESSGE(1:IEND)
  END IF
!
  IEND = LEN_TRIM ( OUTTXT )
  ISTA = 1
  DO
    IBRK = ISTA + 79
    IF ( IEND .LE. IBRK ) EXIT
! Look for suitable point for line break
    ISPC = INDEX ( OUTTXT(ISTA:IBRK), ' ', .TRUE. ) ! T=reverse search
    IF ( ISPC .GT. 50 ) THEN                        ! insert line break at space
      IBRK = ISTA + ISPC - 2   ! location of last character before space
      WRITE ( LUNLOG, '(A)', IOSTAT=IOS, ERR=900 ) OUTTXT(ISTA:IBRK)  
    ELSE                                            ! no space, break in string
      IBRK = ISTA + 76
      WRITE ( LUNLOG, '(A)', IOSTAT=IOS, ERR=900 ) OUTTXT(ISTA:IBRK)//'...'
    END IF
    ISTA = IBRK + 1
  END DO
  WRITE ( LUNLOG, '(A)', IOSTAT=IOS, ERR=900 ) OUTTXT(ISTA:IEND)
!
900 CONTINUE
  IF ( IOS .NE. 0 ) THEN
    WRITE (*,*) 'F-WRTLOG: Failed to write to rfm.log file. IOSTAT=', IOS
    STOP
  END IF
!
! Normal exit with text written out
!
END SUBROUTINE WRTLOG
END MODULE WRTLOG_SUB

