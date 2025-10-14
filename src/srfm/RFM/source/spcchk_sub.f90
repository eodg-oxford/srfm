MODULE SPCCHK_SUB
CONTAINS
SUBROUTINE SPCCHK ( FAIL, ERRMSG )
!
! VERSION
!   17MAR23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Check spectral range and resolution data
!   Called by DRVSPC for each spectral range.
!   Assumes latest range to be checked is loaded in SPCCOM as NSPC
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE SPCCOM_DAT ! Spectral range data
    USE PHYCON_DAT, ONLY: GHZ2CM ! GHz to cm-1 conv.fac
!
! SUBROUTINES
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL CONSTANTS
    REAL(R8), PARAMETER :: WNUMAX = 50000.0D0 ! Max Wno allowed (50000 for UV)
    REAL(R8), PARAMETER :: WNLMIN = 0.001D0   ! Min Wno allowed
!
! LOCAL VARIABLES
    REAL(R8)      :: DELWID ! Widemesh interval [cm-1] for current range
    REAL(R8)      :: RPT    ! (No.spectral points - 1) as a real number
    CHARACTER(31) :: ERRSTR ! Common text component of error message
    CHARACTER(80) :: MESSGE ! Info message for Log file
    CHARACTER(3)  :: WNOGHZ ! Either 'Wno' (GHZFLG=F) or 'GHz' (GHZFLG=T)
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( GHZFLG ) THEN 
    WNOGHZ = 'GHz'
  ELSE
    WNOGHZ = 'Wno'
  END IF
!
  FAIL = .TRUE.                      ! Set temporarily
  ERRSTR = 'F-SPCCHK: Range Label=' // SPC(NSPC)%LAB //' '  ! C*31
  IF ( SPC(NSPC)%WNL .GT. SPC(NSPC)%WNU ) THEN
    WRITE ( ERRMSG, '(A,G12.5,A,G12.5)' ) ERRSTR // &
            ' Lower ' // WNOGHZ // '=', SPC(NSPC)%WNL, &
            ', > Upper ' // WNOGHZ // '=', SPC(NSPC)%WNU
  ELSE IF ( SPC(NSPC)%WNL .LT. WNLMIN ) THEN
    WRITE ( ERRMSG, '(A,G12.5,A,G12.5)' ) ERRSTR // ' Lower ' // & 
            WNOGHZ // '=', SPC(NSPC)%WNL, ', < Minimum=', WNLMIN
  ELSE IF ( SPC(NSPC)%WNU .GT. WNUMAX ) THEN
    WRITE ( ERRMSG, '(A,G12.5,A,G12.5)' ) ERRSTR // ' Upper ' // &
            WNOGHZ // '=', SPC(NSPC)%WNU, ', > Maximum=', WNUMAX
  ELSE IF ( SPC(NSPC)%WNR .LT. 0.0 ) THEN
    WRITE ( ERRMSG, '(A,G12.5,A)' ) ERRSTR // & 
           ' Reqd Resln =', SPC(NSPC)%WNR, ' is negative'
  ELSE IF ( .NOT. GHZFLG .AND. SPC(NSPC)%WNR .GT. 1.0 .AND. &
             MOD ( SPC(NSPC)%WNR, 1.D0 ) .NE. 0.0 ) THEN 
    WRITE ( ERRMSG, '(A,G12.5)' ) ERRSTR // &
            ' non-integer for No.pts/Wno,=', SPC(NSPC)%WNR
  ELSE IF ( SPC(NSPC)%WNR .EQ. 0.0D0 .AND. ( ILSFLG .OR. AVGFLG ) ) THEN
    ERRMSG = ERRSTR // ' cannot use irreg output grid with ILS,AVG flags'
  ELSE 
    FAIL = .FALSE.
  END IF
  IF ( FAIL ) RETURN
!
! If WNR expressed as points/cm-1 convert to wavenumber spacing
  IF ( .NOT. GHZFLG .AND. SPC(NSPC)%WNR .GT. 1.0D0 ) &
       SPC(NSPC)%WNR = 1.0D0 / SPC(NSPC)%WNR 
!
! Check that upper wavenumber complies with lower wavenumber and resolution:
! adjust upper wavenumber if outside 0.01 of resolution and send warning.
  IF ( SPC(NSPC)%WNR .GT. 0.0D0 ) THEN      ! regular grid       
    RPT = ( SPC(NSPC)%WNU - SPC(NSPC)%WNL ) / SPC(NSPC)%WNR
    IF ( ABS ( RPT - DBLE(NINT(RPT)) ) .GT. 0.01D0 * SPC(NSPC)%WNR ) THEN
      SPC(NSPC)%WNU = SPC(NSPC)%WNL + NINT(RPT) * SPC(NSPC)%WNR
      WRITE ( MESSGE, '(A,F12.6)' ) 'W-SPCCHK: Label=' // &
              SPC(NSPC)%LAB // ' Adjusting upper ' // WNOGHZ // & 
              ' boundary to ', SPC(NSPC)%WNU
      CALL WRTLOG ( MESSGE  )
    END IF
  END IF
!
! Basic range, resolution parameters OK
!
! Send log message indicating current range identified
  WRITE ( MESSGE, '(A,3F12.5)' ) 'I-SPCCHK: Label=' // SPC(NSPC)%LAB // &
         ' Range, Resln=', SPC(NSPC)%WNL, SPC(NSPC)%WNU, SPC(NSPC)%WNR
  CALL WRTLOG ( MESSGE )
!
! Convert from GHz to Wno if necessary
  IF ( GHZFLG ) THEN
    SPC(NSPC)%WNL = SPC(NSPC)%WNL * GHZ2CM 
    SPC(NSPC)%WNU = SPC(NSPC)%WNU * GHZ2CM  
    SPC(NSPC)%WNR = SPC(NSPC)%WNR * GHZ2CM  
  END IF
!
  IF ( AVGFLG ) THEN
    SPC(NSPC)%WXL = SPC(NSPC)%WNL - SPC(NSPC)%WNR
    SPC(NSPC)%WXU = SPC(NSPC)%WNU + SPC(NSPC)%WNR
  ELSE IF ( ILSFLG ) THEN  ! Temporarily estimate ILS width=widemesh interval
    IF ( SPC(NSPC)%WNR .EQ. 0.0D0 ) THEN
      DELWID = 1.0D0
    ELSE
      DELWID = SPC(NSPC)%WNR * NINT ( 1.0D0 / SPC(NSPC)%WNR )
    END IF
    SPC(NSPC)%WXL = SPC(NSPC)%WNL - DELWID
    SPC(NSPC)%WXU = SPC(NSPC)%WNU + DELWID
  ELSE
    SPC(NSPC)%WXL = SPC(NSPC)%WNL 
    SPC(NSPC)%WXU = SPC(NSPC)%WNU
  END IF
!
! Update Min/Max wavenumbers required for any spectral range
  WMNSPC = MINVAL ( SPC%WXL )
  WMXSPC = MAXVAL ( SPC%WXU )
!
END SUBROUTINE SPCCHK
END MODULE SPCCHK_SUB
