MODULE WRTSPC_SUB
CONTAINS
SUBROUTINE WRTSPC ( LUN, NPT, IRREG, WNO, SPC, FAIL, ERRMSG )
!
! VERSION
!   04JAN23 AD Checked.
!   01JUL17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Write spectral data to file
!   Called by SPCWRT for each output file.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE PHYCON_DAT, ONLY: GHZ2CM ! GHz to cm-1 conv.fac
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUN    ! Logical unit number
    INTEGER(I4),   INTENT(IN)  :: NPT    ! No. points to be written
    LOGICAL,       INTENT(IN)  :: IRREG  ! T=irreg grid, F=regular grid
    REAL(R8),      INTENT(IN)  :: WNO(:) ! Spectral axis [cm-1]
    REAL(R8),      INTENT(IN)  :: SPC(:) ! Output spectrum
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error occurs
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    LOGICAL, PARAMETER :: OLDFMT = .FALSE. ! T=write old-format binary files
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IOS    ! Value of IOSTAT saved for error messages.
    INTEGER(I4)   :: IPT    ! Counter for spectral points
    REAL(R8)      :: WNOFAC ! Scaling to convert [cm-1] to [GHz], or not
    CHARACTER(10) :: WNOSTR ! String containing output Wno.
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( IRREG ) THEN
    IF ( GHZFLG ) THEN
      WNOFAC = 1.0D0 / GHZ2CM 
    ELSE
      WNOFAC = 1.0D0
    END IF
!
    DO IPT = 1, NPT
      IF ( BINFLG ) THEN
        IF ( DBLFLG ) THEN
          WRITE ( LUN, IOSTAT=IOS, ERR=900 )  WNO(IPT)*WNOFAC, SPC(IPT)
        ELSE
          WRITE ( LUN, IOSTAT=IOS, ERR=900 ) WNO(IPT)*WNOFAC, SNGL(SPC(IPT))
        END IF
      ELSE   
        WRITE ( WNOSTR, '(F10.4)' ) WNO(IPT)*WNOFAC
        IF ( DBLFLG ) THEN
          WRITE ( LUN, *, IOSTAT=IOS, ERR=900 ) WNOSTR, SPC(IPT)
        ELSE
          WRITE ( LUN, *, IOSTAT=IOS, ERR=900 ) WNOSTR, SNGL(SPC(IPT))
        END IF
      END IF
    END DO
!
  ELSE
    IF ( BINFLG ) THEN
      IF ( DBLFLG ) THEN
        IF ( OLDFMT ) THEN
          WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NPT, SPC(1:NPT)
        ELSE
          WRITE ( LUN, IOSTAT=IOS, ERR=900 ) SPC(1:NPT)
        END IF
      ELSE
        IF ( OLDFMT ) THEN
          WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NPT, SNGL ( SPC(1:NPT) )
        ELSE
          WRITE ( LUN, IOSTAT=IOS, ERR=900 ) SNGL ( SPC(1:NPT) ) 
        END IF
      END IF
    ELSE   
      IF ( DBLFLG ) THEN
        WRITE ( LUN, *, IOSTAT=IOS, ERR=900 ) SPC(1:NPT)
      ELSE
        WRITE ( LUN, *, IOSTAT=IOS, ERR=900 ) SNGL ( SPC(1:NPT) ) 
      END IF
    END IF
  END IF
!
900  CONTINUE
  FAIL = IOS .NE. 0
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-WRTSPC: Error writing output spectrum. IOSTAT=', IOS
!
END SUBROUTINE WRTSPC
END MODULE WRTSPC_SUB
