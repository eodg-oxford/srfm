MODULE DRVFIN_SUB
CONTAINS
SUBROUTINE DRVFIN ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   19NOV23 AD Checked.
!   01MAY17 AD F90 conversion of inpfin.for. Tested.
!
! DESCRIPTION
!   Read RFM driver table *FIN section
!   Called by RFMDRV following *FIN marker in Driver table and FINFLG=TRUE
!   Note that this is only required if coupled with an ILS convolution -
!   otherwise the fine mesh resolution is the same as the output resolution. 
!   If a finemesh resolution is not specified, a default value is used.
!
!   Note that any user-specified resolution >1 is interpreted as pts/cm-1,
!   while values <1 are interpreted as spacing in cm-1.
!
!   If GHZ flag is enabled, then any value <1 is interpreted as GHz spacing
!   but values >1 are still pts/cm-1.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE FINCOM_DAT, ONLY: NOMFIN ! Nominal No.pts/cm-1 for fine mesh
    USE FLGCOM_DAT, ONLY: GHZFLG ! T = Use GHz spectral axis, F = Wavenumber
    USE PHYCON_DAT, ONLY: GHZ2CM ! GHz to cm-1 conv.fac (~1/30)
!
! SUBROUTINE
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE C9REAL_GEN ! Write real number as C*9 string
    USE ENDCHK_SUB ! Check end of Driver Table section has been reached
    USE NXTFLD_SUB ! Load next field from section of driver file
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for Driver File
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)       :: IOS    ! Saved value of IOSTAT for error message
    INTEGER(I4)       :: LENGTH ! No. of characters in FIELD
    REAL(R8)          :: GHZFIN ! Original resln in GHz before conversion
    REAL(R8)          :: USRFIN ! User-specified fine mesh resolution
    CHARACTER(LENREC) :: FIELD  ! Field read from driver table
    CHARACTER(80)     :: MESSGE ! Message sent to log file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
  IF ( LENGTH .EQ. 0 ) THEN
    FAIL = .TRUE. 
    ERRMSG = 'F-DRVFIN: No data in *FIN section of Driver Table'
    RETURN
  END IF
!
  READ ( FIELD, *, IOSTAT=IOS ) USRFIN
  IF ( IOS .NE. 0 ) THEN
    FAIL = .TRUE.
    WRITE ( ERRMSG, * ) 'F-DRVFIN: Error reading fine-mesh resolution.' &
      //' IOSTAT=', IOS
    RETURN
  END IF      
!
! If value < 1 and GHZ flag, read as GHz spacing and convert to cm-1 spacing
  IF ( USRFIN .LT. 1.0 .AND. USRFIN .GT. 0.0 .AND. GHZFLG ) THEN
    GHZFIN = USRFIN
    USRFIN = GHZFIN * GHZ2CM
    MESSGE = 'I-DRVFIN: Converting *FIN spacing from ' // &
      TRIM(C9REAL(GHZFIN)) // ' GHz to ' // TRIM(C9REAL(USRFIN)) // ' cm-1'
    CALL WRTLOG ( MESSGE )
  END IF
!
  IF ( USRFIN .GT. 1.0D0 ) THEN
    NOMFIN = NINT ( USRFIN )
  ELSE
    NOMFIN = NINT ( 1.D0 / USRFIN ) 
  END IF
!
  MESSGE = 'I-DRVFIN: Using fine mesh calculation @' // TRIM(C11INT(NOMFIN)) &
           // ' pts/cm-1'
  CALL WRTLOG ( MESSGE )
!
! Check no more data in *FIN section
  CALL ENDCHK ( LUNDRV, '*FIN', FAIL, ERRMSG )
!
END SUBROUTINE DRVFIN
END MODULE DRVFIN_SUB
