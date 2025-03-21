MODULE DRVREJ_SUB
CONTAINS
SUBROUTINE DRVREJ ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   13JUL24 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Read RFM driver table *REJ section
!   Called by RFMDRV following *REJ marker in Driver table.
!   If REJFLG is TRUE read modified line rejection parameters
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data
    USE LENREC_DAT ! Max length of input text record
    USE REJCOM_DAT ! Minimum line strength limits
!
! SUBROUTINES
    USE CHKGAS_SUB ! Check for valid molecule name, isotope, Vib.Level
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
    USE NXTREC_SUB ! Load next record from input file
    USE TXTFLD_SUB ! Identify start and end points of text field in record
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
    LOGICAL           :: SETALL = .FALSE.  ! T=wildcard '*' found
    LOGICAL           :: ENDSEC ! T=end of section reached
    INTEGER(I4)       :: IDXISO ! HITRAN index of isotope
    INTEGER(I4)       :: IDXMOL ! HITRAN/RFM index of molecule
    INTEGER(I4)       :: IEND   ! End position of field in RECORD
    INTEGER(I4)       :: IGAS   ! Index of molecule,isotope in GASCOM
    INTEGER(I4)       :: IOS    ! Saved value of IOSTAT for error message
    INTEGER(I4)       :: ISTA   ! Start position of field in RECORD
    REAL(R4)          :: REJ    ! Min line strength for inclusion
    CHARACTER(LENGAS) :: MOLEC  ! Molecule name, or '*'
    CHARACTER(LENREC) :: RECORD ! Record read from driver table *REJ section
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  USEREJ = .TRUE.
  ALLOCATE ( STRREJ(NGAS) )
  STRREJ = -1.0               ! flag for 'not set'
!
  DO                          ! loop over each record in section
    CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( ENDSEC ) EXIT        ! end of section reached
! Field#1 is molec, or '*'
    CALL TXTFLD ( RECORD, 1, ISTA, IEND )    
    IF ( RECORD(ISTA:IEND) .EQ. '*' ) THEN
      IF ( SETALL ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-DRVREJ: Repeated ''*'' in *REJ section'
        RETURN
      END IF
      MOLEC = '*'
      SETALL = .TRUE.
    ELSE
      CALL CHKGAS ( RECORD(ISTA:IEND), IDXMOL, IDXISO ) 
      IF ( IDXMOL .EQ. 0 ) THEN
        FAIL = .TRUE.
        IEND = MIN ( ISTA+20, IEND ) ! Limit to fit into error message
        ERRMSG = 'F-DRVREJ: Unrecognised molecule=' // RECORD(ISTA:IEND)
        RETURN
      ELSE IF ( IDXMOL .GT. MAXHLN ) THEN  ! MAXHLN in GASCOM
        FAIL = .TRUE.
        ERRMSG = 'F-DRVREJ: ' // RECORD(ISTA:IEND) // ' not a line molecule'
        RETURN
      END IF
      MOLEC = RECORD(ISTA:IEND) 
      IGAS = IDXGAS ( IDXMOL, IDXISO ) 
      IF ( IGAS .EQ. 0 ) THEN
        CALL WRTLOG ( 'W-DRVREJ: molecule ' // MOLEC // ' not used' ) 
        EXIT
      END IF
    END IF
! Field#2 is minimum line strength for inclusion
    CALL TXTFLD ( RECORD, 2, ISTA, IEND ) 
    READ ( RECORD(ISTA:IEND), *, IOSTAT=IOS ) REJ
    IF ( IOS .NE. 0 ) THEN
      FAIL = .TRUE.
      WRITE ( ERRMSG, * ) 'F-DRVREJ: Failed to read rej.limit for molec' // &
        MOLEC // ' IOSTAT=', IOS
      RETURN
    ELSE IF ( REJ .LT. 0 ) THEN
      FAIL = .TRUE.
      WRITE ( ERRMSG, * ) 'F-DRVREJ: Negative values not allowed, =', REJ
      RETURN
    END IF
!
    IF ( MOLEC .EQ. '*' ) THEN
      WHERE ( STRREJ .EQ. -1 ) STRREJ = REJ
    ELSE
      STRREJ(IGAS) = REJ
    END IF
!
  END DO
!
! If no wildcard, set value 0.0 for other molecules, ie include all lines.
  WHERE ( STRREJ .EQ. -1 ) STRREJ = 0.0
!
END SUBROUTINE DRVREJ
END MODULE DRVREJ_SUB

