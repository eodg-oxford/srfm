MODULE DRVDIM_SUB
CONTAINS
SUBROUTINE DRVDIM ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   25OCT23 AD Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Read RFM driver table *DIM section
!   Called by RFMDRV once if TABFLG set.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE TABCOM_DAT ! Axes for tabulated absorption coefficients
    USE LENREC_DAT ! Max length of input text record
!
! SUBROUTINES
    USE DIMGRD_SUB ! Read TAB axis grid parameters from driver file
    USE DIMFIL_SUB ! Read list of TAB axis values from file
    USE DIMPCG_SUB ! Set TAB p-axis to CG values of internal ATM profile
    USE DIMPLV_SUB ! Set TAB p-axis to internal PRE ATM profile
    USE DIMTAB_SUB ! Complete TAB header information
    USE ENDCHK_SUB ! Check end of Driver Table section has been reached
    USE LEXIST_FNC ! Check if file exists
    USE NXTFLD_SUB ! Load next field from section of driver file
    USE UPCASE_FNC ! Convert text string to upper case
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for Driver File
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)       :: IFLD   ! Counter for fields in Driver Table section
    INTEGER(I4)       :: LENGTH ! Length of field read from driver table
    CHARACTER(LENREC) :: FIELD  ! Field read from driver table
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! First field is either NP or name of file containing list of pressures or
! one of 'PCG' or 'PLV' (case insensitive)
  IFLD = 1
  CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
  IF ( LENGTH .EQ. 0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVDIM: No data in *DIM/*TAN section'  
    RETURN
  END IF
!
  IF ( UPCASE ( FIELD ) .EQ. 'PCG' ) THEN
    CALL DIMPCG
  ELSE IF ( UPCASE ( FIELD ) .EQ. 'PLV' ) THEN
    CALL DIMPLV
  ELSE IF ( LEXIST ( FIELD ) ) THEN
    CALL DIMFIL ( 'P', FIELD, NPTAB, PAXTAB, FAIL, ERRMSG )
    IFLD = IFLD + 1
  ELSE
! Read 2 more fields containing lower/upper pressures
    CALL DIMGRD ( LUNDRV, 'P', IFLD, FIELD, NPTAB, PAXTAB, FAIL, ERRMSG )
  END IF
  IF ( FAIL ) RETURN
!
! Next field is either NT or name of file containing list of temperatures
  CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
  IF ( LENGTH .EQ. 0 ) THEN   ! No temperature data incl in tabulation axes
    NTTAB = 1
    ALLOCATE ( TAXTAB(1) )
    TAXTAB(1) = 0.0
  ELSE IF ( LEXIST ( FIELD ) ) THEN
    CALL DIMFIL ( 'T', FIELD, NTTAB, TAXTAB, FAIL, ERRMSG )
    IFLD = IFLD + 1
  ELSE
    CALL DIMGRD ( LUNDRV, 'T', IFLD, FIELD, NTTAB, TAXTAB, FAIL, ERRMSG )
  END IF
  IF ( FAIL ) RETURN
!
! Next field is either NQ or name of file containing list of VMR values
  CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
  IF ( LENGTH .EQ. 0 ) THEN      ! No VMR data included in tabulation axes
    NQTAB = 1
    ALLOCATE ( QAXTAB(1) )
    QAXTAB(1) = 100.0
  ELSE IF ( LEXIST ( FIELD ) ) THEN
    CALL DIMFIL ( 'Q', FIELD, NQTAB, QAXTAB, FAIL, ERRMSG )
    IFLD = IFLD + 1
  ELSE
    CALL DIMGRD ( LUNDRV, 'Q', IFLD, FIELD, NQTAB, QAXTAB, FAIL, ERRMSG )
  END IF
  IF ( FAIL ) RETURN
!
! Set up additional header info
  CALL DIMTAB
!
! Check no more fields in *DIM section
  CALL ENDCHK ( LUNDRV, '*DIM', FAIL, ERRMSG )                   
!
END SUBROUTINE DRVDIM
END MODULE DRVDIM_SUB
