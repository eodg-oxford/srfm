MODULE DRVJAC_SUB
CONTAINS
SUBROUTINE DRVJAC ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   22APR24 AD Checked.
!   05MAR19 AD Bug fix - ensure JDX=0 is defined at start. Checked.
!   01MAY17 AD F90 conversion of inpjac.for. Checked.
!
! DESCRIPTION
!   Read RFM driver table *JAC section
!   Called once by RFMDRV if JACFLG set TRUE.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE FLGCOM_DAT, ONLY: JTPFLG ! T = Jacobians for tan.pt. pertubrations only
!
! SUBROUTINES
    USE JACALT_SUB ! Set Jacobian perturbation altitudes
    USE JACTAN_SUB ! Pass list of tangent heights to Jacobian check
    USE JACTGT_SUB ! Check valid Jacobian target
    USE LOCASE_FNC ! Convert text string to lower case
    USE NXTFFL_SUB ! Load next field from rfm.drv, expanding filenames
    USE UPCASE_FNC ! Convert text string to upper case
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
    LOGICAL           :: GOTTAN ! T=using TAN list of altitudes
    INTEGER(I4)       :: IOS    ! Saved value of IOSTAT for error messages
    INTEGER(I4)       :: JDX=0  ! Index of Jacobian target 
    INTEGER(I4)       :: LENGTH ! Length of field from driver table or jac file
    INTEGER(I4)       :: NFLD   ! Counter for number of altitude fields
    REAL(R4)          :: ALT    ! Perturbation altitude [km] for Jacobians
    CHARACTER(LENREC) :: FIELD  ! Field extracted from driver table record
    CHARACTER(LENREC) :: TARGET ! Target species for Jacobian
!
! EXECUTABLE CODE ------------------------------------------------------------
!
! Each record should start with the name of target (eg molecule, 'TEM' etc)
  DO
    CALL NXTFFL ( LUNDRV, TARGET, LENGTH, FAIL, ERRMSG, .TRUE. ) 
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT      ! end if *JAC section reached
!
! Check target and add to list of Jacobians
    TARGET = LOCASE ( TARGET ) 
    CALL JACTGT ( TARGET, JDX, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    CALL WRTLOG ( 'I-DRVJAC: ' // TARGET, .TRUE. ) 
!
    NFLD = 0
    GOTTAN = .FALSE.
    DO                              ! altitude fields (if any) from same record
      CALL NXTFFL ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG, .TRUE. )
      IF ( FAIL ) RETURN
      IF ( LENGTH .EQ. 0 ) EXIT   ! end of record reached
      NFLD = NFLD + 1
      IF ( UPCASE ( FIELD ) .EQ. 'TAN' ) THEN
        CALL JACTAN ( JDX, TARGET, NFLD, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
        GOTTAN = .TRUE.
      ELSE 
        READ ( FIELD, *, IOSTAT=IOS ) ALT
        IF ( IOS .NE. 0 ) THEN
          FAIL = .TRUE.
          LENGTH = MIN ( LENGTH, 15 )     ! Restrict for error message
          ERRMSG = 'F-DRVJAC: Error reading *JAC altitude for ' // & 
                   TRIM(TARGET) // ' from field=' // TRIM(FIELD(1:LENGTH))
        ELSE
          CALL JACALT ( JDX, TARGET, .FALSE., ALT, FAIL, ERRMSG ) 
        END IF
        IF ( FAIL ) RETURN
        CALL WRTLOG ( FIELD, .TRUE. ) 
      END IF
    END DO
    CALL JACALT ( JDX, TARGET, .TRUE., 0.0, FAIL, ERRMSG )  ! Finish alts.
    IF ( FAIL ) RETURN
! If using 'TAN' for list of altitudes, must be middle of 3 fields
    IF ( GOTTAN .AND. NFLD .NE. 3 ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-DRVJAC: Jacobian ''TAN'' altitudes not ' // &
               'of form ''zmin TAN zmax'''
      RETURN
    END IF
! If using JTP Flag, then 'TAN' for list of altitudes is mandatory
    IF ( JTPFLG .AND. .NOT. GOTTAN ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-DRVJAC: JTP flag requires *JAC section ' // &
               'altitudes of the form ''zmin TAN zmax'''
      RETURN
    END IF
    CALL WRTLOG ( '', .FALSE. )
  END DO
!
  IF ( JDX .EQ. 0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVJAC: No data found in *JAC section of driver table'
    RETURN
  END IF
!
END SUBROUTINE DRVJAC
END MODULE DRVJAC_SUB

