MODULE DRVNTE_SUB
CONTAINS
SUBROUTINE DRVNTE ( LUNDRV, FAIL, ERRMSG )
!
! VERSION
!   20SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Read RFM driver table *NTE section
!   Called once by RFMDRV.
!   Expects to find set of filenames of Vibrational Temperatures.
!   This may add isotopic profiles which have not already been loaded in the 
!   *ATM section.
!
! VARIABLE KINDS
    USE KIND_DAT
! 
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE FLGCOM_DAT, ONLY: GRAFLG ! T=2D Atmosphere
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE ATMGRA_SUB ! Interpolate to fill 2-D atmospheric profile field
    USE ATMPSI_SUB ! Extract Psi angle from brackets following .atm filename
    USE NTEDEF_SUB ! Use default .nte filename to find any missing files
    USE NTEFIL_SUB ! Read NTE Vibrational Temperature Data from file
    USE NXTFLD_SUB ! Load next field from section of driver file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: LUNDRV ! LUN for Driver File
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)       :: LENGTH      ! Length of filename read from driver file
    CHARACTER(LENREC) :: NAMDEF = '' ! Default filename template
    CHARACTER(LENREC) :: NAMNTE      ! Name of NTE file read from Driver Table
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO                        ! repeat for each .nte filename  in *NTE section
    CALL NXTFLD ( LUNDRV, NAMNTE, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT 
!
    IF ( INDEX ( NAMNTE(1:LENGTH), '*' ) .EQ. 0 ) THEN
      CALL ATMPSI ( NAMNTE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      CALL NTEFIL ( NAMNTE(1:LENGTH), FAIL, ERRMSG )
      IF ( FAIL ) RETURN
    ELSE IF ( NAMDEF .EQ. '' ) THEN
      NAMDEF = NAMNTE(1:LENGTH)
    ELSE
      FAIL = .TRUE.
      ERRMSG = 'F-DRVNTE: *NTE section has more than one filename template'
      RETURN
    END IF
  END DO
!
  IF ( NAMDEF .NE. '' ) THEN
    CALL ATMPSI ( NAMNTE, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    CALL NTEDEF ( NAMDEF, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
  END IF
!
  IF ( GRAFLG ) CALL ATMGRA 
!
END SUBROUTINE DRVNTE
END MODULE DRVNTE_SUB
