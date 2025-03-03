MODULE DRVGAS_SUB
CONTAINS
SUBROUTINE DRVGAS ( LUNDRV, FAIL, ERRMSG )
!
! VERSION 
!   15APR24 AD Bug#45 Set fatal error rather than warning if NGAS = 0
!   19APR22 AD Bug#33 Add 'air' if REX flag enabled. Checked.
!   24JUN19 AD Remove setting of 'air' as absorber with REX flag. Checked.
!   03OCT17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Read RFM driver table *GAS section
!   Called by RFMDRV.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE FLGCOM_DAT, ONLY: REXFLG ! T = incorporate Rayleigh extinction
    USE GASCOM_DAT, ONLY: NGAS   ! No. different absorbers
    USE IDXCON_DAT, ONLY: IDXAIR ! RFM index for 'air'
!
! SUBROUTINES
    USE ADDGAS_SUB ! Add new molecule/isotope to list of absorbers
    USE CHKCO2_SUB ! Check if CO2 being used with CHI,MIX flags
    USE GASALL_SUB ! Add wildcard gases to list of absorbers
    USE GASCHK_SUB ! Check Gas name and set indices
    USE NAMGAS_FNC ! Return molecule name + (iso) associated with GASCOM index
    USE NXTFFL_SUB ! Load next field from rfm.drv, expanding filenames
    USE REAQAL_SUB ! Decode Qualifier strings '(IQAL)' or '(IQAL:JQAL)'
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
    LOGICAL     :: NOTGAS           ! Set TRUE if gas not recognised
    LOGICAL     :: USEALL = .FALSE. ! Set TRUE if wildcard found in section
    INTEGER(I4) :: IGAS             ! Gas counter
    INTEGER(I4) :: IQAL             ! Qualifier string attached to '*'
    INTEGER(I4) :: LENGTH           ! Length of field read from Driver file
    CHARACTER(LENREC) :: FIELD      ! Molecule or name of .gas file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO
! Load next field (if the first field is a filename, also open file using
! LUNTMP and get the first field in the file)
    CALL NXTFFL ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT
!
    IF ( FIELD(1:1) .EQ. '*' ) THEN
      IF ( USEALL ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-DRVGAS: ''*'' appears twice in *GAS section'
        RETURN
      END IF
      USEALL = .TRUE.
      IQAL = 0
      IF ( LENGTH .GT. 1 ) THEN    ! eg '*(10)' with appended qualifier string
        CALL REAQAL ( FIELD(2:LENGTH), IQAL, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END IF
!
    ELSE
      CALL GASCHK ( FIELD(1:LENGTH), NOTGAS, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      IF ( NOTGAS ) THEN
        LENGTH = MIN ( LENGTH, 15 ) 
        FAIL = .TRUE.
        ERRMSG = 'F-DRVGAS: *GAS section contains unrecognised ' // &
                 'molecule or file: '//FIELD(1:LENGTH)
        RETURN
      END IF
    END IF
  END DO
!
! Add any gases included by wildcard 
  IF ( USEALL ) THEN
    CALL GASALL ( IQAL, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
  END IF
!
! If using Rayleigh Extinction, add 'air' as an absorber 
  IF ( REXFLG ) CALL ADDGAS ( IDXAIR )
!
! Check at least one absorber specified
  IF ( NGAS .EQ. 0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-DRVGAS: No gases specified' 
    RETURN
  END IF
!
! Construct message for log file listing all identified species
  CALL WRTLOG ( 'I-DRVGAS: Using gases:', .TRUE. )
  DO IGAS = 1, NGAS
    CALL WRTLOG ( ' '//NAMGAS(IGAS), .TRUE. )
  END DO
  CALL WRTLOG ( '', .FALSE. )        ! Output all gases so far
!
! Check if CO2 listed in case of MIX or CHI flags
  CALL CHKCO2
!
END SUBROUTINE DRVGAS
END MODULE DRVGAS_SUB

