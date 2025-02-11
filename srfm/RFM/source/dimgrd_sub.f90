MODULE DIMGRD_SUB
CONTAINS
SUBROUTINE DIMGRD ( LUNDRV, AXTYP, IFLD, FIELD, NTAB, AXTAB, FAIL, ERRMSG ) 
!
! VERSION
!   20OCT23 AD Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Read TAB axis grid parameters from driver file
!   Called by DRVDIM for p, T and q-axis values.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE NXTFLD_SUB ! Load next field from section of driver file
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),    INTENT(IN) :: LUNDRV ! LUN for RFM driver file
    CHARACTER(1),   INTENT(IN) :: AXTYP  ! Axis type ('P', 'T' or 'Q')
    INTEGER(I4), INTENT(INOUT) :: IFLD   ! Counter for field# in *DIM section
    CHARACTER(*),INTENT(INOUT) :: FIELD  ! Current field read from *DIM section
    INTEGER(I4),   INTENT(OUT) :: NTAB   ! No elements of AXTAB written
    REAL(R4), ALLOCATABLE, &
                   INTENT(OUT) :: AXTAB(:) ! Array for axis values
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IOS    ! Saved value of IOSTAT for error message
    INTEGER(I4)   :: ITAB   ! Counter for axis points
    INTEGER(I4)   :: LENGTH ! Length of field read from driver table
    REAL(R4)      :: AX1    ! 1st axis value
    REAL(R4)      :: AX2    ! last axis value
    REAL(R4)      :: DX     ! Axis increment
    REAL(R4)      :: LOGAX1 ! Ln(AX1)
    REAL(R4)      :: RDUM   ! Real variable
    CHARACTER(12) :: FSTR   ! Field counter/contents for error message
!
! EXECUTABLE CODE --------------------------------------------------------------
!         123456789012
  FSTR = 'Field#x (xx)'
  WRITE ( FSTR(7:7), '(I1)' ) IFLD
  WRITE ( FSTR(10:11), '(2A1)' ) 'N', AXTYP
!
  CALL WRTLOG ( AXTYP // '-axis parameters: ', .TRUE. )
  READ ( FIELD, *, IOSTAT=IOS, ERR=900 ) RDUM
  NTAB = NINT ( RDUM )
  FAIL = .TRUE.
  IF ( RDUM - NTAB .NE. 0.0 ) THEN
    ERRMSG = 'F-DIMGRD: ' // FSTR // &
             ' was not an integer, actual value=' // FIELD
  ELSE IF ( NTAB .EQ. 0 .AND. AXTYP .EQ. 'P' ) THEN
    ERRMSG = 'F-DIMGRD: At least one value required for P-axis'
  ELSE IF ( NTAB .LE. 0 .AND. AXTYP .EQ. 'T' ) THEN
    ERRMSG = 'F-DIMGRD: NPts for T-axis must be +ve'
  ELSE
    FAIL = .FALSE.
  END IF
  IF ( FAIL ) RETURN
  CALL WRTLOG ( '', .FALSE. ) 
!
  IFLD = IFLD + 1
  WRITE ( FSTR(7:7), '(I1)' ) IFLD
  WRITE ( FSTR(10:11), '(2A1)' ) AXTYP, '1'
!
  CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
  IF ( LENGTH .EQ. 0 ) THEN
    ERRMSG = 'F-DIMGRD: No value for ' // FSTR
    FAIL = .TRUE.
    RETURN
  END IF
   READ ( FIELD, *, IOSTAT=IOS, ERR=900 ) AX1 
  CALL WRTLOG ( FIELD, .TRUE. )
  IF ( ABS ( NTAB ) .EQ. 1 ) THEN                       ! Single axis value
    AX2 = AX1
  ELSE 
!
    IFLD = IFLD + 1
    WRITE ( FSTR(7:7), '(I1)' ) IFLD
    WRITE ( FSTR(10:11), '(2A1)' ) AXTYP, '2'
!
    CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) THEN
      ERRMSG = 'F-DIMGRD: No value for ' // FSTR
      FAIL = .TRUE.
      RETURN
    END IF
    CALL WRTLOG ( FIELD, .TRUE. ) 
    READ ( FIELD, *, IOSTAT=IOS, ERR=900 ) AX2 
!
  END IF
  CALL WRTLOG ( '', .FALSE. ) 
!
  FAIL = .TRUE.
  IF ( MIN ( AX1, AX2 ) .LE. 0.0 .AND. &
       ( AXTYP .EQ. 'P' .OR. AXTYP .EQ. 'Q' ) ) THEN
    ERRMSG = 'F-DIMGRD: ' // AXTYP // '-axis values must all be +ve'
  ELSE IF ( MAX ( AX1, AX2 ) .LT. 0.0 .AND. AXTYP .EQ. 'T' ) THEN
    ERRMSG = 'F-DIMGRD: Some T-axis values must be 0 or +ve'
  ELSE
    FAIL = .FALSE.
  END IF
  IF ( FAIL ) RETURN
!
! Calculate axis values
  ALLOCATE ( AXTAB(ABS(NTAB)) )
  IF ( NTAB .EQ. 1 ) THEN
    AXTAB(1) = AX1
  ELSE IF ( NTAB .LT. 0 ) THEN
    NTAB = ABS ( NTAB )
    LOGAX1 = LOG ( AX1 )
    DX = ( LOG ( AX2 ) - LOGAX1 ) / ( NTAB - 1 )
    DO ITAB = 1, NTAB 
      AXTAB(ITAB) = EXP ( LOGAX1 + ( ITAB - 1 ) * DX ) 
    END DO
  ELSE
    DX = ( AX2 - AX1 ) / ( NTAB - 1 )
    DO ITAB = 1, NTAB
      AXTAB(ITAB) = AX1 + ( ITAB - 1 ) * DX 
    END DO
  END IF
!
 900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL )  WRITE ( ERRMSG, * ) & 
    'F-DIMGRD: I/O error reading ' // FSTR // '. IOSTAT=', IOS
!
END SUBROUTINE DIMGRD
END MODULE DIMGRD_SUB
