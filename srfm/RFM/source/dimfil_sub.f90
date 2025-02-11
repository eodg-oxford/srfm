MODULE DIMFIL_SUB
CONTAINS
SUBROUTINE DIMFIL ( AXTYP, FILDIM, NTAB, AXTAB, FAIL, ERRMSG )
!
! VERSION
!   16NOV23 AD Checked.
!   01MAY17 AD F90 conversion. Tested
!
! DESCRIPTION
!   Read list of TAB axis values from file
!   Called by DRVDIM for p, T and q-axis files.
!   Also checks that at least one value is contained in file and that
!   array values are monotonic.
!
! VARIABLE KINDS
    USE KIND_DAT
! 
! GLOBAL DATA
    USE LENREC_DAT ! Max length of input text record
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE NXTFL2_SUB ! Load next field from secondary file
    USE OPNFIL_SUB ! Open input file
    USE SGNARR_GEN ! Determine if array is ascending or descending 
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(1),  INTENT(IN)  :: AXTYP  ! Tabulation axis 'P', 'T' or 'Q'
    CHARACTER(*),  INTENT(IN)  :: FILDIM ! Name of file containing axis values
    INTEGER(I4),  INTENT(OUT)  :: NTAB   ! No elements of AXTAB written
    REAL(R4), ALLOCATABLE, &
                   INTENT(OUT) :: AXTAB(:) ! Array for axis values
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)           :: IOS       ! Saved value of IOSTAT for error message
    INTEGER(I4)           :: LENGTH    ! Length of field read from driver table
    CHARACTER(LENREC)     :: FIELD     ! Tangent height or Tan.Hgt filename
    REAL(R4), ALLOCATABLE :: TABSAV(:) ! Saved value of AXTAB
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  CALL WRTLOG ( 'I-DIMFIL: Reading ' // AXTYP // '-axis values from file' )
!
! Open file containing axis values
  CALL OPNFIL ( LUNTMP, FILDIM, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
! Read in each field from file
  NTAB = 0
  LENGTH = 1
  CALL WRTLOG ( AXTYP//'-axis values:', .TRUE. )
  DO 
    CALL NXTFL2 ( LUNTMP, FIELD, LENGTH, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( LENGTH .EQ. 0 ) EXIT       ! LENGTH=0 denotes end of file
    NTAB = NTAB + 1 
    IF ( NTAB .EQ. 1 ) THEN
      ALLOCATE ( AXTAB(NTAB) ) 
    ELSE 
      CALL MOVE_ALLOC ( AXTAB, TABSAV ) 
      ALLOCATE ( AXTAB(NTAB) ) 
      AXTAB(1:NTAB-1) = TABSAV
    END IF
    READ ( FIELD, *, IOSTAT=IOS, ERR=900 ) AXTAB(NTAB) 
    CALL WRTLOG ( FIELD, .TRUE. ) 
    IF ( FAIL ) RETURN
  END DO
!
  CLOSE ( LUNTMP, IOSTAT=IOS, ERR=900 )
  CALL WRTLOG ( '', .FALSE. ) 
!
  FAIL = .TRUE.
  IF ( NTAB .EQ. 0 ) THEN
    ERRMSG = 'F-DIMFIL: No values found in file'
! NTAB=1 is valid but SGNARR will return 0 for this case, so exclude
  ELSE IF ( NTAB .GE. 2 .AND. SGNARR ( AXTAB ) .EQ. 0 ) THEN
    ERRMSG = 'F-DIMFIL: axis values are not monotonic'
! P,Q axes must all be positive numbers
  ELSE IF ( MIN ( AXTAB(1), AXTAB(NTAB) ) .LE. 0.0 .AND. &
            ( AXTYP .EQ. 'P' .OR. AXTYP .EQ. 'Q' )        ) THEN
    ERRMSG = 'F-DIMFIL: '//AXTYP//'-axis values must all be +ve'
  ELSE IF ( MAX ( AXTAB(1), AXTAB(NTAB) ) .LT. 0.0 .AND. AXTYP .EQ. 'T' ) THEN
    ERRMSG = 'F-DIMFIL: Some T-axis values must be 0 or +ve'
  ELSE
    FAIL = .FALSE.
  END IF
!
! Normal exit, may be with FAIL=TRUE from above tests
  RETURN
!
! Exit with fatal I/O error
900 CONTINUE
  FAIL = .TRUE.
  WRITE ( ERRMSG, * ) &
    'F-DIMFIL: I/O error on ' // AXTYP // '-axis file. IOSTAT=', IOS
!
END SUBROUTINE DIMFIL
END MODULE DIMFIL_SUB

