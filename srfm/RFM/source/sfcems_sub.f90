MODULE SFCEMS_SUB
CONTAINS
SUBROUTINE SFCEMS ( EMSSTR, FAIL, ERRMSG )
!
! VERSION
!   22DEC23 AD Checked.
!   21JUN17 AD F90 adpated from sfcfil.for. Checked.
!
! DESCRIPTION
!   Read surface emissivity data
!   Called by DRVSFC if any surface emissivity data in *SFC section.
!   The emissivity data in EMSSTR can either be the name of a spectral file
!   containing the emissivity spectrum, or just a single value.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE SFCCOM_DAT ! Surface parameters
    USE SPCCOM_DAT, ONLY: WMNSPC, WMXSPC ! Min/Max spec.range [cm-1]
    USE RFMLUN_DAT, ONLY: LUNTMP ! LUN for temporarily open files
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
    USE LEXIST_FNC ! Check if file exists
    USE OPNFIL_SUB ! Open input file
    USE SGNARR_GEN ! Determine if array is ascending or descending
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE 
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: EMSSTR ! Name of SFC data file or value string
    LOGICAL,       INTENT(OUT) :: FAIL   ! T=A fatal error was detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES      
    INTEGER(I4) :: IOS    ! Saved value of IOSTAT for error message
    INTEGER(I4) :: ISFC   ! Counter for SFC tabulation points
    REAL(R8)    :: DWSFC  ! Wno [cm-1] increment for emissivity values
    REAL(R8)    :: W1SFC  ! Wno [cm-1] of 1st emissivity value in file
    REAL(R8)    :: W2SFC  ! Wno [cm-1] of last emissivity value in file
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  IF ( LEXIST ( EMSSTR ) ) THEN                    ! file of emissivity values
    CALL OPNFIL ( LUNTMP, EMSSTR, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
!      
    READ ( LUNTMP, *, IOSTAT=IOS, ERR=900 ) NSFC, W1SFC, DWSFC, W2SFC
    IF ( NSFC .LE. 0 ) THEN
      FAIL = .TRUE.
      WRITE ( ERRMSG, * ) 'F-SFCEMS: No.tabulated SFC pts .le. 0, =', NSFC
      RETURN
    END IF
!
    ALLOCATE ( WNOSFC(NSFC), EMSSFC(NSFC) )
!
    IF ( DWSFC .EQ. 0.0D0 ) THEN        ! irregularly gridded data
      READ ( LUNTMP, *, IOSTAT=IOS, END=900, ERR=900 ) &
           (  WNOSFC(ISFC), EMSSFC(ISFC), ISFC = 1, NSFC ) 
      IF ( NSFC .GT. 1 ) THEN           ! check irregular grid is monotonic
        IF ( SGNARR ( WNOSFC ) .NE. 1 ) THEN 
          ERRMSG = 'F-SFCEMS: Emissivity spectral grid is not monotonic'
          FAIL = .TRUE.
          RETURN
        END IF
      END IF
    ELSE                                ! regularly gridded data
      READ ( LUNTMP, *, IOSTAT=IOS, END=900, ERR=900 ) EMSSFC
      DO ISFC = 1, NSFC 
        WNOSFC(ISFC) = W1SFC + ( ISFC - 1 ) * DWSFC 
      END DO
    END IF
!
    CLOSE ( LUNTMP, IOSTAT=IOS, ERR=900 ) 
!
    RFLSFC = .TRUE.
  ELSE    ! read value directly from field
    NSFC = 1
    ALLOCATE ( EMSSFC(NSFC), WNOSFC(NSFC) )
    WNOSFC(1) = 0.0D0
    READ ( EMSSTR, *, IOSTAT=IOS, ERR=900 ) EMSSFC(1)
    CALL WRTLOG ( 'I-SFCEMS: Setting Surface Emissivity = ' // &
                  C9REAL ( EMSSFC(1) ) )
    RFLSFC = EMSSFC(1) .LT. 1.0D0
  END IF    

! Check valid SFC values
  IF ( MINVAL ( EMSSFC ) .LT. 0.0D0 .OR. MAXVAL ( EMSSFC ) .GT. 1.0D0 ) THEN
    FAIL = .TRUE.
    ERRMSG = 'F-SFCEMS: Surface Emissivity outside range 0:1'
    RETURN
  END IF
!
! Check emissivity values span spectral range of all microwindows
  IF ( NSFC .GT. 1 .AND. &
       ( WMNSPC .LT. MINVAL ( WNOSFC ) .OR. WMXSPC .GT. MAXVAL ( WNOSFC ) ) ) &
    CALL WRTLOG ( 'W-SFCEMS: Emissivity not specified over full spectral range')
!
900 CONTINUE
  FAIL = IOS .NE. 0
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-SFCEMS: I/O failure reading SFC emissivity. IOSTAT=', IOS
!
END SUBROUTINE SFCEMS
END MODULE SFCEMS_SUB
