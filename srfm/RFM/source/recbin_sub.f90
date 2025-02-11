MODULE RECBIN_SUB
CONTAINS
SUBROUTINE RECBIN ( LUNHIT, HIT, USEIDM, IFP, EOF, FAIL, ERRMSG ) 
!
! VERSION
!   23AUG24 AD Checked.
!   11AUG23 AD Add arguments
!   31MAY23 AD Simplified version of HITREC for binary file only.
!
! DESCRIPTION
!   Read record from HITRAN line data binary file
!   Called by INIHIT and HITREC.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE HITCOM_DAT, ONLY: HITTYP ! HITRAN line data structure
    USE IDXCON_DAT, ONLY: IDXH2O ! RFM/HITRAN index for H2O
!
! SUBROUTINES
    USE C11INT_FNC ! Write integer as left-adjusted string
!    
  IMPLICIT NONE
!
! ARGUMENTS      
    INTEGER(I4),   INTENT(IN)  :: LUNHIT    ! LUN for HITRAN file
    TYPE(HITTYP),  INTENT(OUT) :: HIT       ! HITRAN line data
    LOGICAL,       INTENT(IN)  :: USEIDM(:) ! T=use molecule#
    INTEGER(I4), INTENT(INOUT) :: IFP(:)    ! Forward pointers 
    LOGICAL,       INTENT(OUT) :: EOF       ! Set TRUE if end-of-file reached
    LOGICAL,       INTENT(OUT) :: FAIL      ! T= a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG    ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)  :: IFWDPT ! Forward pointer for molecule
    INTEGER(I4)  :: IOS    ! Saved value of IOSTAT for error messages
    INTEGER(I4)  :: IREC   ! Record# to be read
    INTEGER(I4)  :: LSTAT  ! Status of transition information.
    REAL(R4)     :: TPROB  ! Transition probability [Debyes2].
    CHARACTER(9) :: SPARE9 ! Spare bytes in binary file record
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
  IREC = MINVAL ( IFP, MASK=USEIDM )
  EOF = IREC .EQ. MAXVAL ( IFP )  ! Assume MAX IFP set to end-of-file record
  IF ( EOF ) RETURN

  READ ( LUNHIT, REC=IREC, IOSTAT=IOS, ERR=900 ) &     ! Load /HITCOM/
    LSTAT, HIT%IDM, HIT%IDI, HIT%WNO, HIT%STR,  TPROB, &
    HIT%HWA, HIT%HWS, HIT%ELS, HIT%TCA, HIT%PSA, HIT%IUS, HIT%ILS, &
    HIT%ULQ, HIT%BLQ, SPARE9, IFWDPT
! Older data may have self-broadened half-width values set to zero
  IF ( HIT%HWS .EQ. 0.0 ) THEN
    IF ( HIT%IDM .EQ. IDXH2O ) THEN
      HIT%HWS = 5.0 * HIT%HWA
    ELSE
      HIT%HWS = HIT%HWA
    END IF 
  END IF
  IFP(HIT%IDM) = IREC + IFWDPT
!
900 CONTINUE
  FAIL = IOS .NE. 0 
  IF ( FAIL ) WRITE ( ERRMSG, * ) 'F-RECBIN: Failed to read Rec#' // &
    TRIM(C11INT(IREC)) // ' in HITRAN File. IOSTAT=', IOS
!
END SUBROUTINE RECBIN
END MODULE RECBIN_SUB
