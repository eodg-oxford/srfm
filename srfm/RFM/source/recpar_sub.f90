MODULE RECPAR_SUB
CONTAINS
SUBROUTINE RECPAR ( LUNHIT, HIT, USEIDM, WNOMAX, EOF, FAIL, ERRMSG ) 
!
! VERSION
!   24AUG24 AD Checked.
!   11AUG23 AD Pass data via arguments.
!   31MAY23 AD Original. Simplified from HITREC
!
! DESCRIPTION
!   Read next record from HITRAN line data .par file
!   Called by HITREC, INIHIT.
!   Reads through successive records until one of the required molecules,
!   or WNOMAX or EOF is reached.
!   EOF is also set TRUE if WNO > WNOMAX.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE HITCOM_DAT, ONLY: HITTYP ! HITRAN line data structure
    USE IDXCON_DAT, ONLY: IDXH2O ! RFM/HITRAN index for H2O
    USE PHYCON_DAT, ONLY: AVOG   ! Avogradro's number [kmol/cm2]
!
  IMPLICIT NONE
!
! ARGUMENTS      
    INTEGER(I4),   INTENT(IN)  :: LUNHIT    ! LUN for HITRAN file
    TYPE(HITTYP),  INTENT(OUT) :: HIT       ! HITRAN Line parameters
    LOGICAL,       INTENT(IN)  :: USEIDM(:) ! T=use molecule#
    REAL(R8),      INTENT(IN)  :: WNOMAX    ! Highest wno [cm-1] reqd
    LOGICAL,       INTENT(OUT) :: EOF       ! Set TRUE if end-of-file reached
    LOGICAL,       INTENT(OUT) :: FAIL      ! TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG    ! Error message if FAIL is TRUE
!
! LOCAL CONSTANTS
    CHARACTER(*), PARAMETER :: HITFMT = &
      '( I2, Z1, F12.6, F10.3, E10.3, F5.2, F5.2, F10.4, F4.1, F8.5 )'
    INTEGER(I4), PARAMETER :: FIXIDI(0:15) = &
              (/ 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16 /)
!
! LOCAL VARIABLES
    INTEGER(I4) :: IDI    ! Isotope ID read from .par file
    INTEGER(I4) :: IOS    ! Saved value of IOSTAT for error messages
    REAL(R4)    :: TPROB  ! Transition probability [Debyes2] (dummy)
    REAL(R8)    :: DSTR   ! Original HITRAN linestrength
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
  EOF = .FALSE.
!
  DO         ! Continue reading until a record with a required molecule found
    READ ( LUNHIT, HITFMT, IOSTAT=IOS, ERR=900, END=800 ) HIT%IDM, IDI, &
          HIT%WNO, DSTR, TPROB, HIT%HWA, HIT%HWS, HIT%ELS, HIT%TCA, HIT%PSA
    IF ( HIT%WNO .GT. WNOMAX ) GOTO 800 ! no more records within required range
    IF ( USEIDM(HIT%IDM) ) THEN 
      HIT%STR = SNGL ( DSTR * AVOG ) ! Scale intensity by Avogadro
      HIT%IDI = FIXIDI(IDI)          ! Convert to correct numbering
! Older data may have self-broadened half-width values set to zero
      IF ( HIT%HWS .EQ. 0.0 ) THEN
        IF ( HIT%IDM .EQ. IDXH2O ) THEN
          HIT%HWS = 5.0 * HIT%HWA
        ELSE
          HIT%HWS = HIT%HWA
        END IF 
      END IF
      RETURN                         ! found record with a required molecule
    END IF
  END DO
!
! Either reached end-of-file or beyond highest required wavenumber
800 CONTINUE
  EOF = .TRUE.
  RETURN
!
900 CONTINUE
  FAIL = .TRUE.
  WRITE ( ERRMSG, * ) &
    'F-RECPAR: Failed to read record in HITRAN File. IOSTAT=', IOS
!
END SUBROUTINE RECPAR
END MODULE RECPAR_SUB
