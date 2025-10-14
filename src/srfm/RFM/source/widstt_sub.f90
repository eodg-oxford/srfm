MODULE WIDSTT_SUB 
CONTAINS
SUBROUTINE WIDSTT ( WNO, IGAS, JWIDL, JEXCL, JEXCU, JWIDU ) 
!
! VERSION
!   31OCT23 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION  
!   Add to line statistics 
!   Called by SPCWID for each line if WID flag enabled.
!
! VARIABLE KINDS
    USE KIND_DAT  
!
! GLOBAL DATA
    USE STTCOM_DAT ! Widemesh statistics
!
! SUBROUTINES
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8),    INTENT(IN) :: WNO   ! Wavenumber [cm-1] of line
    INTEGER(I4), INTENT(IN) :: IGAS  ! Index of absorber
    INTEGER(I4), INTENT(IN) :: JWIDL ! Lower WM Intvl for inclusion
    INTEGER(I4), INTENT(IN) :: JEXCL ! Lower WM Intvl for inclusion as non-local
    INTEGER(I4), INTENT(IN) :: JEXCU ! Upper WM Intvl for inclusion as non-local
    INTEGER(I4), INTENT(IN) :: JWIDU ! Upper WM Intvl for inclusion
!
! LOCAL VARIABLES
    INTEGER(I4) :: JWID ! WM interval containing line
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Although JWIDL,JWIDU limited to 1:NWID, JEXCL<1 and JEXCU>NWID are possible
  DO JWID = JWIDL, JWIDU
    IF ( JWID .LE. JEXCL .OR. JWID .GE. JEXCU ) THEN
      STT(JWID,IGAS)%NLW = STT(JWID,IGAS)%NLW + 1
    ELSE
      STT(JWID,IGAS)%NLF = STT(JWID,IGAS)%NLF + 1
    END IF
  END DO 
!
  JWID = FLOOR ( ( WNO - WN1STT ) / WNDSTT ) + 1
  IF ( JWID .GE. 1 .AND. JWID .LE. NSTW ) &
    STT(JWID,IGAS)%NLI = STT(JWID,IGAS)%NLI + 1
!
  IF ( WNO .LT. WNLSTT .OR. WNO .GT. WNUSTT ) THEN
    NEXSTT(IGAS) = NEXSTT(IGAS) + 1
  ELSE
    NINSTT(IGAS) = NINSTT(IGAS) + 1
  END IF
!
END SUBROUTINE WIDSTT
END MODULE WIDSTT_SUB
