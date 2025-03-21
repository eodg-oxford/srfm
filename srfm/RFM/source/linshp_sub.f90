MODULE LINSHP_SUB
CONTAINS
SUBROUTINE LINSHP ( ISHP, WNOGRD, ABSLIN, SUBWNG )
!
! VERSION
!   29AUG24 AD Checked.
!   24JUN19 AD Allow for CIA 'lineshape'. Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Apply spectral lineshape
!   Called by SPCWID, SPCFIN for each line.
!   Assumes line parameters are stored in ADJCOM.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE SHPCON_DAT ! Line-shape indices
!
! SUBROUTINES
    USE CHISHP_SUB ! Calculate Voigt Line shape allowing for chi-factor
    USE DOPSHP_SUB ! Calculate Doppler Line shape
    USE LORSHP_SUB ! Calculate Lorentz Line shape
    USE MIXSHP_SUB ! Calculate Voigt Line shape allowing for line mixing
    USE SUB25W_SUB ! Subtract absorption value at 25cm-1 from line centre
    USE VOISHP_SUB ! Calculate Voigt Line shape
    USE VVWCOR_SUB ! Apply Van Vleck-Weisskopf correction to line shape
    USE VVWSHP_SUB ! Calculate Van Vleck-Weisskopf Line shape
!
  IMPLICIT NONE
!                  
! ARGUMENTS
    INTEGER(I4), INTENT(IN)  :: ISHP      ! Index of lineshape function
    REAL(R8),    INTENT(IN)  :: WNOGRD(:) ! Array of Wavenumbers [cm-1]
    REAL(R4),    INTENT(OUT) :: ABSLIN(:) ! Absorption
    LOGICAL,     INTENT(IN)  :: SUBWNG    ! T = Subtract abs.coeff at 25cm-1
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  SELECT CASE ( ISHP ) 
  CASE ( SHPVOI ) ; CALL VOISHP ( WNOGRD, ABSLIN ) 
  CASE ( SHPMIX ) ; CALL MIXSHP ( WNOGRD, ABSLIN ) 
  CASE ( SHPLOR ) ; CALL LORSHP ( WNOGRD, ABSLIN )
  CASE ( SHPDOP ) ; CALL DOPSHP ( WNOGRD, ABSLIN )
  CASE ( SHPCHI ) ; CALL CHISHP ( WNOGRD, ABSLIN )
  CASE ( SHPVVW ) ; CALL VVWSHP ( WNOGRD, ABSLIN )
  CASE ( SHPCIA ) ; ABSLIN = 0.0   ! CIA-only also allocated for LBL calc
  CASE ( SHPCTM ) ; ABSLIN = 0.0   ! CTM-only also allocated for LBL calc
  CASE DEFAULT    ; STOP 'F-LINSHP: Unrecognised lineshape value'
  END SELECT
!
  IF ( VVWFLG .AND. ISHP .NE. SHPVVW ) CALL VVWCOR ( WNOGRD, ABSLIN )
!
! If applying CKD H2O continuum need to subtract abs.coeff at 25cm
  IF ( SUBWNG ) CALL SUB25W ( ABSLIN )
!
END SUBROUTINE LINSHP
END MODULE LINSHP_SUB

