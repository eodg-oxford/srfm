MODULE BRIGHT_FNC
CONTAINS
PURE FUNCTION BRIGHT ( RAD, WNOLST )
!
! VERSION
!   01MAR23 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Brightness Temperature calculation
!   General purpose module.
!   Inverse of the Planck function.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE PHYCON_DAT, ONLY: C1,C2 ! Radiation constants
!
  IMPLICIT NONE
!
! ARGUMENTS 
    REAL(R8), INTENT(IN) :: RAD(:)    ! Radiance [ nW/(cm2.sr.cm-1) ]
    REAL(R8), INTENT(IN) :: WNOLST(:) ! List of wavenumbers [cm-1]
!
! FUNCTION TYPE
    REAL(R8) :: BRIGHT ( SIZE(WNOLST) ) ! Function array same size as WNOLST
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  WHERE ( RAD .GT. TINY(1.0D0) ) 
    BRIGHT = C2 * WNOLST / LOG ( 1.0D0 + C1 * WNOLST**3 / RAD ) 
  ELSEWHERE
    BRIGHT = 0.0D0
  ENDWHERE
!
END FUNCTION BRIGHT
END MODULE BRIGHT_FNC
