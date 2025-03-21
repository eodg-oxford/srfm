MODULE PLANCK_FNC
CONTAINS
PURE FUNCTION PLANCK ( TEM, WNOLST )
!
! VERSION
!   16JUL24 AD Checked.
!   23MAY20 AD Check for small Wavenumber limit - WNOMIN. Checked.
!   01MAY17 AD F90 conversion of F77 subroutine. Checked.
!
! DESCRIPTION
!   Planck Function
!   General purpose module.
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
    REAL(R4), INTENT(IN) :: TEM       ! Temperature [K]
    REAL(R8), INTENT(IN) :: WNOLST(:) ! List of wavenumbers [cm-1]
!
! FUNCTION TYPE
    REAL(R8) :: PLANCK ( SIZE(WNOLST) ) ! Function array same size as WNOLST
!
! LOCAL CONSTANTS
    REAL(R8), PARAMETER :: WNOMIN = 0.001D0  ! Min Wno for evaluation
!
! LOCAL VARIABLES
    INTEGER(I4) :: NLOW ! No. spectral points .LE. WNOMIN
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( TEM .GT. TINY(1.0) .AND. WNOLST(1) .GT. WNOMIN ) THEN
    PLANCK = C1 * WNOLST**3 / ( EXP ( C2 * WNOLST / TEM ) - 1.0D0 ) 
  ELSE IF ( TEM .LE. TINY(1.0) ) THEN
    PLANCK = 0.0D0
  ELSE
    NLOW = COUNT ( WNOLST .LE. WNOMIN )
    PLANCK(1:NLOW)  = C1 * WNOLST(1:NLOW)**2 * TEM / C2
    PLANCK(NLOW+1:) = C1 * WNOLST(NLOW+1:)**3 / &
                     ( EXP ( C2 * WNOLST(NLOW+1:) / TEM ) - 1.0D0 ) 
  END IF
!
END FUNCTION PLANCK
END MODULE PLANCK_FNC
