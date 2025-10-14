MODULE HITREQ_SUB
CONTAINS
SUBROUTINE HITREQ ( WNOREQ, IDMREQ )
!
! VERSION
!   15AUG24 AD Checked.
!   11AUG23 AD Original.
!
! DESCRIPTION
!   Set Wno range and list of molecules for HITRAN line data
!   Called once by DRVHIT.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GASCOM_DAT, ONLY: IGSMOL        ! IGAS for each molec#
    USE RFMCON_DAT, ONLY: FWIND         ! Window [cm-1] for widemesh calc
    USE SPCCOM_DAT, ONLY: WMNSPC,WMXSPC ! Upper Wavenumber reqd for any range
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(OUT) :: WNOREQ(2) ! Lower,Upper wavenumber [cm-1]
    LOGICAL,  INTENT(OUT) :: IDMREQ(:) ! T=use this line molecule
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IDMREQ = .FALSE.
!
! Establish max spectral range required to span all RFM calcs
  WNOREQ(1) = MAX ( WMNSPC - FWIND, 1.0D0 ) 
  WNOREQ(2) = WMXSPC + FWIND
!
! Flag which molecules are required 
  IDMREQ = IGSMOL(1:SIZE(IDMREQ)) .GT. 0 
!
END SUBROUTINE HITREQ
END MODULE HITREQ_SUB
