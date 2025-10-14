MODULE INIQAD_SUB
CONTAINS
SUBROUTINE INIQAD 
!
! VERSION
!   26APR24 AD Checked.
!   05MAR19 AD Original. Checked.
!
! DESCRIPTION
!   Initialise Gaussian quadrature for flux calculations
!   Called once by CHKSFC or TANFLX.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE QADCOM_DAT ! Gaussian quadrature data
    USE FLGCOM_DAT, ONLY: VRTFLG ! Option flags
    USE PHYCON_DAT, ONLY: PI
!
! SUBROUTINES
    USE GAUQAD_SUB ! Values and Weights for Gaussian First Moment Quadrature
!
  IMPLICIT NONE
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( VRTFLG ) THEN    
    NQAD = 1            ! alter default value of NQAD in QADCOM_DAT
    RPIQAD = 1.0D0
    ALLOCATE ( XQAD(1), WQAD(1) ) 
    XQAD(1) = 1.0D0
    WQAD(1) = 1.0D0
  ELSE                  ! use default value of NQAD in QADCOM_DAT
    RPIQAD = 1.0D0 / PI
    ALLOCATE ( XQAD(NQAD), WQAD(NQAD) ) 
    CALL GAUQAD ( NQAD, XQAD, WQAD ) 
    WQAD = WQAD * 2 * PI
  END IF
!
END SUBROUTINE INIQAD
END MODULE INIQAD_SUB
