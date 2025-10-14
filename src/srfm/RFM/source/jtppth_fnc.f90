MODULE JTPPTH_FNC
CONTAINS
LOGICAL FUNCTION JTPPTH ( IPTH, IJAC ) 
!
! VERSION
!   28NOV23 AD Checked.
!   31MAY18 AD Ensure JTAN=0 is handled. Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   T=path corresponds to Jacobian Tan.Pt perturbation
!   Called by JACPTH
!   Set TRUE if tangent point of ray for this path segment lies within the
!   perturbation range of the Jacobian. 
!
! GLOBAL DATA
    USE JACCOM_DAT ! Jacobian data
    USE PTHCOM_DAT ! Path segment data
    USE FOVCOM_DAT, ONLY: ITNFOV ! Indices of tan paths for FOV convolution
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IPTH ! Index of path segment
    INTEGER(I4), INTENT(IN) :: IJAC ! Index of Jacobian
!
! LOCAL VARIABLES
    INTEGER(I4) :: ITAN ! Index of tangent ray for PTH
    INTEGER(I4) :: JTAN ! Index of nominal tangent ray for Jacobian
!
! EXECUTABLE CODE ------------------------------------------------------------
!
  JTAN = JAC(IJAC)%ITN
  ITAN = PTH(IPTH)%ITN
  JTPPTH = JTAN .EQ. ITAN
  IF ( JTPPTH ) RETURN       ! Path along nominal tan.ray
  IF ( JTAN .EQ. 0 ) RETURN  ! If called with IJAC not a tan.pt Jacobian
!
! Also TRUE if IPTH lies along any tan.rays required to construct FOV for
! nominal tan.ray
  IF ( ALLOCATED(ITNFOV) ) JTPPTH =  ANY ( ITNFOV(JTAN,:) .EQ. ITAN ) 
!
END FUNCTION JTPPTH
END MODULE JTPPTH_FNC
