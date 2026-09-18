MODULE COOLAY_SUB
CONTAINS
SUBROUTINE COOLAY ( LATM, FLXRAD, ITNOFF ) 
!
! VERSION
!   06JUL26 AD Checked.
!   01AUG25 AD Original. 
!
! DESCRIPTION
!   Cooling rate calculation for atmospheric level
!   Called by FLXATM if COO flag enabled.
!   This calculates the contribution from a particular 
!   atmospheric level to all cooling rate output levels.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FULCOM_DAT ! Full grid data
    USE ATMCOM_DAT ! Atmospheric profile data
    USE TANCOM_DAT ! Tangent path data
!
! SUBROUTINES
    USE COOWGT_FNC ! Calculate cooling rate weights
    USE IDXARR_GEN ! Index of value within an array, else 0
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: LATM      ! Atmos level for cooling rate
    REAL(R8),    INTENT(IN) :: FLXRAD(:) ! Radiance flux
    INTEGER(I4), INTENT(IN) :: ITNOFF    ! Offset in TAN arrays
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: MAXWGT = 3 !  Max.no o/p levels affected by atm#
!                                           Quadratic fit implies a maximum=3
! LOCAL VARIABLES
    LOGICAL     :: FIRST = .TRUE. ! T=first call
    INTEGER(I4) :: IATM  ! Index of atm level below JATM
    INTEGER(I4) :: ITAN  ! Counter for output levels
    INTEGER(I4) :: IWGT  ! Counter for cooling rate weights
    INTEGER(I4) :: JATM  ! Index of mominal atm level
    INTEGER(I4) :: JTAN  ! Counter for cooling rate dependencies
    INTEGER(I4) :: KATM  ! Index of atm level above JATM
    REAL(R8),    ALLOCATABLE :: WGTATM(:,:) ! Weights for cooling rates
    SAVE FIRST, WGTATM
!       
! EXECUTABLE CODE -------------------------------------------------------------
!
! Calculate weights for each atm level on first call and save as WGTATM
  IF ( FIRST ) THEN
    ALLOCATE ( WGTATM(MAXWGT,NATM) ) 
    DO JATM = 1, NATM
      IATM = JATM - 1
      KATM = JATM + 1
      IF ( JATM .EQ. 1 ) IATM = 3
      IF ( JATM .EQ. NATM ) KATM = NATM - 2
      WGTATM(:,JATM) = &
        COOWGT ( (/ HGTATM(IATM), HGTATM(JATM), HGTATM(KATM) /), DNSATM(JATM) ) 
    END DO
    FIRST = .FALSE.
  END IF
!
! latm is the atm level for which flxrad has been calculated
! iatm is atm level influenced by flxrad at latm
  DO IWGT = 1, 3
    IATM = LATM + 2 - IWGT    ! latm+1, latm, latm-1, unless edges
    IF ( IATM .EQ. 0 ) IATM = 3
    IF ( IATM .EQ. NATM+1 ) IATM = NATM - 2
    ITAN = IDXARR ( TAN(1:NTAN)%IAT, IATM )
    IF ( ITAN .GT. 0 ) THEN
      JTAN = ITAN + ITNOFF
      COOFUL(IFUL1:IFUL2,JTAN) = COOFUL(IFUL1:IFUL2,JTAN) + & 
                                 FLXRAD * WGTATM(IWGT,IATM)
    END IF
  END DO
!
END SUBROUTINE COOLAY
END MODULE COOLAY_SUB
