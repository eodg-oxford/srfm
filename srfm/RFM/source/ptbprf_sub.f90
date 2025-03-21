MODULE PTBPRF_SUB
CONTAINS
SUBROUTINE PTBPRF ( ILOW, IPTB, IUPP, PTB, LOGPTB, HGT, PRF ) 
!
! VERSION
!   04DEC23 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Apply perturbation to a profile
!   Called by PTBATM.
!   If ILOW,IUPP within range of PRF index values, 
!      adds triangular perturbation (0:PTB:0) at PRF levels (ILOW:IPTB:IUPP).
!   If ILOW=0
!      adds flat perturbation (PTB:PTB) to PRF levels (1:IPTB)
!   If IUPP=NLEV+1 where NLEV=SIZE(PRF)
!      adds flat perturbation (PTB:PTB) to PRF levels (IPTB:NLEV)
!   If ILOW=0 and IUPP=NLEV+1
!      adds flat perturbation (PTB) to the whole profile
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN)    :: ILOW   ! Prf level# just below ptbd levels
    INTEGER(I4), INTENT(IN)    :: IPTB   ! Prf level# of max perturbation
    INTEGER(I4), INTENT(IN)    :: IUPP   ! Prf level# just above ptbd levels
    REAL(R4),    INTENT(IN)    :: PTB    ! Size of maximum perturbation
    LOGICAL,     INTENT(IN)    :: LOGPTB ! T=fract.ptb, F=linear ptb
    REAL(R4),    INTENT(IN)    :: HGT(:) ! Height profile [km]
    REAL(R4),    INTENT(INOUT) :: PRF(:) ! Profile
!
! LOCAL VARIABLES
    LOGICAL     :: LOWCOL ! True=all lower levels treated as column perturbation
    LOGICAL     :: UPPCOL ! True=all upper levels treated as column perturbation
    INTEGER(I4) :: IATM   ! Atmospheric level counter
    REAL(R4)    :: FACTOR ! Scaling factor for perturbation at atmos.level
!
! EXECUTABLE CODE ------------------------------------------------------------
!
  IF ( ILOW .LT. 0 ) STOP 'F-PTBPRF: Logical error#1'
  IF ( IUPP .GT. SIZE(PRF)+1 ) STOP 'F-PTBPRF: Logical error#2'
  IF ( ILOW .GT. IUPP-2 ) STOP 'F-PTBPRF: Logical error#3'
!
  LOWCOL = ILOW .EQ. 0
  UPPCOL = IUPP .EQ. SIZE(PRF)+1
!
  DO IATM = ILOW+1, IUPP-1
    IF ( IATM .LT. IPTB .AND. .NOT. LOWCOL ) THEN
      FACTOR = ( HGT(IATM) - HGT(ILOW) ) / ( HGT(IPTB) - HGT(ILOW) )
    ELSE IF ( IATM .GT. IPTB .AND. .NOT. UPPCOL ) THEN
      FACTOR = ( HGT(IUPP) - HGT(IATM) ) / ( HGT(IUPP) - HGT(IPTB) )
    ELSE
      FACTOR = 1.0
    END IF
    IF ( LOGPTB ) THEN
      PRF(IATM) = PRF(IATM) * (1.0 + FACTOR * PTB )
    ELSE
      PRF(IATM) = PRF(IATM) + FACTOR * PTB
    END IF
  END DO
!
END SUBROUTINE PTBPRF
END MODULE PTBPRF_SUB

