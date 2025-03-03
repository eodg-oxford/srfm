MODULE INTGRA_SUB
CONTAINS
SUBROUTINE INTGRA ( PSI, FLG, LININT, PRF )
!
! VERSION
!   07SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Interpolate to fill 2-D atmospheric profile field
!   Called by ATMGRA.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE INTERP_GEN ! Interpolate array
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN)    :: PSI(:)   ! Horizontal angle of each location
    LOGICAL,  INTENT(IN)    :: FLG(:)   ! T =location already filled
    LOGICAL,  INTENT(IN)    :: LININT   ! T = linear interpolation, F=log
    REAL(R4), INTENT(INOUT) :: PRF(:,:) ! Profiles
!
! LOCAL VARIABLES
    INTEGER(I4) :: IATM ! Profile level counter
    INTEGER(I4) :: IFIX ! Counter for fixed points
    INTEGER(I4) :: IPSI ! Index of horizontal locations
    INTEGER(I4) :: NFIX ! No. fixed points
    INTEGER(I4), ALLOCATABLE :: IDX(:)  ! Indices of fixed points
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( ALL ( FLG ) ) RETURN   ! All points already set
!
  NFIX = COUNT ( FLG ) 
! Should always be at least one set profile location
  IF ( NFIX .EQ. 0 ) STOP 'F-INTGRA: Logical error' 
!
  ALLOCATE ( IDX(NFIX) )
  IFIX = 0
  DO IPSI = 1, SIZE(FLG)         ! 1, NPSI
    IF ( FLG(IPSI) ) THEN
      IFIX = IFIX + 1
      IDX(IFIX) = IPSI 
    END IF
  END DO
!
  DO IATM = 1, SIZE(PRF,1)       ! 1,NATM 
    PRF(IATM,:) = INTERP ( PSI(IDX), PSI, PRF(IATM,IDX), .NOT. LININT  )     
  END DO
!
END SUBROUTINE INTGRA
END MODULE INTGRA_SUB

