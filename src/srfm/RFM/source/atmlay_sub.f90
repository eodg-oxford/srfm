MODULE ATMLAY_SUB
CONTAINS
SUBROUTINE ATMLAY ( GRDLEV, PREGRD ) 
!
! VERSION
!   21DEC23 AD Checked.
!   21JUN17 AD Allow for pressure grid as well as height grid. Checked.
!   01MAY17 AD F90 conversion. Tested.
!
! DESCRIPTION
!   Subdivide atmospheric profile layers
!   Called once by ATMGRD if LAY flag enabled.
!   Adds extra altitude grid points half-way between user-supplied profile.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! SUBROUTINES
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), ALLOCATABLE, INTENT(INOUT) :: GRDLEV(:) ! Profile altitudes [km]
    LOGICAL,               INTENT(IN)    :: PREGRD    ! T=p grid, F=z grid
!
! LOCAL VARIABLES
    INTEGER(I4)   :: INEW    ! New profile level
    INTEGER(I4)   :: IORG    ! Original profile level
    INTEGER(I4)   :: NNEW    ! No. of new levels in profile
    INTEGER(I4)   :: NORG    ! No. of original levels
    CHARACTER(80) :: MESSGE  ! Text message sent to Log file
    REAL(R4), ALLOCATABLE :: GRDORG(:) ! Original altitudes
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NORG = SIZE ( GRDLEV )
  NNEW = 2 * NORG - 1  ! Calculate the number of new profile levels
  MESSGE = 'I-ATMLAY: Number of atmos.profile levels increased from ' // &
           TRIM(C11INT(NORG)) // ' to ' // C11INT(NNEW) 
  CALL WRTLOG ( MESSGE )
!
  CALL MOVE_ALLOC ( GRDLEV, GRDORG ) 
  ALLOCATE ( GRDLEV(NNEW) )
!
  INEW = 1
  DO IORG = 1, NORG - 1
    GRDLEV(INEW) = GRDORG(IORG)
    IF ( PREGRD ) THEN 
      GRDLEV(INEW+1) = SQRT ( GRDORG(IORG) * GRDORG(IORG+1) ) 
    ELSE
      GRDLEV(INEW+1) = 0.5 * ( GRDORG(IORG) + GRDORG(IORG+1) ) 
    END IF
    INEW = INEW + 2
  END DO
  GRDLEV(NNEW) = GRDORG(NORG)
!
END SUBROUTINE ATMLAY
END MODULE ATMLAY_SUB

