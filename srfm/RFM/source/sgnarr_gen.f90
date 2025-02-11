MODULE SGNARR_GEN
!
! VERSION
!   06SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Determine if array is ascending or descending
!   General purpose module
!   Returns value +1 if ascending, -1 if descending, 0= anything else, eg mixed
!
! VARIABLE KINDS
    USE KIND_DAT
!
INTERFACE SGNARR
  MODULE PROCEDURE SGNARR_R, SGNARR_D
END INTERFACE
!
CONTAINS
!
INTEGER(I4) PURE FUNCTION SGNARR_R ( X )
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: X(:) ! Array to be tested  
!
! LOCAL VARIABLES
    INTEGER(I4) :: ISIGN ! Local assignment of sign
    INTEGER(I4) :: IX    ! Array counter
    INTEGER(I4) :: NX    ! Size of array
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! Initially assume sign is indeterminate
  SGNARR_R = 0
!
  NX = SIZE(X)
!
! For scalar, return SGNARR=0
  IF ( NX .LE. 1 ) RETURN
!
! Determine sign from first 2 array elements
  IF ( X(2) .GT. X(1) ) THEN
    ISIGN = 1
  ELSE IF ( X(2) .LT. X(1) ) THEN
    ISIGN = -1
  ELSE     
    RETURN               ! If equal, return with SGNARR=0
  END IF
!
  DO IX = 3, NX
    IF ( ( ISIGN .EQ. 1 .AND. X(IX) .LE. X(IX-1) ) .OR. &
         ( ISIGN .EQ.-1 .AND. X(IX) .GE. X(IX-1) )         ) RETURN
  END DO
!
  SGNARR_R = ISIGN
!
END FUNCTION SGNARR_R

INTEGER(I4) PURE FUNCTION SGNARR_D ( X )
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(IN) :: X(:) ! Array to be tested  
!
! LOCAL VARIABLES
    INTEGER(I4) :: ISIGN ! Local assignment of sign
    INTEGER(I4) :: IX    ! Array counter
    INTEGER(I4) :: NX    ! Size of array
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! Initially assume sign is indeterminate
  SGNARR_D = 0
!
  NX = SIZE(X)
!
! For scalar, return SGNARR=0
  IF ( NX .LE. 1 ) RETURN
!
! Determine sign from first 2 array elements
  IF ( X(2) .GT. X(1) ) THEN
    ISIGN = 1
  ELSE IF ( X(2) .LT. X(1) ) THEN
    ISIGN = -1
  ELSE     
    RETURN               ! If equal, return with SGNARR=0
  END IF
!
  DO IX = 3, NX
    IF ( ( ISIGN .EQ. 1 .AND. X(IX) .LE. X(IX-1) ) .OR. &
         ( ISIGN .EQ.-1 .AND. X(IX) .GE. X(IX-1) )         ) RETURN
  END DO
!
  SGNARR_D = ISIGN
!
END FUNCTION SGNARR_D
END MODULE SGNARR_GEN
