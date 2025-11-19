MODULE IBRAKT_GEN
!
! VERSION
!   18JAN24 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Lower index of array interpolation
!   General purpose function
!   Return lower index I of the two points in array ARRAY(N)
!   whose values bracket X.  
!   i.e. L:  ARRAY(I)>X>ARRAY(I+1)  or  ARRAY(I)<X<ARRAY(I+1)
!   Special cases:
!       If X<ARRAY(1)<ARRAY(N)  or  X>ARRAY(1)>ARRAY(N), then I = 0
!       If X<ARRAY(N)<ARRAY(1)  or  X>ARRAY(N)>ARRAY(1), then I = N
!   Routine does not check parameters or if the array is monotonic
!
! VARIABLE KINDS
    USE KIND_DAT 
!
INTERFACE IBRAKT
  MODULE PROCEDURE IBRAKT_R, IBRAKT_D
END INTERFACE
!
CONTAINS
!
INTEGER(I4) PURE FUNCTION IBRAKT_R ( ARRAY, X, IGUESS, LIMIT )
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4),          INTENT(IN) :: ARRAY(:) ! Data Array
    REAL(R4),          INTENT(IN) :: X        ! Value to be bracketted
    INTEGER(I4), OPTIONAL, &
                       INTENT(IN) :: IGUESS   ! Guess of IBRAKT element below X
    LOGICAL, OPTIONAL, INTENT(IN) :: LIMIT    ! T=limit output to 1:N-1
!
! LOCAL VARIABLES
    LOGICAL     :: ASCEND ! TRUE=ascending array, FALSE=descending (assumed)
    INTEGER(I4) :: ISTEP  ! Step size, in array indices
    INTEGER(I4) :: I      ! Ongoing estimate of IBRAKT
    INTEGER(I4) :: J,M    ! High, Mid indices
    INTEGER(I4) :: N      ! Size of ARRAY
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  N = SIZE ( ARRAY )
!
! Ensure I is well-defined if N=1
  IF ( N .EQ. 1 ) THEN
    IBRAKT_R = 0
    RETURN
  END IF
!
! Determine if ascending or descending array (based on two end values only)
  ASCEND = ARRAY(1) .LT. ARRAY(N)
!
! If initial guess outside array index range, set to mid-point
  IF ( PRESENT ( IGUESS ) ) THEN
    I = IGUESS
    IF ( I .LE. 0 .OR. I .GT. N ) I = N / 2
  ELSE
    I = N / 2
  END IF
!
! Starting with an initial step size of one, find pair of indices I,J
! where J = I + ISTEP, which bracket the required X value I < X < J (asc.array)
  ISTEP = 1
  J = I
  IF ( X .GE. ARRAY(I) .EQV. ASCEND ) THEN     ! I < X in ascend.array
    DO WHILE ( X .GE. ARRAY(J) .EQV. ASCEND )  ! I < J < X
      I = J  
      J = I + ISTEP
      ISTEP = ISTEP * 2
      IF ( J .GT. N ) THEN
        J = N + 1
        EXIT
      END IF
    END DO
  ELSE                                          ! X < I in ascend.array 
    DO WHILE ( X .LT. ARRAY(I) .EQV. ASCEND )   ! X < I < J
      J = I
      I = J - ISTEP 
      ISTEP = ISTEP * 2
      IF ( I .LT. 1 ) THEN
        I = 0
        EXIT
      END IF
    END DO
  END IF
!
! Then progressively divide interval by 2 until J=I+1 with I .LE. X .LE. J 
  DO WHILE ( J - I .GT. 1 ) 
    M =  ( J + I ) / 2 
    IF ( X .GT. ARRAY(M) .EQV. ASCEND ) THEN
      I = M
    ELSE
      J = M
    END IF
  END DO
!
  IF ( PRESENT ( LIMIT ) ) THEN
    IF ( LIMIT ) THEN
      IF ( I .EQ. 0 ) I = 1
      IF ( I .EQ. N ) I = N-1
    END IF
  END IF
!
  IBRAKT_R = I
!
END FUNCTION IBRAKT_R

INTEGER(I4) PURE FUNCTION IBRAKT_D ( ARRAY, X, IGUESS, LIMIT )
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8),          INTENT(IN) :: ARRAY(:) ! Data Array
    REAL(R8),          INTENT(IN) :: X        ! Value to be bracketted
    INTEGER(I4), OPTIONAL, &
                       INTENT(IN) :: IGUESS   ! Guess of IBRAKT element below X
    LOGICAL, OPTIONAL, INTENT(IN) :: LIMIT    ! T=limit output to 1:N-1
!
! LOCAL VARIABLES
    LOGICAL     :: ASCEND ! TRUE=ascending array, FALSE=descending (assumed)
    INTEGER(I4) :: ISTEP  ! Step size, in array indices
    INTEGER(I4) :: I      ! Ongoing estimate of IBRAKT
    INTEGER(I4) :: J,M    ! High, Mid indices
    INTEGER(I4) :: N      ! Size of ARRAY
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  N = SIZE ( ARRAY )
!
! Ensure I is well-defined if N=1
  IF ( N .EQ. 1 ) THEN
    IBRAKT_D = 0
    RETURN
  END IF
!
! Determine if ascending or descending array (based on two end values only)
  ASCEND = ARRAY(1) .LT. ARRAY(N)
!
! If initial guess outside array index range, set to mid-point
  IF ( PRESENT ( IGUESS ) ) THEN 
    I = IGUESS
    IF ( I .LE. 0 .OR. I .GT. N ) I = N / 2
  ELSE
    I = N / 2
  END IF
!
! Starting with an initial step size of one, find pair of indices I,J
! where J = I + ISTEP, which bracket the required X value I < X < J (asc.array)
  ISTEP = 1
  J = I
  IF ( X .GE. ARRAY(I) .EQV. ASCEND ) THEN     ! I < X in ascend.array
    DO WHILE ( X .GE. ARRAY(J) .EQV. ASCEND )  ! I < J < X
      I = J  
      J = I + ISTEP
      ISTEP = ISTEP * 2
      IF ( J .GT. N ) THEN
        J = N + 1
        EXIT
      END IF
    END DO
  ELSE                                          ! X < I in ascend.array 
    DO WHILE ( X .LT. ARRAY(I) .EQV. ASCEND )   ! X < I < J
      J = I
      I = J - ISTEP 
      ISTEP = ISTEP * 2
      IF ( I .LT. 1 ) THEN
        I = 0
        EXIT
      END IF
    END DO
  END IF
!
! Then progressively divide interval by 2 until J=I+1 with I .LE. X .LE. J 
  DO WHILE ( J - I .GT. 1 ) 
    M =  ( J + I ) / 2 
    IF ( X .GT. ARRAY(M) .EQV. ASCEND ) THEN
      I = M
    ELSE
      J = M
    END IF
  END DO
!
  IF ( PRESENT ( LIMIT ) ) THEN
    IF ( LIMIT ) THEN
      IF ( I .EQ. 0 ) I = 1
      IF ( I .EQ. N ) I = N-1
    END IF
  END IF
!
  IBRAKT_D = I
!
END FUNCTION IBRAKT_D

END MODULE IBRAKT_GEN
