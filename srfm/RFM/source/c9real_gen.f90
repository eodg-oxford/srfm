MODULE C9REAL_GEN
!
! VERSION
!   11FEB24 AD Checked.
!   01MAY17 AD F90 conversion of c9flt.for. Checked.
!
! DESCRIPTION
!   Write real number as C*9 string
!   General purpose module
!   C*9 accommodates at least 3 sig figures.
!
! VARIABLE KINDS
    USE KIND_DAT
!
INTERFACE C9REAL
  MODULE PROCEDURE C9REAL_R, C9REAL_D
END INTERFACE

CONTAINS

CHARACTER(9) PURE FUNCTION C9REAL_R ( RNUM ) 
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: RNUM ! Number to be written
!
! LOCAL VARIABLES
    REAL(R4)     :: ANUM ! Absolute value of RNUM  
    CHARACTER(9) :: C9   ! Internal version of eventual output 
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  ANUM = ABS ( RNUM ) 
  IF ( ANUM .GT. 999.5 ) THEN
    WRITE ( C9, '(ES9.2)' ) RNUM
  ELSE IF ( ANUM .GT. 99.95 ) THEN
     WRITE ( C9, '(F9.0)' ) RNUM
   ELSE IF ( ANUM .GT. 9.995 ) THEN
    WRITE ( C9, '(F9.1)' ) RNUM
  ELSE IF ( ANUM .GT. 0.9995 ) THEN
    WRITE ( C9, '(F9.2)' ) RNUM
  ELSE IF ( ANUM .GT. 0.09995 ) THEN
    WRITE ( C9, '(F9.3)' ) RNUM
  ELSE IF ( ANUM .GT. 0.009995 ) THEN
    WRITE ( C9, '(F9.4)' ) RNUM
  ELSE 
    WRITE ( C9, '(ES9.2)' ) RNUM
  END IF 
!
  C9REAL_R = ADJUSTL ( C9 ) 
!
END FUNCTION C9REAL_R

CHARACTER(9) PURE FUNCTION C9REAL_D ( DNUM ) 
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(IN) :: DNUM ! Number to be written
!
! LOCAL VARIABLES
    REAL(R4)     :: ANUM ! Absolute value of DNUM  
    CHARACTER(9) :: C9   ! Internal version of eventual output 
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  ANUM = ABS ( SNGL ( DNUM ) ) 
  IF ( ANUM .GT. 999.5 ) THEN
    WRITE ( C9, '(ES9.2)' ) DNUM
  ELSE IF ( ANUM .GT. 99.95 ) THEN
     WRITE ( C9, '(F9.0)' ) DNUM
   ELSE IF ( ANUM .GT. 9.995 ) THEN
    WRITE ( C9, '(F9.1)' ) DNUM
  ELSE IF ( ANUM .GT. 0.9995 ) THEN
    WRITE ( C9, '(F9.2)' ) DNUM
  ELSE IF ( ANUM .GT. 0.09995 ) THEN
    WRITE ( C9, '(F9.3)' ) DNUM
  ELSE IF ( ANUM .GT. 0.009995 ) THEN
    WRITE ( C9, '(F9.4)' ) DNUM
  ELSE 
    WRITE ( C9, '(ES9.2)' ) DNUM
  END IF 
!
  C9REAL_D = ADJUSTL ( C9 ) 
!
END FUNCTION C9REAL_D

END MODULE C9REAL_GEN

