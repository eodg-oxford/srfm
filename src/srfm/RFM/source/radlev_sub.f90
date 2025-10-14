MODULE RADLEV_SUB
CONTAINS
SUBROUTINE RADLEV ( ITAN, IATM, IDIR ) 
!
! VERSION
!   07JAN24 AD Checked.
!   01JUL17 AD F90 version. Checked.
!
! DESCRIPTION
!   Save radiances at intermediate output levels
!   Called by SPCRAD for each widemesh interval for each atm level.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FULCOM_DAT ! Full grid data
    USE LEVCOM_DAT ! Intermediate output levels
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN)  :: ITAN ! Index of full ray path
    INTEGER(I4), INTENT(IN)  :: IATM ! Index of atmos.profile level
    INTEGER(I4), INTENT(IN)  :: IDIR ! Direction: -1=downward, +1=upward
!
! LOCAL VARIABLES
    INTEGER(I4) :: ILEV ! Index of output level
    INTEGER(I4) :: JTAN ! Index of secondary ray path
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  JTAN = 0
  IF ( .NOT. ANY ( LEV%IAT .EQ. IATM ) ) RETURN
!
  DO ILEV = 1, NLEV
    IF ( LEV(ILEV)%IAT .EQ. IATM .AND. LEV(ILEV)%IDR .EQ. IDIR ) THEN
      JTAN = ITNLEV(ITAN,ILEV) 
      IF ( JTAN .EQ. 0 ) STOP 'F-RADLEV; Logical error'
      OPTFUL(IFUL1:IFUL2,JTAN) = OPTFUL(IFUL1:IFUL2,ITAN)
      TRAFUL(IFUL1:IFUL2,JTAN) = TRAFUL(IFUL1:IFUL2,ITAN)
      RADFUL(IFUL1:IFUL2,JTAN) = RADFUL(IFUL1:IFUL2,ITAN)
      RETURN
    END IF
  END DO
  STOP 'F-RADLEV: Logical error'
!
END SUBROUTINE RADLEV
END MODULE RADLEV_SUB

