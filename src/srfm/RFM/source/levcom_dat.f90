MODULE LEVCOM_DAT
!
! VERSION
!   19JUL26 AD Checked.
!   01AUG25 AD Add BOALEV, TOALEV, %STR, LENHGT
!   19FEB25 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Intermediate output levels
!   Loaded by ADDLEV.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL CONSTANTS
    USE LENHGT_DAT ! Max length of height component of RFM output filenames
!
  IMPLICIT NONE
  SAVE
  PUBLIC :: LEVCOM_RESET
!
  TYPE :: LEVTYP
    INTEGER(I4)       :: IAT ! Atmospheric profile level
    INTEGER(I4)       :: IDR ! Ray direction -1=toa downwards +1=upwards
    REAL(R4)          :: HGT ! Altitude [km]
    CHARACTER(LENHGT) :: STR ! Hgt as string in output filename
  END TYPE LEVTYP
!
! GLOBAL VARIABLES
    TYPE(LEVTYP), ALLOCATABLE :: LEV(:)
!
    LOGICAL     :: BOALEV = .FALSE.         ! T=BOA selected for output 
    LOGICAL     :: TOALEV = .FALSE.         ! T=TOA selected for output 
    INTEGER(I4) :: NLEV = 0                 ! No. output levels
    INTEGER(I4), ALLOCATABLE :: ITNLEV(:,:) ! [MTAN,NLEV] Indices of lev. rays
!
CONTAINS
!
  SUBROUTINE LEVCOM_RESET()
    IF ( ALLOCATED ( LEV ) ) DEALLOCATE ( LEV )
    IF ( ALLOCATED ( ITNLEV ) ) DEALLOCATE ( ITNLEV )
    BOALEV = .FALSE.
    TOALEV = .FALSE.
    NLEV   = 0
  END SUBROUTINE LEVCOM_RESET
!
END MODULE LEVCOM_DAT
