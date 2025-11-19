MODULE SPCCOM_DAT
!
! VERSION
!   15MAR23 AD Checked.
!   01MAY17 AD F90 Conversion. Checked.
!
! DESCRIPTION
!   Spectral range data.
!   NSPC is incremented by SPCLAB and can be cleared via SPCCOM_RESET when the
!   host application needs to reuse the Fortran library within the same process.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
  PUBLIC :: SPCCOM_RESET
!
! GLOBAL CONSTANTS
    INTEGER(I4), PARAMETER :: LENSPC = 8 ! Max length of spectral label
!
  TYPE :: SPCTYP
    INTEGER(I4)       :: IGD ! Index of grid file in GFLCOM for spectral range
    INTEGER(I4)       :: ILS ! Index of ILS function in ILSCOM for spec. range
    INTEGER(I4)       :: NPT ! No.output grid points
    REAL(R8)          :: WNL ! Specified lower wavenumber limit [/cm]
    REAL(R8)          :: WNR ! Calcuated Wavenumber resolution [/cm]
    REAL(R8)          :: WNU ! Specified upper wavenumber limit [/cm]
    REAL(R8)          :: WXL ! Extended lower wno including convol.
    REAL(R8)          :: WXU ! Extended upper wno including convol.
    CHARACTER(LENSPC) :: LAB ! Spectral range labels
  END TYPE SPCTYP
!
! GLOBAL VARIABLES
    TYPE(SPCTYP), ALLOCATABLE :: SPC(:) 
!
    INTEGER(I4) :: NSPC = 0 ! No.tabulated Spc.ranges
    REAL(R8)    :: WMNSPC   ! Lower Wavenumber reqd for any range
    REAL(R8)    :: WMXSPC   ! Upper Wavenumber reqd for any range
!
CONTAINS

  SUBROUTINE SPCCOM_RESET()
    IF ( ALLOCATED ( SPC ) ) DEALLOCATE ( SPC )
    NSPC   = 0
    WMNSPC = 0.0_R8
    WMXSPC = 0.0_R8
  END SUBROUTINE SPCCOM_RESET

END MODULE SPCCOM_DAT
