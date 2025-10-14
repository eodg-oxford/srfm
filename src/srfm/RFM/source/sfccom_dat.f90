MODULE SFCCOM_DAT
!
! VERSION
!   03MAY24 AD Checked.
!   05MAR19 AD Add DIFSFC. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Surface parameters
!   Read from Driver Table if SFC flag is enabled.
!   Loaded/checked by SFCFIL
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL VARIABLES
    LOGICAL     :: DIFSFC = .FALSE.    ! T=diffuse sfc.scattering, F=specular
    LOGICAL     :: RFLSFC = .FALSE.    ! T=reflective surface, F=no reflection
    INTEGER(I4) :: NSFC = 0            ! No. spec pts defining sfc emissivity
    REAL(R4)    :: TEMSFC              ! Surface Temperature [K]
    REAL(R8), ALLOCATABLE :: EMSSFC(:) ! Surface Emissivity (0:1)
    REAL(R8), ALLOCATABLE :: WNOSFC(:) ! Spectral axis [cm-1]
!
END MODULE SFCCOM_DAT
