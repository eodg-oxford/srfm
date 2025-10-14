MODULE PRFGRA_SUB
CONTAINS
SUBROUTINE PRFGRA ( TYPE, IPRF )
!
! VERSION
!   18SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Copy atm profile from 1D to 2D field
!   General purpose module.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE GRACOM_DAT ! Atmospheric 2-D field
!
! SUBROUTINES
    USE ADDGRA_SUB ! Set horiz location in 2D field
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(3),          INTENT(IN) :: TYPE ! Type of profile
    INTEGER(I4), OPTIONAL, INTENT(IN) :: IPRF ! Profile index 
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL ADDGRA ! Sets MPSI for current PSI value
!
  SELECT CASE ( TYPE ) 
  CASE ( 'PRE' )      ! Pressure profile
    PREGRA(:,MPSI) = PREATM
    FLPGRA(MPSI) = .TRUE.
  CASE ( 'TEM' )      ! Temperature profile
    TEMGRA(:,MPSI) = TEMATM
    FLTGRA(MPSI) = .TRUE.
  CASE ( 'VMR' )      ! VMR profile
    VMRGRA(:,MPSI,IPRF) = VMRATM(:,IPRF)
    FLVGRA(MPSI,IPRF) = .TRUE.
  CASE ( 'VIB' )      ! Vib Tem profile
    VIBGRA(:,MPSI,IPRF) = VIBATM(:,IPRF)
    FLNGRA(MPSI,IPRF) = .TRUE.
  CASE ( 'AER' )      ! Aerosol extinction profile
    EXTGRA(:,MPSI) = EXTATM     ! converted to VMR in ATMAUX
    FLVGRA(MPSI,IPRF) = .TRUE.
  CASE DEFAULT
    STOP 'F-PRFGRA: Logical error'     
  END SELECT
!
END SUBROUTINE PRFGRA
END MODULE PRFGRA_SUB

