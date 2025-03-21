MODULE VALATM_SUB
CONTAINS
SUBROUTINE VALATM ( HGT, TEM, PRE, DNS, VMR )
!
! VERSION
!   11JAN24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Interpolate values from atmospheric profiles
!   Called by RAYSUM
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
!
! SUBROUTINES
    USE VAL1DI_GEN ! Interpolate value from 1D array
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN)  :: HGT    ! Altitude [km]
    REAL(R4), INTENT(OUT) :: TEM    ! Temperature [K]
    REAL(R4), INTENT(OUT) :: PRE    ! Pressure [mb]
    REAL(R4), INTENT(OUT) :: DNS    ! Density [molec/cm^3]
    REAL(R4), INTENT(OUT) :: VMR(:) ! VMR [ppmv]
!
! LOCAL VARIABLES
    INTEGER(I4) :: IVMR ! Counter for different species
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  TEM = VAL1DI ( HGTATM, HGT, TEMATM, .FALSE. ) 
  PRE = VAL1DI ( HGTATM, HGT, PREATM, .TRUE. ) 
  DNS = VAL1DI ( HGTATM, HGT, DNSATM, .TRUE. ) 
  DO IVMR = 1, NVMR
    VMR(IVMR) = VAL1DI ( HGTATM, HGT, VMRATM(:,IVMR), .NOT. LINVMR(IVMR) )
  END DO
!
END SUBROUTINE VALATM
END MODULE VALATM_SUB
