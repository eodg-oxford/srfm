MODULE VALGRA_SUB
CONTAINS
SUBROUTINE VALGRA ( HGT, PSI, TEM, PRE, DNS, VMR )
!
! VERSION
!   10JAN24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Interpolate values from 2D atmospheric field
!   Called by GRASUM
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GRACOM_DAT ! Atmospheric 2-D field
    USE ATMCOM_DAT ! Atmospheric profile data
!
! SUBROUTINES
    USE VAL2DI_FNC ! Interpolate value from 2D array
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN)  :: HGT    ! Altitude [km]
    REAL(R4), INTENT(IN)  :: PSI    ! LOS angle [deg]
    REAL(R4), INTENT(OUT) :: TEM    ! Temperature [K]
    REAL(R4), INTENT(OUT) :: PRE    ! Pressure [mb]
    REAL(R4), INTENT(OUT) :: DNS    ! Density [molec/cm^3]
    REAL(R4), INTENT(OUT) :: VMR(:) ! VMR [ppmv]
!
! LOCAL VARIABLES
    INTEGER(I4) :: IVMR   ! Counter for different VMR profiles
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  TEM = VAL2DI ( HGTATM, HGT, PSIGRA, PSI, TEMGRA, .FALSE. ) 
  PRE = VAL2DI ( HGTATM, HGT, PSIGRA, PSI, PREGRA, .TRUE. ) 
  DNS = VAL2DI ( HGTATM, HGT, PSIGRA, PSI, DNSGRA, .TRUE. ) 
  DO IVMR = 1, NVMR
    VMR(IVMR) = VAL2DI ( HGTATM, HGT, PSIGRA, PSI, VMRGRA(:,:,IVMR), &
                         .NOT. LINVMR(IVMR) ) 
  END DO
!
END SUBROUTINE VALGRA
END MODULE VALGRA_SUB

