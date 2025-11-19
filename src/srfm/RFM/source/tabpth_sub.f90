MODULE TABPTH_SUB
CONTAINS
SUBROUTINE TABPTH 
!
! VERSION
!   25JUN24 AD Bug#46: Write path data to PTHCOM rather than CLCCOM
!   01JUL23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Set paths for TABulated abs.coeff. calculations
!   Called by RFMPTH if TAB flag enabled.
!   A path is defined for each absorber, pressure and temperature.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE PTHCOM_DAT ! Path segment data
    USE TABCOM_DAT ! Axes for tabulated absorption coefficients
    USE PHYCON_DAT, ONLY: ATMB ! Std Atmos. Pressure [mb]
!
! SUBROUTINES
   USE ICLPTH_FNC ! Index of corresponding (new) calculated path
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IPTH  ! Counter for calculated paths
    INTEGER(I4) :: IPTAB ! Counter for P-axis points 
    INTEGER(I4) :: IQTAB ! Counter for Q-axis points 
    INTEGER(I4) :: ITTAB ! Counter for T-axis points 
    INTEGER(I4) :: IVMR  ! Gas counter
    REAL(R4)    :: PRE   ! Pressure axis value converted to [atm]
    REAL(R4)    :: TEM   ! Temperature axis value [K]
    REAL(R4)    :: VMR   ! Volume mixing ratio [ppv]
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! Calculate the total number of different paths required
  NPTH = NATAB * NTAB
  ALLOCATE ( PTH(NPTH) ) 
!
  IPTH = 0
  DO IVMR = 1, NTAB                     ! Loop over molecules
    DO IPTAB = 1, NPTAB
      PRE = PAXTAB(IPTAB) / ATMB        ! Convert [mb] to [atm]
      TEM = TEMTAB(IPTAB)               
      VMR = VMRTAB(IPTAB,IVMR) * 1.0E-6 ! Convert [ppmv] to [ppv]
      DO ITTAB = 1, NTTAB
        TEM = TAXTAB(ITTAB)             
        IF ( OFFTAB ) TEM = TEM + TEMTAB(IPTAB)
        DO IQTAB = 1, NQTAB
          IPTH = IPTH + 1
          PTH(IPTH)%IGS = IVMR
          PTH(IPTH)%PRE = PRE
          PTH(IPTH)%TEM = TEM
          PTH(IPTH)%PPA = VMR * PRE * QAXTAB(IQTAB) * 0.01  ! 0.01 for %
          PTH(IPTH)%AMT = 1.0E-4 ! 10^-4 kmol/cm^2 = 1 kmol/m^2 in output
          PTH(IPTH)%ICL = ICLPTH ( PTH(IPTH), .TRUE. ) 
! The following are not used, but ensure sensibly defined for PTH flag output
          PTH(IPTH)%ITN = 1   
          PTH(IPTH)%IDR = 1   
          PTH(IPTH)%PSI = 0.0
          PTH(IPTH)%RAY = 1.0
        END DO
      END DO
    END DO
  END DO
!
END SUBROUTINE TABPTH
END MODULE TABPTH_SUB
