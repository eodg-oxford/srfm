MODULE HOMPTH_SUB
CONTAINS
SUBROUTINE HOMPTH 
!
! VERSION
!   20SEP21 AD Replace ADDCLC with ICLPTH. Checked.
!   05MAR19 AD Reduce to single TAN path. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Set Homogeneous paths
!   Called by RFMPTH if HOM flag enabled.
!   Called by JACPTH if HOM+JAC flags enabled.
!   A path is defined for each combination of absorber and atmos layer.
!   This sets absorber amounts for 1km path, assumes TAN%SEC scales length.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE PTHCOM_DAT ! Path segment data
    USE PHYCON_DAT, ONLY: ATMB, AVOG ! mb/atm, Avogadro's number
!
! SUBROUTINES
    USE ICLPTH_FNC ! Index of corresponding (new) calculated path
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER :: IPTH ! Path segment counter
    INTEGER :: IVMR ! Absorber counter
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NPTH = NVMR 
  ALLOCATE ( PTH(NPTH) )
!
  DO IVMR = 1, NVMR
    IPTH = IVMR
    PTH(IPTH)%IGS = IVMR
    PTH(IPTH)%ITN = 1
    PTH(IPTH)%NTE = NTEVMR(IVMR)
    PTH(IPTH)%IAT = 1
    PTH(IPTH)%TEM = TEMATM(1)
    PTH(IPTH)%PRE = PREATM(1) / ATMB       ! mb to atm
    PTH(IPTH)%PPA = VMRATM(1,IVMR) * 1.0E-6 * PREATM(1) / ATMB ! ppmv to ppv
    PTH(IPTH)%RAY = 1.0                    ! assume 1km path
! 0.1 is from cm/km (10^5) * ppv/ppmv (10^-6)
    PTH(IPTH)%AMT = 0.1 / AVOG * VMRATM(1,IVMR) * DNSATM(1) ! assume 1km path
    PTH(IPTH)%PSI = 0.0
    PTH(IPTH)%ICL = ICLPTH ( PTH(IPTH) ) 
  END DO
!
END SUBROUTINE HOMPTH
END MODULE HOMPTH_SUB
