MODULE LIMPTH_SUB
CONTAINS
SUBROUTINE LIMPTH 
!
! VERSION
!   29APR24 AD Checked.
!   05MAR19 AD Simplified to exclude 2D paths. Checked.
!   01JUL17 AD F90 conversion of tanpth.for. 
!
! DESCRIPTION
!   Set Limb-viewing paths
!   Called by RFMPTH, JACPTH
!   A path segment is defined for each combination of absorber and atmos layer.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE PTHCOM_DAT ! Path segment data
    USE TANCOM_DAT ! Tangent path data
    USE PHYCON_DAT, ONLY: DG2RAD ! [rad]/[deg] conversion factor
!
! SUBROUTINES
    USE RAYSUM_SUB ! CG integrals in ray path segment
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IATM ! Atmospheric layer counter
    INTEGER(I4) :: IPTH ! Path segment coutner
    INTEGER(I4) :: ITAN ! Tangent path counter
    INTEGER(I4) :: IVMR ! Gas counter
    REAL(R4)    :: LEN  ! Path length [km]
    REAL(R4)    :: ZBOT ! Altitude [km] at base of layer
    REAL(R4)    :: ZTOP ! Altitude [km] at top of layer
    REAL(R8)    :: DSNBOT ! Sine of zenith angle at base of layer
    REAL(R8)    :: DSNTOP ! Sine of zenith angle at top of layer
    REAL(R4), ALLOCATABLE :: AMT(:) ! Absorber amounts [kmol/cm2]
    REAL(R4), ALLOCATABLE :: ECG(:) ! CG Pressures [atm]
    REAL(R4), ALLOCATABLE :: PCG(:) ! CG Part.Pres [atm]
    REAL(R4), ALLOCATABLE :: PSU(:) ! Absorber-weighted horiz.angle [deg]
    REAL(R4), ALLOCATABLE :: TCG(:) ! CG Temperatures [K]
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Calculate no of paths required
  NPTH = 0
  DO ITAN = 1, MTAN           
    IF ( TAN(ITAN)%CLC ) NPTH = NPTH + NVMR * ( NATM - TAN(ITAN)%IAT ) 
  END DO 
!   
  ALLOCATE ( PTH(NPTH) )
!
  IPTH = 0
!
  ALLOCATE ( AMT(NVMR), ECG(NVMR), PCG(NVMR), TCG(NVMR), PSU(NVMR) )
  DO ITAN = 1, MTAN
    IF ( .NOT. TAN(ITAN)%CLC ) CYCLE
    DSNTOP = TAN(ITAN)%SZN
    ZTOP = TAN(ITAN)%HGT 
    DO IATM = TAN(ITAN)%IAT, NATM-1
      ZBOT = ZTOP
      ZTOP = HGTATM(IATM+1) 
      DSNBOT = DSNTOP
      CALL RAYSUM ( NVMR, ZBOT, DSNBOT, ZTOP, DSNTOP, &
                    TCG, PCG, ECG, AMT, LEN )
      DO IVMR = 1, NVMR
        IPTH = IPTH + 1
        PTH(IPTH)%IGS = IVMR
        PTH(IPTH)%ITN = ITAN
        PTH(IPTH)%IAT = IATM
        PTH(IPTH)%ICL = 0     ! Calc paths assigned separately
        PTH(IPTH)%IDR = 0
        PTH(IPTH)%NTE = NTEVMR(IVMR)
        PTH(IPTH)%TEM = TCG(IVMR)
        PTH(IPTH)%PRE = PCG(IVMR)
        PTH(IPTH)%PPA = ECG(IVMR)
        PTH(IPTH)%AMT = AMT(IVMR)
        PTH(IPTH)%RAY = LEN
        PTH(IPTH)%PSI = SNGL ( ASIN ( DSNTOP ) ) / DG2RAD
      END DO
    END DO
  END DO
!
  IF ( IPTH .NE. NPTH ) STOP 'F-LIMPTH: Logical error'
!
END SUBROUTINE LIMPTH
END MODULE LIMPTH_SUB
