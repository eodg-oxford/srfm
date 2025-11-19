MODULE GRAPTH_SUB
CONTAINS
SUBROUTINE GRAPTH 
!
! VERSION
!   24APR24 AD Checked.
!   05MAR19 AD Separated out from LIMPTH for 2D atmospheres. Checked.
!   01JUL17 AD F90 conversion of tanpth.for. Checked.
!
! DESCRIPTION
!   Set 2D Limb-viewing paths
!   Called by RFMPTH, JACPTH if GRA flag enabled.
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
    USE SFCCOM_DAT, ONLY: RFLSFC ! Diffuse surface reflection
!
! SUBROUTINES
    USE GRASUM_SUB ! Construct ray paths through 2D atmosphere
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IATM   ! Atmospheric layer counter
    INTEGER(I4) :: IDIR   ! Direction pointer for path 
    INTEGER(I4) :: IDIR1  ! Initial value of IDIR 
    INTEGER(I4) :: IDIR2  ! Final value of IDIR 
    INTEGER(I4) :: IPTH   ! Path segment coutner
    INTEGER(I4) :: ITAN   ! Tangent path counter
    INTEGER(I4) :: IVMR   ! Gas counter
    REAL(R4)    :: LEN    ! Path length [km]
    REAL(R4)    :: PBOT   ! Horiz.angle [deg] at base of layer
    REAL(R4)    :: PTOP   ! Horiz.angle [deg] at top of layer
    REAL(R4)    :: TBOT   ! Zenith angle [deg] at base of layer
    REAL(R4)    :: TTOP   ! Zenith angle [deg] at top of layer
    REAL(R4)    :: ZBOT   ! Altitude [km] at base of layer
    REAL(R4)    :: ZTOP   ! Altitude [km] at top of layer
    REAL(R4), ALLOCATABLE :: AMT(:) ! Absorber amounts [kmol/cm2]
    REAL(R4), ALLOCATABLE :: ECG(:) ! CG Pressures [atm]
    REAL(R4), ALLOCATABLE :: PCG(:) ! CG Part.Pres [atm]
    REAL(R4), ALLOCATABLE :: PSU(:) ! Absorber-weighted horiz.angle [deg]
    REAL(R4), ALLOCATABLE :: TCG(:) ! CG Temperatures [K]
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  USEDIR = .TRUE.   ! Define paths as direction-dependent
!
! Calculate no of paths required
  NPTH = 0
  DO ITAN = 1, MTAN           
    IF ( TAN(ITAN)%CLC ) THEN
      IF ( TAN(ITAN)%SKY .OR. &
           TAN(ITAN)%ISK .GT. 0 .OR. &
           ( TAN(ITAN)%SFC .AND. .NOT. RFLSFC ) ) THEN ! End path at surface
        NPTH = NPTH + NVMR * ( NATM - TAN(ITAN)%IAT ) 
      ELSE                                    ! path either side of tan.pt/sfc
        NPTH = NPTH + 2 * NVMR * ( NATM - TAN(ITAN)%IAT ) 
      END IF
    END IF
  END DO 
!
  ALLOCATE ( PTH(NPTH) )
!
  IPTH = 0
!
  ALLOCATE ( AMT(NVMR), ECG(NVMR), PCG(NVMR), TCG(NVMR), PSU(NVMR) )
  DO ITAN = 1, MTAN
    IF ( .NOT. TAN(ITAN)%CLC ) CYCLE
    IF ( TAN(ITAN)%SKY ) THEN        ! Used for diffuse sky radiance
      IDIR1 = -1                     ! so just define downward path
      IDIR2 = -1
    ELSE IF ( TAN(ITAN)%ISK .GT. 0 ) THEN ! Downward path is diffuse sky
      IDIR1 = 1                      ! so just define upward path
      IDIR2 = 1
    ELSE IF ( TAN(ITAN)%SFC .AND. .NOT. RFLSFC ) THEN ! Path ends at surface
      IDIR1 = 1                      ! so just define upward path
      IDIR2 = 1
    ELSE         ! Default is to consider slant down and upward paths         
      IDIR1 = -1   
      IDIR2 = 1
    END IF 
    DO IDIR = IDIR1, IDIR2, 2
! RFM ray paths from space towards observer
! If IDIR=+1, ray from t.p. towards observer, so T is 0:-90 (a/c from zen)
! If IDIR=-1, ray from t.p. away from observer, so T is 0:90 (c/w from zen)
      TTOP = + SIGN ( SNGL ( ASIN ( TAN(ITAN)%SZN ) ) / DG2RAD, FLOAT(IDIR) )
      PTOP = TAN(ITAN)%PSI
      ZTOP = TAN(ITAN)%HGT 
      DO IATM = TAN(ITAN)%IAT, NATM-1
        ZBOT = ZTOP
        ZTOP = HGTATM(IATM+1) 
        PBOT = PTOP
        TBOT = TTOP
        CALL GRASUM ( NVMR, TAN(ITAN)%HGT, ZBOT, PBOT, TBOT, ZTOP, &
                      PTOP, TTOP, TCG, PCG, ECG, AMT, LEN, PSU ) 
        DO IVMR = 1, NVMR
          IPTH = IPTH + 1
          PTH(IPTH)%IGS = IVMR
          PTH(IPTH)%ITN = ITAN
          PTH(IPTH)%IAT = IATM
          PTH(IPTH)%ICL = 0
          PTH(IPTH)%IDR = IDIR
          PTH(IPTH)%NTE = NTEVMR(IVMR)
          PTH(IPTH)%TEM = TCG(IVMR)
          PTH(IPTH)%PRE = PCG(IVMR)
          PTH(IPTH)%PPA = ECG(IVMR)
          PTH(IPTH)%AMT = AMT(IVMR)
          PTH(IPTH)%RAY = LEN
          IF ( IDIR .EQ. -1 ) THEN
            PTH(IPTH)%PSI = PTOP ! Near=upper boundary for down path
          ELSE
            PTH(IPTH)%PSI = PBOT ! Near=lower boundary for up path
          END IF
          PTH(IPTH)%PSU = PSU(IVMR)
        END DO
      END DO
    END DO
  END DO
!
  IF ( IPTH .NE. NPTH ) STOP 'F-GRAPTH: Logical error'
!
END SUBROUTINE GRAPTH
END MODULE GRAPTH_SUB
