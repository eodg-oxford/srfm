MODULE FLXPTH_SUB
CONTAINS
SUBROUTINE FLXPTH
!
! VERSION
!   20SEP21 AD Replace ADDCLC with ICLPTH. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Set up paths for flux calculations
!   Called by JACPTH.
!   All paths are assigned to a single, vertical TAN ray ITN=1.
!   Combination of MTX+TRA flags is a special case where extra paths are 
!   required to store parameters at output levels for escape functions.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE PTHCOM_DAT ! Path segment data
    USE TANCOM_DAT ! Tangent path data
    USE FLGCOM_DAT, ONLY: MTXFLG, TRAFLG ! Option flags
    USE PHYCON_DAT, ONLY: ATMB ! [mb]/[atm] conversion factor
!
! SUBROUTINES
    USE ICLPTH_FNC ! Index of corresponding (new) calculated path
    USE VRTSUM_SUB ! Curtis-Godson integrals for vertical path segment
!
  IMPLICIT NONE 
!
! LOCAL VARIABLES
    INTEGER(I4) :: IATM ! Counter for atmospheric profile levels
    INTEGER(I4) :: IPTH ! Counter for path segments
    INTEGER(I4) :: ITAN ! Counter for Flux matrix output levels
    INTEGER(I4) :: IVMR ! Counter for absorbing species
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NPTH = NVMR * ( NATM - 1 ) 
  IF ( MTXFLG .AND. TRAFLG ) NPTH = NPTH + ( NTAN * NVMR )
  ALLOCATE ( PTH(NPTH) )
  IPTH = 0
  DO IATM = 1, NATM-1
    DO IVMR = 1, NVMR
      IPTH = IPTH + 1
      PTH(IPTH)%IGS = IVMR
      PTH(IPTH)%ITN = 1
      PTH(IPTH)%IAT = IATM
      PTH(IPTH)%NTE = NTEVMR(IVMR)
      CALL VRTSUM ( IVMR, IATM, PTH(IPTH)%TEM, PTH(IPTH)%PRE, PTH(IPTH)%PPA, &
                    PTH(IPTH)%AMT, PTH(IPTH)%RAY )
      PTH(IPTH)%PSI = 1.0
      PTH(IPTH)%ICL = ICLPTH ( PTH(IPTH) ) ! Should always add a new Calc path 
    END DO
  END DO
!
! If MTX+TRA flags, also create paths storing parameters at each output *level*
! for escape function FLXEFN which will be output as tra. files 
  IF ( MTXFLG .AND. TRAFLG ) THEN
    DO ITAN = 1, NTAN                ! NTAN= no. output levels
      IATM = TAN(ITAN)%IAT 
      DO IVMR = 1, NVMR
        IPTH = IPTH + 1
        PTH(IPTH)%IDR = 2            ! 2=flag for particular set of paths
        PTH(IPTH)%IGS = IVMR
        PTH(IPTH)%ITN = 1
        PTH(IPTH)%IAT = IATM
        PTH(IPTH)%TEM = TEMATM(IATM)
        PTH(IPTH)%PRE = PREATM(IATM) / ATMB
        PTH(IPTH)%PPA = VMRATM(IATM,IVMR) * PREATM(IATM) / ATMB
        PTH(IPTH)%AMT = VMRATM(IATM,IVMR)  * DNSATM(IATM)  ! set AMT=molec/cm3
        PTH(IPTH)%RAY = 1.0E-5                             ! Set path = 1cm
        PTH(IPTH)%ICL = ICLPTH ( PTH(IPTH) ) 
      END DO
    END DO
  END IF     
!
END SUBROUTINE FLXPTH
END MODULE FLXPTH_SUB
