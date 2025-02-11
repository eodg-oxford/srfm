MODULE VRTPTH_SUB
CONTAINS
SUBROUTINE VRTPTH 
!
! VERSION
!   20SEP21 AD Use ICLPTH instead of ADDCLC. Checked.
!   05MAR19 AD Use single set of paths for all ray directions. Checked.
!   02MAY18 AD Bug#2 - correct SECANG if user-specified elevation angle
!   15DEC17 AD F90 conversion of nadpth.for. Checked.
!
! DESCRIPTION
!   Determine path segments along vertical path
!   Called by RFMPTH, JACPTH for plane-parallel geometry
!   A path segment is defined for each combination of absorber and atmos layer.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE FLGCOM_DAT ! Option flags
    USE PTHCOM_DAT ! Path segment data
    USE OBSCOM_DAT, ONLY: IATOBS ! Observer level (or 0 if unspecified)
!
! SUBROUTINES
    USE ICLPTH_FNC ! Index of corresponding (new) calculated path
    USE VRTSUM_SUB ! Curtis-Godson integrals for vertical path segment
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IATM   ! Atmospheric layer counter
    INTEGER(I4) :: IATM1  ! Start for atmospheric layer integration
    INTEGER(I4) :: IDIR   ! Ray direction -1,0, or 1
    INTEGER(I4) :: IPTH   ! Path segment coutner
    INTEGER(I4) :: IVMR   ! Gas counter
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Calculate no of paths required
  IF ( ZENFLG .AND. OBSFLG ) THEN
    NPTH = NVMR * ( NATM - IATOBS ) 
    IATM1 = IATOBS
  ELSE
    NPTH = NVMR * ( NATM - IATSFC )
    IATM1 = IATSFC
  END IF
!
  ALLOCATE ( PTH(NPTH) )
!
  IDIR = 0                   ! direction of radiation undefined (HOMFLG)
  IF ( NADFLG ) IDIR = 1     ! radiation moving upwards through atmosphere
  IF ( ZENFLG ) IDIR = -1    ! radiation moving downwards through atmosphere
  IPTH = 0
  DO IATM = IATM1, NATM-1
    DO IVMR = 1, NVMR
      IPTH = IPTH + 1
      PTH(IPTH)%IDR = IDIR
      PTH(IPTH)%IGS = IVMR
      PTH(IPTH)%ITN = 1
      PTH(IPTH)%IAT = IATM
      PTH(IPTH)%NTE = NTEVMR(IVMR)
      CALL VRTSUM ( IVMR, IATM, PTH(IPTH)%TEM, PTH(IPTH)%PRE, PTH(IPTH)%PPA, &
                    PTH(IPTH)%AMT, PTH(IPTH)%RAY )
      IF ( ZENFLG ) THEN
        PTH(IPTH)%PSI = 180.0 
      ELSE
        PTH(IPTH)%PSI = 0.0
      ENDIF
      PTH(IPTH)%ICL = ICLPTH ( PTH(IPTH) ) 
    END DO
  END DO
!
END SUBROUTINE VRTPTH
END MODULE VRTPTH_SUB
