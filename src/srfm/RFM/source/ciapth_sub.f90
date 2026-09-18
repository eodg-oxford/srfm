MODULE CIAPTH_SUB
CONTAINS
SUBROUTINE CIAPTH 
!
! VERSION
!   15AUG26 AD Checked
!   18AUG25 AD Bug#51: Rewritten to properly allow for multiple CIA molec.
!   08OCT24 AD Checked.
!   20SEP21 AD Use ICLPTH instead of ADDCLC. Checked.
!   01MAY17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   Set CIA paths
!   Called once by RFMPTH if CIA flag enabled.
!   Sets up a list of additional path data in CIPCOM for all molecules modelled
!   by collision-induced absorption data (.cia files).
!   Unlike other absorption data, CIA data is a depends on the product of the
!   absorber amounts of two molecules.
!   It is possible that several paths using the same absorber amount for a 
!   molecule share the same calculated path. However, since there is no test
!   on whether molec#2 is similar for such cases, this module also ensures
!   that each path has its own associated calc.pth.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CIACOM_DAT ! Collision-induced absorption data
    USE CIPCOM_DAT ! CIA path data
    USE PTHCOM_DAT ! Path segment data
    USE CLCCOM_DAT, ONLY: NCLC       ! No. calculated paths
    USE PHYCON_DAT, ONLY: ATMB, RGAS ! [mb/atm], Gas Constant
!
! SUBROUTINES
    USE ICLPTH_FNC ! Index of corresponding (new) calculated path
    USE IDXPTH_FNC ! Index in PTHCOM of tan/atm/gas/dir
!                  
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: ICIA ! Counter for CIA tabulations
    INTEGER(I4) :: ICIP ! Counter for CIA paths
    INTEGER(I4) :: ICLC ! Index of equivalent calc path
    INTEGER(I4) :: IG2  ! GASCOM index of molec#2
    INTEGER(I4) :: IGAS ! Absorber index in GASCOM of molec#1
    INTEGER(I4) :: IPTH ! Index of path segement for IGAS
    INTEGER(I4) :: JCIA ! Secondary counter for CIA tabulations
    INTEGER(I4) :: JGAS ! Absorber index in GASCOM of molec#2
    INTEGER(I4) :: JPTH ! Index of equiv. path segment for molec#2
    INTEGER(I4) :: ID2LST(NCIA) ! List of HITRAN/RFM IDs for molec#2
    INTEGER(I4) :: IG2LST(NCIA) ! List of GASCOM indices for molec#2
    INTEGER(I4) :: NIG2 ! No. of different molec#2 for this molec#1 
    LOGICAL,      ALLOCATABLE :: USECLC(:) ! Flags marking used Calc paths    
    TYPE(CIPTYP), ALLOCATABLE :: CIPSAV(:) ! Saved CIP during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Ensure that just one calc path associated with each path segment containing
! CIA molecule since scaling by molec#2 uses absorber amount
  ALLOCATE ( USECLC(NCLC) ) ; USECLC = .FALSE.

  DO ICIA = 1, NCIA
    IGAS = CIA(ICIA)%IG1
    IF ( ANY ( CIA(1:ICIA-1)%IG1 .EQ. IGAS ) ) CYCLE ! already checked
!
! Create list of different CIA 2nd abs. associated with this 1st abs
! NB ID2LST,IG2LST declared size NCIA to ensure always large enough.
    NIG2 = 0
    DO JCIA = 1, NCIA
      IF ( CIA(JCIA)%IG1 .NE. IGAS ) CYCLE
      IG2 = CIA(JCIA)%IG2
      IF ( .NOT. ANY ( IG2LST(1:NIG2) .EQ. IG2 ) ) THEN
        NIG2 = NIG2 + 1
        IG2LST(NIG2) = IG2
        ID2LST(NIG2) = CIA(JCIA)%ID2
      END IF
    END DO
!
! Check that any PTH with CIA 1st abs gas has a unique CLC path
    DO IPTH = 1, NPTH
      IF ( PTH(IPTH)%IGS .EQ. IGAS ) THEN
        ICLC = PTH(IPTH)%ICL
        IF ( USECLC(ICLC) ) THEN   ! Used already, so assign a new CLC path
          ICLC = ICLPTH ( PTH(IPTH), .TRUE. )
          PTH(IPTH)%ICL = ICLC
        ELSE
          USECLC(ICLC) = .TRUE.    ! Flag as used
        END IF

        IF ( ALLOCATED ( CIP ) ) CALL MOVE_ALLOC ( CIP, CIPSAV )
        ALLOCATE ( CIP(NCIP+NIG2) ) 
        IF ( ALLOCATED ( CIPSAV ) ) CIP(1:NCIP) = CIPSAV
        ICIP = NCIP
        NCIP = NCIP + NIG2
! Loop through all 2nd CIA absorbers associated with this 1st abs.
        DO IG2 = 1, NIG2
          JGAS = IG2LST(IG2)
          IF ( JGAS .EQ. IGAS ) THEN    
            JPTH = IPTH
          ELSE
            JPTH = IDXPTH ( PTH(IPTH)%ITN, PTH(IPTH)%IAT, JGAS, PTH(IPTH)%IDR )
          END IF
          ICIP = ICIP + 1
          CIP(ICIP)%ICL = ICLC
          CIP(ICIP)%ID1 = CIA(ICIA)%ID1
          CIP(ICIP)%ID2 = ID2LST(IG2)
          CIP(ICIP)%TEM = PTH(IPTH)%TEM
! Factor 100 to convert from mb to Pa, 1E-6 to convert from /m3 to /cm3
! Use the JPTH amount to include the path length integration, hence the need
! for individual calc paths for each segment.
! Since CLC path is later multiplied by IPTH AMT, divide by it here
          CIP(ICIP)%AM2 = PTH(IPTH)%PPA * ATMB * 100.0 * 1.0E-6 / RGAS &
                          / PTH(IPTH)%TEM * PTH(JPTH)%AMT / PTH(IPTH)%AMT
        END DO
      END IF
    END DO
  END DO
!      
END SUBROUTINE CIAPTH
END MODULE CIAPTH_SUB
