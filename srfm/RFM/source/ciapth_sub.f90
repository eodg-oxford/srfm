MODULE CIAPTH_SUB
CONTAINS
SUBROUTINE CIAPTH 
!
! VERSION
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
    INTEGER(I4) :: ICLC ! Index of equivalent calc path
    INTEGER(I4) :: ID1  ! HITRAN/RFM index of molec#1
    INTEGER(I4) :: ID2  ! HITRAN/RFM index of molec#2
    INTEGER(I4) :: IGAS ! Absorber index in GASCOM of molec#1
    INTEGER(I4) :: IGG  ! Counter for molec-molec pairs in CIA data
    INTEGER(I4) :: IPTH ! Index of path segement for IGAS
    INTEGER(I4) :: JGAS ! Absorber index in GASCOM of molec#2
    INTEGER(I4) :: JPTH ! Index of equiv. path segment for molec#2
    LOGICAL,      ALLOCATABLE :: USECLC(:) ! Flags marking used Calc paths    
    TYPE(CIPTYP), ALLOCATABLE :: CIPSAV(:) ! Saved CIP during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Ensure that just one calc path associated with each path segment containing
! CIA molecule since scaling by molec#2 uses absorber amount
  ALLOCATE ( USECLC(NCLC) ) ; USECLC = .FALSE.
!
  DO IGG = 1, NGGCIA
    IGAS = IGGCIA(1,IGG) 
    JGAS = IGGCIA(2,IGG) 
    ID1  = IDDCIA(1,IGG)
    ID2  = IDDCIA(2,IGG)
    DO IPTH = 1, NPTH
      IF ( PTH(IPTH)%IGS .NE. IGAS ) CYCLE
      ICLC = PTH(IPTH)%ICL
      IF ( USECLC(ICLC) ) THEN
        ICLC = ICLPTH ( PTH(IPTH), .TRUE. )
      ELSE
        USECLC(ICLC) = .TRUE.
      END IF
      IF ( JGAS .EQ. IGAS ) THEN    
        JPTH = IPTH
      ELSE
        JPTH = IDXPTH ( PTH(IPTH)%ITN, PTH(IPTH)%IAT, JGAS, PTH(IPTH)%IDR )
      END IF
      IF ( ALLOCATED ( CIP ) ) CALL MOVE_ALLOC ( CIP, CIPSAV )
      NCIP = NCIP + 1
      ALLOCATE ( CIP(NCIP) ) 
      IF ( ALLOCATED ( CIPSAV ) ) CIP(1:NCIP-1) = CIPSAV
      CIP(NCIP)%ICL = ICLC
      CIP(NCIP)%ID1 = ID1
      CIP(NCIP)%ID2 = ID2
      CIP(NCIP)%TEM = PTH(IPTH)%TEM
! Factor 100 to convert from mb to Pa, 1E-6 to convert from /m3 to /cm3
! Use the JPTH amount to include the path length integration, hence the need
! for individual calc paths for each segment.
! Since CLC path is later multiplied by IPTH AMT, divide by it here
      CIP(NCIP)%AM2 = PTH(IPTH)%PPA * ATMB * 100.0 * 1.0E-6 / RGAS &
                      / PTH(IPTH)%TEM * PTH(JPTH)%AMT / PTH(IPTH)%AMT
    END DO
  END DO
!      
END SUBROUTINE CIAPTH
END MODULE CIAPTH_SUB
