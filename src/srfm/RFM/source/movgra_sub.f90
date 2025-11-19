MODULE MOVGRA_SUB
CONTAINS
SUBROUTINE MOVGRA ( IPSI, JPSI )
!
! VERSION
!   19SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Move profiles within GRACOM arrays
!   Called by ATMGRA, RFMPRF (not sure if this is actually necessary).
!   An index IPSI or JPSI =0 means use profiles in ATMCOM
!   Other indices refer to position in GRACOM
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE GRACOM_DAT ! Atmospheric 2-D field
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IPSI ! Old index of profile data
    INTEGER(I4), INTENT(IN) :: JPSI ! New index of profile data
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( JPSI .EQ. 0 ) THEN                      ! Move IPSI to ATM
    TEMATM = TEMGRA(:,IPSI)
    PREATM = PREGRA(:,IPSI)
    IF ( IAXVMR .GT. 0 ) EXTATM = EXTGRA(:,IPSI)
    VMRATM = VMRGRA(:,IPSI,:)
    IF ( NVIB .GT. 0 ) VIBATM = VIBGRA(:,IPSI,:)
  ELSE IF ( IPSI .EQ. 0 ) THEN                 ! Move ATM to JPSI
    TEMGRA(:,JPSI) = TEMATM
    PREGRA(:,JPSI) = PREATM
    IF ( IAXVMR .GT. 0 ) EXTGRA(:,JPSI) = EXTATM
    VMRGRA(:,JPSI,:) = VMRATM
    IF ( NVIB .GT. 0 ) VIBGRA(:,IPSI,:) = VIBATM
  END IF
!
END SUBROUTINE MOVGRA
END MODULE MOVGRA_SUB

