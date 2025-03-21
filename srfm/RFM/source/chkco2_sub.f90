MODULE CHKCO2_SUB
CONTAINS
SUBROUTINE CHKCO2
!
! VERSION
!   18JUL24 AD Checked.
!   01MAY17 AD F90 original. Tested.
!
! DESCRIPTION
!   Check if CO2 being used with CHI,MIX flags
!   Called once by DRVGAS
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE GASCOM_DAT ! Molecule and isotope data
    USE SHPCON_DAT ! Line-shape indices
    USE IDXCON_DAT, ONLY: IDXCO2 ! HITRAN/RFM index for CO2
!
! SUBROUTINES
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IGAS ! Index of CO2 in GASCOM
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IGAS = IGSMOL(IDXCO2)
  IF ( IGAS .GT. 0 ) THEN
    IF ( MIXFLG ) THEN
      CALL WRTLOG ( 'I-CHKCO2: Applying line-mixing for CO2 lines' )
      GAS(IGAS)%SHP = SHPMIX
    END IF
    IF ( CHIFLG ) THEN
      CALL WRTLOG ( 'I-CHKCO2: Applying Chi-Factor to CO2 lines' ) 
      GAS(IGAS)%SHP = SHPCHI  ! Chi-shape also allows for MIX flag
    END IF
  ELSE
    IF ( CHIFLG ) CALL WRTLOG ( 'W-CHKCO2: CHI flag selected but CO2 not used' )
    IF ( MIXFLG ) CALL WRTLOG ( 'W-CHKCO2: MIX flag selected but CO2 not used' )
  END IF 
!
END SUBROUTINE CHKCO2
END MODULE CHKCO2_SUB

