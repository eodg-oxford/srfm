MODULE CHKSFC_SUB
CONTAINS
SUBROUTINE CHKSFC
!
! VERSION
!   20APR24 AD Checked.
!   05MAR19 AD F90 original. Checked.
!
! DESCRIPTION
!   Check vertical path defined for surface reflection
!   Called once by DRVCHK if SFC flag is enabled. 
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE TANCOM_DAT ! Tangent path data
    USE ATMCOM_DAT, ONLY: IATSFC, HGTSFC ! Surface profile level#, height.
    USE FLGCOM_DAT, ONLY: NADFLG         ! T=nadir viewing in plane par.atm.
    USE SFCCOM_DAT, ONLY: DIFSFC, RFLSFC ! Surface parameters
!
! SUBROUTINES
    USE ADDTAN_SUB ! Add tangent ray path
    USE INIQAD_SUB ! Initialise Gaussian quadrature for flux calculations
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE 
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( .NOT. RFLSFC ) RETURN      ! Non-reflecting surface
!
  IF ( .NOT. ANY ( TAN%SFC ) ) RETURN ! No paths intersect surface
!
  IF ( DIFSFC ) THEN
    CALL WRTLOG ( 'I-CHKSFC: Assuming surface reflectance: Diffuse' )
    CALL INIQAD
    CALL ADDTAN ( 1, .TRUE. ) 
    TAN(MTAN)%SFC = .TRUE.
    TAN(MTAN)%SKY = .TRUE.
    TAN(MTAN)%IAT = IATSFC
    IF ( NADFLG ) THEN
      TAN(MTAN)%ITN = 1
    ELSE
      TAN(MTAN)%ITN = 0
    END IF
    TAN(MTAN)%JDX = 0
    TAN(MTAN)%SZN = 0.0D0
    TAN(MTAN)%HGT = HGTSFC
    TAN(MTAN)%PSI = 0.0
    TAN(MTAN)%SZN = 0.0
! Set all surface intersecting paths to use this for diffuse sky radiance
    WHERE ( TAN%SFC ) TAN%ISK = MTAN
  ELSE
    CALL WRTLOG ( 'I-CHKSFC: Assuming surface reflectance: Specular' ) 
  END IF
!
END SUBROUTINE CHKSFC
END MODULE CHKSFC_SUB
