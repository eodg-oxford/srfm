MODULE DIMTAB_SUB
CONTAINS
SUBROUTINE DIMTAB 
!
! VERSION
!   01AUG24 AD Checked.
!   24JUN19 AD Interpolate LOG(PREATM) since LNPATM may not yet be set. Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Complete TAB header information
!   Called by DRVDIM if TAB flag enabled.
!   Assumes axes have been defined, sets p, T, VMR profiles
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE FLGCOM_DAT ! Option flags
    USE TABCOM_DAT ! Axes for tabulated absorption coefficients
    USE PHYCON_DAT, ONLY: GHZ2CM ! GHz to cm-1 conv.fac (~1/30)
    USE RFMLUN_DAT, ONLY: LUNNXT ! Next free LUN
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
    USE INTERP_GEN ! Interpolate array
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IVMR          ! Gas counter
    REAL(R4)      :: LNPTAB(NPTAB) ! Log ( P-axis [mb] values )
    REAL(R4)      :: PMAX          ! Max value [mb] of P-axis 
    REAL(R4)      :: PMIN          ! Min value [mb] of P-axis
    CHARACTER(80) :: WRNMSG        ! Warning message sent to log file
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! Calculate the total number of tabulated points
  NATAB = NPTAB * NTTAB * NQTAB
!
! Compare tabulated pressures with atmospheric profile range and send warning
! message if P-axis extends outside profile (unless HOMFLG enabled, which would
! always generate a warning - assume the user knows what they're doing!)
  IF ( NATM .GT. 1 ) THEN
    PMIN = MIN ( PAXTAB(1), PAXTAB(NPTAB) )
    IF ( PMIN .LT. PREATM(NATM) ) THEN
      WRNMSG = 'W-DIMTAB: Min P-axis value=' // C9REAL(PMIN) // & 
               '[mb] < Min P profile value=' // C9REAL(PREATM(NATM)) // '[mb]'
      CALL WRTLOG ( WRNMSG )
    END IF
    PMAX = MAX ( PAXTAB(1), PAXTAB(NPTAB) )
    IF ( PMAX .GT. PREATM(1) ) THEN
      WRNMSG = 'W-DIMTAB: Max P-axis value=' // C9REAL(PMAX) // & 
               '[mb] > Max P profile value=' // C9REAL(PREATM(1)) // '[mb]'
      CALL WRTLOG ( WRNMSG )
    END IF
  END IF
!
! Set up log(p) array for p-axis values
  LNPTAB = LOG ( PAXTAB )
!
! Interpolate atmospheric T profile values to lnp-axis
  ALLOCATE ( TEMTAB(NPTAB) )
  TEMTAB = INTERP ( LOG(PREATM), LNPTAB, TEMATM ) 
!
! Determine if temperature axis is absolute or relative to TEMTAB
  OFFTAB = MIN ( TAXTAB(1), TAXTAB(NTTAB) ) .LE. 0.0
!
! Interpolate VMR [ppmv] profile to lnp axis, linearly or log in VMR
  ALLOCATE ( VMRTAB(NPTAB,NVMR) )
  DO IVMR = 1, NVMR
    VMRTAB(:,IVMR) = INTERP ( LOG(PREATM), LNPTAB, VMRATM(:,IVMR), &
                              .NOT. LINVMR(IVMR) ) 
  END DO
!
! Reserve separate LUN for each molecule
  LUNTAB = LUNNXT
  NTAB = NVMR
  LUNNXT = LUNNXT + NTAB
!
! Scaling factor if output spectra in GHz
  IF ( GHZFLG ) THEN
    GHZTAB = 1.0D0 / GHZ2CM
  ELSE
    GHZTAB = 1.0D0
  END IF
!
  BINTAB = BINFLG
!
END SUBROUTINE DIMTAB
END MODULE DIMTAB_SUB
