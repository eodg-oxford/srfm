MODULE ATMLEV_SUB
CONTAINS
SUBROUTINE ATMLEV ( LEV, USEHGT, IDXATM )
!
! VERSION
!   25SEP24 AD Checked.
!   02AUG20 AD Correction: use GRDLEV rather than LEV for IDXATM,DELGRD
!   05AUG19 AD Add USEHGT argument, allow for pressure as well as altitude.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Find/insert atmospheric level for given altitude/pressure
!   General purpose module.
!   If the altitude is above the top of the atmosphere or below the base of
!     the atmosphere a fatal error message results
!   If the altitude exactly matches an existing atmospheric level, the index
!     IATM of that level is returned. 
!   If the altitude matches an existing atmospheric level within a given 
!     tolerance (set by TOLLEV) the altitude is adjusted to that level and 
!     the index returned.
!   If the altitude does not match, a new level is inserted and its index is
!     returned
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
!
! SUBROUTINES
    USE ADDATM_SUB ! Add extra level to atm profiles in ATMCOM
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE C9REAL_GEN ! Write real number as C*9 string
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4),    INTENT(INOUT) :: LEV    ! Altitude [km] or pressure [hPa]
    LOGICAL,     INTENT(IN)    :: USEHGT ! T=LEV is altitude, F=pressure
    INTEGER(I4), INTENT(OUT)   :: IDXATM ! Index of new atmospheric level
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: TOLLEV = 0.001 ! Fraction of layer thickness which 
!                                 can be ignored if HGT close to layer boundary
! LOCAL VARIABLES
    REAL(R4)      :: ALPHA  ! Profile interpolation fraction (0:1)
    REAL(R4)      :: DELGRD ! Difference between HGT and nearest profile level
    REAL(R4)      :: GRDLEV ! LEV in vertical grid coordinates (km or lnp)
    CHARACTER(80) :: WRNMSG ! Message sent to log file
    REAL(R4), POINTER :: GRDATM(:)     ! Vertical grid of profiles
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! First some checks of things that should already have been checked elsewhere
  IF ( USEHGT ) THEN
    IF ( .NOT. SETHGT ) STOP 'F-ATMLEV: Logical error#1' 
    IF ( LEV .GT. HGTTOA .OR. LEV .LT. HGTSFC ) STOP 'F-ATMLEV: Logical error#2'
  ELSE
    IF ( LEV .LT. PREATM(NATM) .OR. LEV .GT. PRESFC ) &
      STOP 'F-ATMLEV: Logical error#3'
  END IF
!
  IF ( USEHGT ) THEN
    GRDATM => HGTATM
    GRDLEV = LEV
  ELSE 
    GRDATM => LNPATM     ! LnP should be already set, doesn't require ATMAUX
    GRDLEV = LOG(LEV)
  END IF
!
  IDXATM = MINLOC ( ABS ( GRDATM - GRDLEV ), 1 ) 
  DELGRD = MINVAL ( ABS ( GRDATM - GRDLEV ) ) 
  IF ( DELGRD .EQ. 0.0 ) RETURN              ! Normal exit with exact match
!
  IF ( IDXATM .EQ. 1 ) THEN 
    ALPHA = DELGRD / ( GRDATM(2) - GRDATM(1) )
  ELSE IF ( IDXATM .EQ. NATM ) THEN
    ALPHA = DELGRD / ( GRDATM(NATM) - GRDATM(NATM-1) ) 
  ELSE
    ALPHA = 2.0 * DELGRD / ( GRDATM(IDXATM+1) - GRDATM(IDXATM-1) )
  END IF
!
  IF ( ALPHA .LE. TOLLEV ) THEN  ! Adjust lev to match profile level
    IF ( USEHGT ) THEN
      WRNMSG = 'W-ATMLEV: Change altitude from ' // C9REAL(LEV) // ' km to' &
               // C9REAL(HGTATM(IDXATM)) // ' km to match profile.'
      LEV = HGTATM(IDXATM)
    ELSE
      WRNMSG = 'W-ATMLEV: Change pressure from ' // C9REAL(LEV) // ' hPa to' &
               // C9REAL(PREATM(IDXATM)) // ' hPa to match profile.'
      LEV = PREATM(IDXATM)
    END IF
  ELSE                           ! Insert new profile level
    CALL ADDATM ( LEV, USEHGT, IDXATM )
    WRNMSG = 'W-ATMLEV: Inserting extra profile level#' // C11INT(IDXATM)
  END IF
  CALL WRTLOG ( WRNMSG )
!
END SUBROUTINE ATMLEV
END MODULE ATMLEV_SUB
