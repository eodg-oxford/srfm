MODULE OBSCHK_SUB
CONTAINS
SUBROUTINE OBSCHK ( RVAL, USEHGT, FAIL, ERRMSG )
!
! VERSION
!   21MAY24 AD Checked.
!   29APR21 AD Bug#30 Don't set TAN%IAT = IATOBS with NADFLG, just ZENFLG
!   05AUG19 AD Add USEHGT, RVAL arguments. Change ALTOBS to HGTOBS. Checked.
!   02OCT17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Check and set observer altitude
!   Called by DRVOBS if OBS flag enabled
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE TANCOM_DAT ! Tangent path data
    USE OBSCOM_DAT, ONLY: IATOBS ! Observer location data
    USE FLGCOM_DAT, ONLY: ZENFLG ! T = nadir,zenith-viewing plane atm.
!
! SUBROUTINES
    USE ATMLEV_SUB ! Find/insert atmospheric level for given altitude
    USE C9REAL_GEN ! Write real number as C*9 string
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4),      INTENT(INOUT) :: RVAL   ! Observer altitude or pressure
    LOGICAL,       INTENT(IN)    :: USEHGT ! T=Check HGTOBS, F=Check PREOBS
    LOGICAL,       INTENT(OUT)   :: FAIL   ! T= fatal error is detected
    CHARACTER(80), INTENT(OUT)   :: ERRMSG ! Error message written if FAIL
!
! LOCAL VARIABLES 
    LOGICAL     :: LABOVE ! T=observer above top of atmosphere
    INTEGER(I4) :: IATM   ! Index of atmospheric level for observer
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .TRUE.
!
  IF ( USEHGT ) THEN
    IF ( .NOT. SETHGT ) THEN
      ERRMSG = 'F-OBSCHK: Cannot set observer altitude without *HGT profile'
    ELSE   
      LABOVE = RVAL .GT. HGTTOA  
      IF ( RVAL .LT. HGTSFC ) THEN
        ERRMSG = 'F-OBSCHK: Observer Altitude below base of atmosphere, =' // &
             C9REAL(HGTSFC) // ' [km]'              
      ELSE IF ( ZENFLG .AND. LABOVE ) THEN ! Observer > atmosphere
        ERRMSG = 'F-OBSCHK: Observer Altitude above top of atmosphere, =' // &
               C9REAL(HGTTOA) // ' [km]'              
      ELSE
        FAIL = .FALSE.
      END IF
    END IF
  ELSE
    LABOVE = RVAL .LT. PREATM(NATM)     
    IF ( RVAL .GT. PREATM(IATSFC) ) THEN
      ERRMSG = 'F-OBSCHK: Observer Pressure below base of atmosphere, =' // &
             C9REAL(PREATM(IATSFC)) // ' [hPa]'              
    ELSE IF ( ZENFLG .AND. LABOVE ) THEN ! Above atmosphere
      ERRMSG = 'F-OBSCHK: Observer Pressure above top of atmosphere, =' // &
             C9REAL(PREATM(NATM)) // ' [hPa]'              
    ELSE
      FAIL = .FALSE.
    END IF
  END IF  
!
  IF ( FAIL ) RETURN  
!
  IF ( LABOVE ) THEN
    IATOBS = NATM
  ELSE
    CALL ATMLEV ( RVAL, USEHGT, IATM )
! Note that ATMLEV itself adjusts IATOBS so cannot use this as argument
    IATOBS = IATM
  END IF
!
  IF ( ZENFLG ) TAN%IAT = IATOBS
!
END SUBROUTINE OBSCHK
END MODULE OBSCHK_SUB

