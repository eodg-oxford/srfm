MODULE CHKLEV_SUB
CONTAINS
SUBROUTINE CHKLEV
!
! VERSION
!   25FEB23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Set up tangent paths for intermediate output levels
!   Called by DRVCHK
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE LEVCOM_DAT ! Intermediate output levels 
    USE TANCOM_DAT ! Tangent path data
    USE ATMCOM_DAT, ONLY: NATM   ! No. Atmospheric levels
    USE FLGCOM_DAT, ONLY: ZENFLG ! T=zenith-viewing
    USE OBSCOM_DAT, ONLY: IATOBS ! Observer location data
!
! SUBROUTINES
    USE ADDTAN_SUB ! Add tangent ray path
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE 
!
! LOCAL VARIABLES
    INTEGER(I4) :: ILEV ! Intermediate output level counter
    INTEGER(I4) :: ILOW ! Lower atmos profile level for tangent path
    INTEGER(I4) :: ITAN ! Counter for output tangent heights
    INTEGER(I4) :: IUPP ! Upper atmos profile level for tangent path
    INTEGER(I4) :: LLEV ! No output levels selected (input value of NLEV)
    INTEGER(I4) :: LTAN ! Saved value of MTAN on input
    TYPE(LEVTYP), ALLOCATABLE :: LEVSAV(:) ! Saved LEV during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!    
  LLEV = NLEV  ! No output levels 
!
! Set paths 1:NLEV for upward viewing, ie downward integration from TOA.
  LEV(1:NLEV)%IDR = -1
!
! Unless ZEN flag, duplicate LEV for downward viewing, upward integ. paths
  IF ( .NOT. ZENFLG ) THEN
    CALL MOVE_ALLOC ( LEV, LEVSAV ) 
    NLEV = 2 * NLEV
    ALLOCATE ( LEV(NLEV) ) 
    LEV(1:LLEV) = LEVSAV
    LEV(LLEV+1:NLEV) = LEVSAV
    LEV(LLEV+1:NLEV)%IDR = +1
  END IF
!
  ALLOCATE ( ITNLEV(MTAN,NLEV) )
  ITNLEV = 0
!
  LTAN = MTAN
  IUPP = NATM
  IF ( IATOBS .GT. 0 ) IUPP = IATOBS
  DO ITAN = 1, MTAN
    ILOW = TAN(ITAN)%IAT
! Downward path
    DO ILEV = LLEV, 1, -1           ! Profile output levels
      IF ( LEV(ILEV)%IAT .LT. ILOW ) EXIT
      CALL ADDTAN ( ITAN, .FALSE. ) 
      ITNLEV(ITAN,ILEV) = MTAN
    END DO
    IF ( ZENFLG ) CYCLE
! Upward path
    DO ILEV = 1, LLEV
      IF ( LEV(ILEV)%IAT .LT. ILOW ) CYCLE
      IF ( LEV(ILEV)%IAT .GT. IUPP ) EXIT
      CALL ADDTAN ( ITAN, .FALSE. ) 
      ITNLEV(ITAN,ILEV+LLEV) = MTAN
    END DO
  END DO
!
  CALL WRTLOG ( 'I-CHKLEV: No.extra tan. paths reqd for intermediate ' // &
                'Lev output=' // C11INT ( MTAN - LTAN ) ) 
!
END SUBROUTINE CHKLEV
END MODULE CHKLEV_SUB

