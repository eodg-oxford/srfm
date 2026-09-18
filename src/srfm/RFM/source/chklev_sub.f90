MODULE CHKLEV_SUB
CONTAINS
SUBROUTINE CHKLEV
!
! VERSION
!   17SEP26 AK Preserve levels at the ray endpoint for Python output capture.
!   05JUL26 AD Checked.
!   01AUG25 AD Rewritten.
!   20FEB25 AD Checked.
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
    USE ATMCOM_DAT, ONLY: IATSFC, HGTSFC   ! No. Atmospheric levels
    USE FLGCOM_DAT, ONLY: ZENFLG ! T=zenith-viewing
!
! SUBROUTINES
    USE ADDTAN_SUB ! Add tangent ray path
    USE C9REAL_GEN ! Write real number as C*9 string
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE 
!
! LOCAL VARIABLES
    INTEGER(I4) :: ILEV ! Intermediate output level counter
    INTEGER(I4) :: ITAN ! Counter for output tangent heights
    INTEGER(I4) :: JLEV ! index of lowest output level that is retained
    INTEGER(I4) :: LLEV ! No output levels selected (input value of NLEV)
    TYPE(LEVTYP), ALLOCATABLE :: LEVSAV(:) ! Saved LEV during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!    
! Check that surface hasn't been moved above lowest output level
  JLEV = 1   ! index of lowest output level that is retained
  DO ILEV = 1, NLEV
    IF ( IATSFC .GT. LEV(ILEV)%IAT ) THEN
      CALL WRTLOG ( 'W-CHKLEV: Removing output level at ' &
                    // TRIM ( C9REAL(LEV(ILEV)%HGT) ) &
                    // ' km: below new surface level' )
      JLEV = JLEV + 1
    END IF
  END DO
!
! If bottom of atmosphere defined as output level
  IF ( BOALEV .AND. JLEV .GT. 1 ) THEN
    IF ( LEV(JLEV)%IAT .GT. IATSFC ) THEN ! need to insert new boa level
      JLEV = JLEV-1
      CALL WRTLOG ( 'W-CHKLEV: Reset BOA level to ' &
                    // TRIM ( C9REAL(HGTSFC) ) // ' km' )  
      LEV(JLEV)%HGT = HGTSFC
      LEV(JLEV)%STR = 'boa_' 
      LEV(JLEV)%IAT = IATSFC 
    END IF
  END IF

  LLEV = NLEV-JLEV+1  ! No. output levels 
  CALL MOVE_ALLOC ( LEV, LEVSAV ) 
!
  IF ( ZENFLG ) THEN  ! Only downward paths required
    NLEV = LLEV 
  ELSE
    NLEV = 2 * LLEV
  END IF
  ALLOCATE ( LEV(NLEV) ) 

! Set paths 1:LLEV for upward viewing, ie downward integration from TOA.
  LEV(1:LLEV)%IDR = -1
  LEV(1:LLEV)%IAT = LEVSAV(JLEV:)%IAT
  LEV(1:LLEV)%HGT = LEVSAV(JLEV:)%HGT
  LEV(1:LLEV)%STR = LEVSAV(JLEV:)%STR
  IF ( .NOT. ZENFLG ) THEN
    LEV(LLEV+1:NLEV)%IDR = 1
    LEV(LLEV+1:NLEV)%IAT = LEV(LLEV:1:-1)%IAT
    LEV(LLEV+1:NLEV)%HGT = LEV(LLEV:1:-1)%HGT
    LEV(LLEV+1:NLEV)%STR = LEV(LLEV:1:-1)%STR
  END IF
!
  ALLOCATE ( ITNLEV(NTAN,NLEV) )
  ITNLEV = 0
!
! Retain levels at the ray endpoint as well as levels above it. RADLEV and the
! Python capture path require a non-zero ITNLEV entry for those endpoint levels.
  DO ITAN = 1, NTAN
    DO ILEV = 1, NLEV 
      IF ( LEV(ILEV)%IAT .LT. TAN(ITAN)%IAT ) CYCLE
      CALL ADDTAN ( ITAN, .FALSE. ) 
      ITNLEV(ITAN,ILEV) = MTAN
    END DO
  END DO
!
  CALL WRTLOG ( 'I-CHKLEV: No.extra tan. paths reqd for intermediate ' // &
                'Lev output=' // C11INT ( MTAN - NTAN ) ) 
!
END SUBROUTINE CHKLEV
END MODULE CHKLEV_SUB
