MODULE SETNTE_SUB 
CONTAINS
SUBROUTINE SETNTE
!
! VERSION
!   08SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Set IUSNTE,ILSNTE in HITCOM
!   Called by REAHIT for each HITRAN line if NTEFLG enabled
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE HITCOM_DAT ! HITRAN line data
    USE NTECOM_DAT ! Non-LTE data
!    
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IDXALS ! Encoded lower state info (all isotopes)
    INTEGER(I4) :: IDXAUS ! Encoded upper state info (all isotopes)
    INTEGER(I4) :: IDXILS ! Encoded lower state info
    INTEGER(I4) :: IDXIUS ! Encoded upper state info
    INTEGER(I4) :: INTE   ! Counter for Non-LTE datasets
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  HIT%IUV = 0
  HIT%ILV = 0
!
! Match encoding in VIBCHK.FOR
  IDXAUS = HIT%IDM * 1000000 + HIT%IUS  ! non-isotopic upper level
  IDXIUS = IDXAUS + HIT%IDI * 1000      ! isotopic upper level
  IDXALS = HIT%IDM * 1000000 + HIT%ILS  ! non-isotopic lower level
  IDXILS = IDXALS + HIT%IDI * 1000      ! isotopic lower level
!
  DO INTE = 1, NNTE
    IF ( IDXIUS .EQ. NTE(INTE)%IDX .OR. &
         IDXAUS .EQ. NTE(INTE)%IDX         ) HIT%IUV = INTE
    IF ( IDXILS .EQ. NTE(INTE)%IDX .OR. &
         IDXALS .EQ. NTE(INTE)%IDX         ) HIT%ILV = INTE
  END DO
!
END SUBROUTINE SETNTE
END MODULE SETNTE_SUB
