MODULE INISTT_SUB
CONTAINS
SUBROUTINE INISTT ( ISPC ) 
!
! VERSION
!   01NOV23 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Initialise widemesh statistics
!   Called by SPCINI at start of each spectral range.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE SPCCOM_DAT ! Spectral range data
    USE STTCOM_DAT ! Widemesh statistics
    USE WIDCOM_DAT ! Widemesh data
    USE GASCOM_DAT, ONLY: NGAS ! No. different molecules
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: ISPC ! Spectral range index
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NSTW = NWID
  NSTG = NGAS
  WN1STT = WN1WID
  WNDSTT = DELWID
!
  IF ( ALLOCATED(STT) ) DEALLOCATE ( STT ) 
  ALLOCATE ( STT(NSTW,NSTG) ) 
  STT%ICT = 0
  STT%IXS = 0
  STT%NLF = 0
  STT%NLI = 0
  STT%NLW = 0
  STT%NQD = 0
!
! No. gases remains constant so no need to reallocate for each spectral range
  IF ( .NOT. ALLOCATED(NEXSTT) ) ALLOCATE ( NEXSTT(NSTG), NINSTT(NSTG) )
  NEXSTT = 0
  NINSTT = 0 
  WNLSTT = SPC(ISPC)%WNL
  WNUSTT = SPC(ISPC)%WNU
!
END SUBROUTINE INISTT
END MODULE INISTT_SUB
