MODULE SPCFIN_SUB
CONTAINS
SUBROUTINE SPCFIN ( NEWSPC, FAIL, ERRMSG )
!
! VERSION
!   18DEC23 AD Bug#43 Change IWID argument to NEWSPC
!   11AUG23 AD Change INIHFL to INIHIT and add WNUWID argument
!   01MAY17 AD F90 conversion of F77 module RFMFIN. Checked.
!
! DESCRIPTION
!   Perform a fine pass over the wavenumber grid
!   Called by RFMSPC for each wide-mesh interval.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE FINCOM_DAT ! Finemesh data
    USE GASCOM_DAT ! Molecule and isotope data
    USE HITCOM_DAT ! HITRAN line data
    USE IDXCON_DAT, ONLY: IDXH2O ! RFM/HITRAN index for H2O
!
! SUBROUTINES
    USE ADJUST_SUB ! Adjust line parameters for path conditions
    USE INIHIT_SUB ! Initialise the HITRAN line data file(s)
    USE LINSHP_SUB ! Apply spectral lineshape
    USE REACYC_SUB ! Read HITRAN data into cyclic buffers
!
  IMPLICIT NONE
!                  
! ARGUMENTS
    LOGICAL,       INTENT(INOUT) :: NEWSPC ! T = start of new spectral range
    LOGICAL,       INTENT(OUT)   :: FAIL   ! Set TRUE if fatal error detected
    CHARACTER(80), INTENT(OUT)   :: ERRMSG ! Error message if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: ICLC   ! Counter for calc paths
    INTEGER(I4) :: ICYC   ! Cyclic buffer index
    INTEGER(I4) :: IGAS   ! Absorber counter
    INTEGER(I4) :: ILIN   ! Line counter
    INTEGER(I4) :: ISHP   ! Lineshape code for gas in path
    LOGICAL     :: SUBWNG ! T = Subtract abs.coeff at 25cm-1
    REAL(R4)    :: ANTE   ! Non-lte factor for k abs
    REAL(R4)    :: CNTE   ! Non-lte factor for c abs
    REAL(R4)    :: ABSLIN(NFIN) ! Single line absorption for current path
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
  IF ( .NOT. ANY ( CLC%LBL ) ) RETURN   ! No line-by-line calcs required
!
! Read lines between WNLFIN and WNUFIN into cyclic buffer
  IF ( NEWSPC ) THEN        
    NLIN = 0                ! Clear cyclic buffer
    ICYC1 = 1
    CALL INIHIT ( WNLFIN, 0.0D0, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    NEWSPC = .FALSE.
  END IF
  CALL REACYC ( WNLFIN, WNUFIN, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN
!
! Loop over paths
  DO ICLC = 1, NCLC
    IF ( .NOT. CLC(ICLC)%LBL ) CYCLE   ! not a line-by-line calc molecule
    IGAS = CLC(ICLC)%IGS
    ISHP = GAS(IGAS)%SHP
    SUBWNG = SUBH2O .AND. GAS(IGAS)%IDM .EQ. IDXH2O
!
! Loop over lines stored in buffer 
    DO ILIN = 1, NLIN
      ICYC = MOD ( ILIN + ICYC1 - 2, NCYC ) + 1
      IF ( CYC(ICYC)%IGS .NE. IGAS ) CYCLE
! NB < WNUFIN allows for any line sitting exactly on a widemesh boundary
! at WNUFIN will already be included in widemesh calc so .exclude from finemesh
      IF ( CYC(ICYC)%WNO .GE. WNUFIN ) CYCLE
      HIT = CYC(ICYC)
      CALL ADJUST ( CLC(ICLC)%TEM, CLC(ICLC)%PRE, CLC(ICLC)%PPA, &
                    CLC(ICLC)%AMT, ANTE, CNTE )
! Apply lineshape
      CALL LINSHP ( ISHP, WNOFIN, ABSLIN, SUBWNG ) 
!
! NB without non-LTE, ANTE=CNTE=1
      IF ( GAS(IGAS)%NTE ) THEN
        ABSFIN(1:NFIN,ICLC) = ABSFIN(1:NFIN,ICLC) + ANTE * ABSLIN 
        CNTFIN(1:NFIN,ICLC) = CNTFIN(1:NFIN,ICLC) + CNTE * ABSLIN
      ELSE
        ABSFIN(1:NFIN,ICLC) = ABSFIN(1:NFIN,ICLC) + ABSLIN
      END IF
!
    END DO  ! end loop over lines
!
  END DO  ! end loop over paths
!
END SUBROUTINE SPCFIN
END MODULE SPCFIN_SUB

