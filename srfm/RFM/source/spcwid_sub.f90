MODULE SPCWID_SUB 
CONTAINS
SUBROUTINE SPCWID ( FAIL, ERRMSG ) 
!
! VERSION
!   26AUG24 AD Checked.
!   11AUG23 AD Remove test for OPNHFL. Checked.
!   22MAY23 AD Change INIHFL to INIHIT and add WNUWID argument
!   29JAN20 AD Check for EOF immediately after REAHIT. Checked.
!   05MAY18 AD Bug#6: exclude any line WNUM=WNUWID
!   01MAY17 AD F90 conversion of rfmwid.for. Checked.
!
! DESCRIPTION  
!   Perform wide mesh absorption calculations.
!   Called by RFMSPC once for each spectral range
!   Make a wide pass over the spectrum frequency grid. Each line is read in 
!   turn and the appropriate paths identified. The absorption in freqency 
!   interval some distance off due to the line wings is calculated at three
!   points within the interval for future quadratic interpolation. 
!   Whether or not the line wing absorption in a given frequency interval is
!   treated in this way is determined by the frequency windows in operation. 
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE GASCOM_DAT ! Molecule and isotope data
    USE HITCOM_DAT ! HITRAN line data
    USE WIDCOM_DAT ! Widemesh data
    USE FLGCOM_DAT, ONLY: WIDFLG ! T = output Wide Mesh Calc diagnostics
    USE IDXCON_DAT, ONLY: IDXH2O ! RFM/HITRAN index for H2O
    USE RFMCON_DAT, ONLY: FWIND, FEXC ! wide,fine mesh limits
!
! SUBROUTINES
    USE ADJUST_SUB ! Adjust line parameters for path conditions
    USE INIHIT_SUB ! Initialise the HITRAN line data file(s)
    USE LINSHP_SUB ! Apply spectral lineshape
    USE REAHIT_SUB ! Read record from HITRAN line data file
    USE WIDSTT_SUB ! Add to line statistics
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL     :: EOF         ! T = end-of-file found
    LOGICAL     :: SUBWNG      ! T = Subtract abs.coeff at 25cm-1
    INTEGER(I4) :: ICLC        ! Index of calc paths
    INTEGER(I4) :: IGAS        ! Absorber counter
    INTEGER(I4) :: ILBL        ! Counter for line-by-line calculated paths
    INTEGER(I4) :: IQAD        ! Parabolic point counter for each WM interval
    INTEGER(I4) :: ISHP        ! Lineshape index for gas in path
    INTEGER(I4) :: IWD2        ! Half-Mesh grid point index
    INTEGER(I4) :: IWD2L,IWD2U ! Low/Upp Half-Mesh points incl.line contrib.
    INTEGER(I4) :: JEXCL,JEXCU ! Closest Low/Upp Wide Mesh Intvl including line 
    INTEGER(I4) :: JWID        ! Wide mesh interval counter
    INTEGER(I4) :: JWIDL,JWIDU ! Low/Upp Wide Mesh Intvls incl.line contrib.
    REAL(R4)    :: ANTE        ! Non-lte factor for k abs
    REAL(R4)    :: CNTE        ! Non-lte factor for c abs
    REAL(R8)    :: WNUM        ! Wavenumber of HITRAN line
    REAL(R4)    :: ABSORP(0:NWD2) ! Absorption
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( NLBL .EQ. 0 )  RETURN                 ! No LBL calcs required
!  
  CALL INIHIT ( WNLWID, WNUWID, FAIL, ERRMSG )
  IF ( FAIL ) RETURN
!
  DO 
    CALL REAHIT ( EOF, FAIL, ERRMSG )       ! Get next line from file
    IF ( FAIL ) RETURN
    IF ( EOF ) EXIT
    IGAS = HIT%IGS
    WNUM = HIT%WNO
    IF ( WNUM .GE. WNUWID ) EXIT
    JWIDL = MAX (    1, INT( 1.0D0+(WNUM-FWIND-WN1WID)/DELWID )  )
    JWIDU = MIN ( NWID, INT( 1.0D0+(WNUM+FWIND-WN1WID)/DELWID )  )
    JEXCL = MIN ( NWID, INT( 1.0D0+(WNUM-FEXC-WN1WID)/DELWID ) -1)
    JEXCU = MAX (    1, INT( 1.0D0+(WNUM+FEXC-WN1WID)/DELWID ) +1)
    IF ( WIDFLG ) CALL WIDSTT ( WNUM, IGAS, JWIDL, JEXCL, JEXCU, JWIDU ) 
    IWD2L = 2 * ( JWIDL - 1 )
    IWD2U = 2 * JWIDU 
    DO ILBL = 1, NLBL                         ! Loop over LBL calc paths 
      ICLC = IDXLBL(ILBL)
      IF ( CLC(ICLC)%IGS .NE. IGAS ) CYCLE
      CALL ADJUST ( CLC(ICLC)%TEM, CLC(ICLC)%PRE, CLC(ICLC)%PPA, &
                    CLC(ICLC)%AMT, ANTE, CNTE )
! Calculate path adjusted linshape parameters 
      SUBWNG = SUBH2O .AND. HIT%IDM .EQ. IDXH2O
      ISHP = GAS(IGAS)%SHP
! Loop over coarse freq. intervals
      CALL LINSHP ( ISHP, WNOWD2(IWD2L:IWD2U), ABSORP(IWD2L:IWD2U), SUBWNG )
!
! Set absorption=0.0 beyond 25cm-1 from line centre 
      IF ( WNOWD2(IWD2L) .LT. WNUM-FWIND )   ABSORP(IWD2L) = 0.0
      IF ( WNOWD2(IWD2L+1) .LT. WNUM-FWIND ) ABSORP(IWD2L+1) = 0.0
      IF ( WNOWD2(IWD2U-1) .GT. WNUM+FWIND ) ABSORP(IWD2U-1) = 0.0
      IF ( WNOWD2(IWD2U) .GT. WNUM+FWIND )   ABSORP(IWD2U) = 0.0
!
      DO JWID = JWIDL, JWIDU
        IF ( JWID .GT. JEXCL .AND. JWID .LT. JEXCU ) CYCLE 
        DO IQAD = 1, 3
          IWD2 = 2 * JWID + IQAD - 3
          IF ( GAS(IGAS)%NTE ) THEN
            ABSWID(IQAD,JWID,ILBL) = ABSWID(IQAD,JWID,ILBL) + &
                                     ANTE * ABSORP(IWD2) 
            CNTWID(IQAD,JWID,ILBL) = CNTWID(IQAD,JWID,ILBL) + &
                                     CNTE * ABSORP(IWD2) 
          ELSE
            ABSWID(IQAD,JWID,ILBL) = ABSWID(IQAD,JWID,ILBL) + ABSORP(IWD2) 
          END IF
        END DO
      END DO
    END DO          ! End loop over paths
  END DO
!
END SUBROUTINE SPCWID
END MODULE SPCWID_SUB
