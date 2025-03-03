MODULE SPCINI_SUB
CONTAINS
SUBROUTINE SPCINI ( ISPC, MWID, FAIL, ERRMSG ) 
!
! VERSION
!   15JUL24 AD Checked.
!   23MAY20 AD Ensure WN1WID is not negative. Checked.
!   20DEC17 AD F90 version. Checked.
!
! DESCRIPTION
!   Initialise widemesh grid
!   Called by RFMSPC once for each spectral range
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE GASCOM_DAT ! Molecule and isotope data
    USE SPCCOM_DAT ! Spectral range data
    USE WIDCOM_DAT ! Widemesh data
    USE FULCOM_DAT, ONLY: WNRFUL ! Spacing of full grid data
    USE RFMCON_DAT, ONLY: FWIND  ! Window [cm-1] for widemesh calc
!
! SUBROUTINES
    USE INILBL_SUB ! Initialise line-by-line calc segments
    USE INILUT_SUB ! Initialise LUT data for each new spectral range
    USE INISTT_SUB ! Initialise widemesh statistics
    USE INISVD_SUB ! Initialise SVD data for each new spectral range
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: ISPC   ! Spectral range index
    INTEGER(I4),   INTENT(OUT) :: MWID   ! No. widemesh intervals
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL     :: FIRST = .TRUE. ! T= first call of this routine
    INTEGER(I4) :: IWD2 ! Widemesh half-interval counter (0:NWD2)
    INTEGER(I4) :: IWID ! Widemesh boundary counter (0:NWID)
    LOGICAL, SAVE, ALLOCATABLE :: HITGAS(:) ! Flags indicating HITRAN gases
    LOGICAL, SAVE, ALLOCATABLE :: XSCGAS(:) ! Flags indicating x/s gases
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DELWID = WNRFUL * NINT ( 1.0D0 / WNRFUL )
!
! Determine no. wide mesh grid points required (NB: IWID = 0,1 .... NWID so
! NWID+1 grid points, NWID intervals).
!
  WN1WID = INT ( MAX ( SPC(ISPC)%WXL, 0.0D0 ) ) ! Possible that %WXL is negative
!
! If upper limit on WNOWID boundary need to add extra WID interval to include it
! hence add 0.50001D0
  NWID   = NINT ( ( SPC(ISPC)%WXU - WN1WID ) / DELWID + 0.50001D0 ) 
  MWID   = NWID
  WN2WID = WN1WID + NWID * DELWID
  NWD2 = 2 * NWID
!
  IF ( ALLOCATED ( WNOWID ) ) DEALLOCATE ( WNOWID, WNOWD2 ) 
  ALLOCATE ( WNOWID(0:NWID) )
  ALLOCATE ( WNOWD2(0:NWD2) )
!
  DO IWID = 0, NWID
    WNOWID(IWID) = WN1WID + IWID * DELWID
  END DO
!
  DO IWD2 = 0, NWD2
    WNOWD2(IWD2) = WN1WID + 0.5 * IWD2 * DELWID
  END DO
!
! Set limits for including lines in widemesh calculations
  WNLWID = WN1WID - FWIND
  WNUWID = WN2WID + FWIND
!
! On first call save GAS%HIT and GAS%XSC flags for each gas
  IF ( FIRST ) THEN
    ALLOCATE ( HITGAS(NGAS), XSCGAS(NGAS) )
    HITGAS = GAS%HIT
    XSCGAS = GAS%XSC
    FIRST = .FALSE.
  END IF
!
  GAS%HIT = HITGAS         ! reload original values
  GAS%XSC = XSCGAS
  IF ( LUTFLG ) THEN
    CALL INILUT ( ISPC, FAIL, ERRMSG ) 
    IF ( FAIL ) RETURN
  END IF
  IF ( SVDFLG ) THEN
    CALL INISVD ( ISPC, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
  END IF
  CALL INILBL
!
  IF ( ALLOCATED ( ABSWID ) ) DEALLOCATE ( ABSWID ) 
  ALLOCATE ( ABSWID(3,NWID,NLBL) ) ; ABSWID = 0.0    ! NLBL in WIDCOM
  USECNT = NTEFLG
  IF ( USECNT ) THEN
    IF ( ALLOCATED ( CNTWID ) ) DEALLOCATE ( CNTWID ) 
    ALLOCATE ( CNTWID(3,NWID,NLBL) ) ; CNTWID = 0.0
  END IF
!
  IF ( WIDFLG ) CALL INISTT ( ISPC ) 
!
END SUBROUTINE SPCINI
END MODULE SPCINI_SUB
