MODULE SAVXSC_SUB
CONTAINS
SUBROUTINE SAVXSC ( IGAS )
!
! VERSION
!   19FEB20 AD Original. Based on part of old REAXSC. Checked.
! 
! DESCRIPTION    
!   Load X/S data into XSCCOM
!   Called by XSCFIL for every required X/S dataset 
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE XFLCOM_DAT ! Contents of .xsc file
    USE XSCCOM_DAT ! Tabulated Cross-Section data
!
! SUBROUTINES
    USE TRIANG_SUB ! 2-D interpolation of irregular grid using triangulation
!
  IMPLICIT NONE
!
! ARGUMENTS      
    INTEGER(I4),   INTENT(IN)  :: IGAS   ! Index of absorber in GASCOM
!
! LOCAL VARIABLES
    INTEGER(I4) :: IOF    ! Offset for each (p,T) tabulation in ABS array
    INTEGER(I4) :: IOFXFL ! Offset for abs.coeff data in ABCXFL
    INTEGER(I4) :: ISPC   ! Counter for spectral ranges
    INTEGER(I4) :: IXFL   ! Counter for tabulated datasets in XFLCOM
    INTEGER(I4) :: IXSC   ! Index of datasets in XSCCOM
    INTEGER(I4) :: IXT    ! Counter for X/S (p,T) tables
    INTEGER(I4) :: NPT    ! No. of X/S data points in current (p,T) table
    INTEGER(I4) :: NSPC   ! No. different spectral ranges in XFLCOM
    INTEGER(I4) :: NTRI   ! No. triangles found
    INTEGER(I4) :: NXP    ! No.abs.coeff datapoints in current spectral range
    INTEGER(I4) :: NXT    ! No. tabulated datasets in current spectral range
    INTEGER(I4),  ALLOCATABLE :: IDXTRI(:,:) ! Triangle coordinates
    TYPE(XSCTYP), ALLOCATABLE :: XSCSAV(:)   ! Saved XSC during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( ALLOCATED ( XSC ) ) CALL MOVE_ALLOC ( XSC, XSCSAV ) 
  IXSC = NXSC 
  NSPC = MAXVAL ( XFL%ISP ) ! No.different spectral ranges in XFL
  NXSC = NXSC + NSPC
  ALLOCATE ( XSC(NXSC) ) 
  IF ( ALLOCATED ( XSCSAV ) ) XSC(1:IXSC) = XSCSAV
!
  DO ISPC = 1, NSPC
    IXSC = IXSC + 1
    XSC(IXSC)%IGS = IGAS

    NXT = COUNT ( XFL%ISP .EQ. ISPC ) ! No.(p,T) tables
    XSC(IXSC)%NXT = NXT                                
    ALLOCATE ( XSC(IXSC)%IOF(NXT) )
    ALLOCATE ( XSC(IXSC)%NPT(NXT) )
    ALLOCATE ( XSC(IXSC)%PRE(NXT) ) 
    ALLOCATE ( XSC(IXSC)%TEM(NXT) )
    ALLOCATE ( XSC(IXSC)%WN1(NXT) )
    ALLOCATE ( XSC(IXSC)%DWN(NXT) )
!
    NXP = 0                           ! No.abs.coefficient data points
    DO IXFL = 1, NXFL
      IF ( XFL(IXFL)%ISP .EQ. ISPC ) NXP = NXP + XFL(IXFL)%NPT
    END DO
    XSC(IXSC)%NXP = NXP
    ALLOCATE ( XSC(IXSC)%ABS(NXP) )

    IOF = 0
    IXT = 0
    DO IXFL = 1, NXFL
      IF ( XFL(IXFL)%ISP .NE. ISPC ) CYCLE
      IXT = IXT + 1
      IF ( IXT .EQ. 1 ) THEN
        XSC(IXSC)%WNL = XFL(IXFL)%WN1
        XSC(IXSC)%WNU = XFL(IXFL)%WN2
      ELSE
        XSC(IXSC)%WNL = MIN ( XSC(IXSC)%WNL, XFL(IXFL)%WN1 ) 
        XSC(IXSC)%WNU = MAX ( XSC(IXSC)%WNU, XFL(IXFL)%WN2 )
      END IF 
      XSC(IXSC)%WN1(IXT) = XFL(IXFL)%WN1
      NPT = XFL(IXFL)%NPT
      XSC(IXSC)%NPT(IXT) = NPT
      XSC(IXSC)%DWN(IXT) = ( XFL(IXFL)%WN2 - XFL(IXFL)%WN1 ) / ( NPT - 1.0D0 )
      XSC(IXSC)%PRE(IXT) = XFL(IXFL)%PRE
      XSC(IXSC)%TEM(IXT) = XFL(IXFL)%TEM
      IOFXFL = XFL(IXFL)%IOF
      XSC(IXSC)%IOF(IXT) = IOF
      XSC(IXSC)%ABS(IOF+1:IOF+NPT) = ABCXFL(IOFXFL+1:IOFXFL+NPT)
      IOF = IOF + NPT
    END DO
!
    CALL TRIANG ( NXT, XSC(IXSC)%TEM, XSC(IXSC)%PRE, NTRI, IDXTRI )
    XSC(IXSC)%NTRI = NTRI
    ALLOCATE ( XSC(IXSC)%ITRI(NTRI,3) )
    XSC(IXSC)%ITRI = IDXTRI
!
  END DO
!
  DEALLOCATE ( XFL, ABCXFL )
!
END SUBROUTINE SAVXSC
END MODULE SAVXSC_SUB
