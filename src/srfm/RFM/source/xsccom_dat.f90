MODULE XSCCOM_DAT
!
! VERSION
!   04APR24 AD Checked.
!   08FEB19 AD Add NTRI.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Tabulated Cross-Section data.
!   The RFM stores all required .xsc datasets simultaneously, read in during
!   reading of the *XSC section of the driver table.
!   Loaded by REAXSC, accessed in SPCXSC. Pointers are deallocated in RFMDAL
!   and via XSCCOM_RESET when the library is reused.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
  PUBLIC :: XSCCOM_RESET
!
  TYPE :: XSCTYP
    INTEGER(I4)          :: IGS       ! Absorber index for XSC file
    INTEGER(I4)          :: NXP       ! Total no. abs.coeff values in dataset
    INTEGER(I4)          :: NXT       ! No.different (p,T) tabulations
    INTEGER(I4)          :: NTRI      ! No. triangles
    REAL(R8)             :: WNL       ! Lower Wno limit [cm-1] of dataset
    REAL(R8)             :: WNU       ! Upper Wno limit [cm-1] of dataset
    INTEGER(I4), POINTER :: IOF(:)    ! [NXT] Offset for each (p,T) tabulation
    INTEGER(I4), POINTER :: ITRI(:,:) ! [NTRI,3] Vertices of triangulated data
    INTEGER(I4), POINTER :: NPT(:)    ! [NXT] No.pts in each (p,T) tabulation
    REAL(R4),    POINTER :: ABS(:)    ! [NXP] Tabulated abs.coeffs [cm^2/molec]
    REAL(R4),    POINTER :: PRE(:)    ! [NXT] p [Torr] of each (p,T) tabulation
    REAL(R4),    POINTER :: TEM(:)    ! [NXT] Temp [K] of each (p,T) tabulation
    REAL(R8),    POINTER :: DWN(:)    ! [NXT] Wno incr. [cm-1] of each (p,T)
    REAL(R8),    POINTER :: WN1(:)    ! [NXT] Lower Wno. [cm-1] of each (p,T)
  END TYPE XSCTYP
!
! GLOBAL VARIABLES
    TYPE(XSCTYP), ALLOCATABLE :: XSC(:)
!
    INTEGER(I4) :: NXSC = 0 ! No. of X/S datasets being used
!
CONTAINS

  SUBROUTINE XSCCOM_RESET()
    INTEGER(I4) :: I

    IF ( ALLOCATED ( XSC ) ) THEN
      DO I = 1, SIZE ( XSC )
        IF ( ASSOCIATED ( XSC(I)%IOF ) )  DEALLOCATE ( XSC(I)%IOF )
        IF ( ASSOCIATED ( XSC(I)%ITRI ) ) DEALLOCATE ( XSC(I)%ITRI )
        IF ( ASSOCIATED ( XSC(I)%NPT ) )  DEALLOCATE ( XSC(I)%NPT )
        IF ( ASSOCIATED ( XSC(I)%ABS ) )  DEALLOCATE ( XSC(I)%ABS )
        IF ( ASSOCIATED ( XSC(I)%PRE ) )  DEALLOCATE ( XSC(I)%PRE )
        IF ( ASSOCIATED ( XSC(I)%TEM ) )  DEALLOCATE ( XSC(I)%TEM )
        IF ( ASSOCIATED ( XSC(I)%DWN ) )  DEALLOCATE ( XSC(I)%DWN )
        IF ( ASSOCIATED ( XSC(I)%WN1 ) )  DEALLOCATE ( XSC(I)%WN1 )
      END DO
      DEALLOCATE ( XSC )
    END IF

    NXSC = 0
  END SUBROUTINE XSCCOM_RESET

END MODULE XSCCOM_DAT
