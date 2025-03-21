MODULE ADDAIR_SUB
CONTAINS
SUBROUTINE ADDAIR
!
! VERSION
!   25JUL24 AD Checked.
!   24JUN19 AD Original. Checked.
!
! DESCRIPTION
!   Add air as absorber and atmospheric profile
!   Called by CIAFIL and/or DRVCHK.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE FLGCOM_DAT, ONLY: GRAFLG ! T=2D atmosphere
    USE IDXCON_DAT, ONLY: IDXAIR ! RFM 'molecular' index for air
!
! SUBROUTINES
    USE ADDGAS_SUB ! Add new molecule/isotope to list of absorbers
    USE ADDVMR_SUB ! Add extra VMR profile to ATMCOM
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
    USE PRFGRA_SUB ! Copy atm profile from 1D to 2D field
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    LOGICAL     :: FIRST = .TRUE. ! T= first call of this subroutine
    LOGICAL     :: NEWGAS         ! T=new molecule added, F=already stored
    INTEGER(I4) :: IVMR           ! Index of 'air' in VMR profiles
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( .NOT. FIRST ) RETURN  ! module has already been used to load 'air'
  FIRST = .FALSE.
!
  CALL ADDGAS ( IDXAIR, 0, NEWGAS ) 
  IF ( .NOT. NEWGAS ) STOP 'F-ADDAIR: Logical error#1'
!
! If a new gas, then 'air' will be last element NGAS, NVMR in GASCOM, ATMCOM 
  CALL ADDVMR
  IVMR = IDXGAS ( IDXAIR ) 
  IF ( IVMR .EQ. 0 ) STOP 'F-ADDAIR: Logical error#2'
!
  VMRATM(:,IVMR) = 1.0E6 ! 'air' mixing ratio = 1ppv.
  SETVMR(IVMR) = .TRUE.
!
! If using horizontal gradients, ensure that Air VMR defined at one (arbitrary)
! location 
  IF ( GRAFLG ) CALL PRFGRA ( 'VMR', IVMR ) 
!
END SUBROUTINE ADDAIR
END MODULE ADDAIR_SUB

