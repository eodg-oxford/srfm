MODULE SETIFP_SUB
CONTAINS
SUBROUTINE SETIFP ( LUNHIT, IREC, IFP, FAIL, ERRMSG )
!
! VERSION
!   26AUG24 AD Checked.
!   11AUG23 AD Add IFP as argument. Remove RECBIN.
!   23MAY23 AD Original. Was part of INIHFL.
!
! DESCRIPTION    
!   Set forward pointers for HITRAN binary file
!   Called by BINFIL and INIBIN
!   Sets the forward pointers for all molecules in file.
!   Intention is to set all FPs to just *before* IREC is read in,
!   so if IREC is a line parameter file it will also be one of the IFP values.
!
! VARIABLE KINDS
    USE KIND_DAT  
!
  IMPLICIT NONE
!
! ARGUMENTS      
    INTEGER(I4),   INTENT(IN)  :: LUNHIT ! LUN for HITRAN binary file
    INTEGER(I4),   INTENT(IN)  :: IREC   ! Record# in file
    INTEGER(I4),   INTENT(OUT) :: IFP(:) ! List of forward pointers
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: LSTFWD = -7 ! LSTAT value for Fwd Ptr record
    INTEGER(I4), PARAMETER :: LSTLIN = 10 ! LSTAT value for line param. record
    INTEGER(I4), PARAMETER :: MAXPTR = 56 ! No. fwd pointers  in HIT .bin file
    INTEGER(I4), PARAMETER :: NFPBIN = 14 ! No. fwd pointers per FP record
!
! LOCAL VARIABLES
    INTEGER(I4)   :: IDM    ! HITRAN Molecule ID
    INTEGER(I4)   :: IDMLOC ! Subset of FPs within FP record
    INTEGER(I4)   :: IDMOFF ! Molec# of 1st molec in IDMLOC list
    INTEGER(I4)   :: IDUMMY ! Dummy integer
    INTEGER(I4)   :: IFWDPT ! Rel.Fwd pointer for next molecule of same type
    INTEGER(I4)   :: IOS    ! Value of IOSTAT for error messages
    INTEGER(I4)   :: JREC   ! Record# in file 
    INTEGER(I4)   :: LSTAT  ! Record type
    REAL(R8)      :: WNO    ! Wavenumber [cm-1] of record just read (dummy)
    CHARACTER(75) :: CDUM75 ! Dummy padding to skip to IFWDPT
    INTEGER(I4)   :: IFPLOC(NFPBIN) ! Forward pointers from HITRAN file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
  IFP = 0
!
! Establish type of record of IREC 
  JREC = IREC
  READ ( LUNHIT, REC=JREC, IOSTAT=IOS, ERR=900 ) LSTAT
!
! Step backwards until last record of a FP block is reached
  DO WHILE ( LSTAT .NE. LSTFWD ) 
    JREC = JREC - 1
    READ ( LUNHIT, REC=JREC, IOSTAT=IOS, ERR=900 ) LSTAT
  END DO
!
! Continue backwards until first record of FP block is reached
  DO WHILE ( LSTAT .EQ. LSTFWD ) 
    JREC = JREC - 1
    READ ( LUNHIT, REC=JREC, IOSTAT=IOS, ERR=900 ) LSTAT
  END DO
! 
! Read forward through FP block setting pointers
  DO 
    JREC = JREC + 1  
    READ ( LUNHIT, REC=JREC, ERR=900, IOSTAT=IOS ) &
           LSTAT, IDMOFF, IDUMMY, WNO, IFPLOC
    IF ( LSTAT .NE. LSTFWD ) EXIT  ! JREC points to first line param record
    DO IDMLOC = 1, NFPBIN
      IDM = IDMOFF + IDMLOC - 1
      IFP(IDM) = JREC + IFPLOC(IDMLOC)
    END DO
  END DO
!
! Read forward through line parameter records up to IREC updating pointers
  DO WHILE ( JREC .LT. IREC ) 
    READ ( LUNHIT, REC=JREC, IOSTAT=IOS, ERR=900 ) LSTAT, IDM, CDUM75, IFWDPT
    IFP(IDM) = JREC + IFWDPT
    JREC = JREC + 1
  END DO
!
900 CONTINUE
  FAIL = IOS .NE. 0
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-SETIFP: Failed to read HITRAN file, rec#:', IREC, '. IOSTAT=', IOS
!
END SUBROUTINE SETIFP
END MODULE SETIFP_SUB
