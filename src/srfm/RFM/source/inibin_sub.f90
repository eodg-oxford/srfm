MODULE INIBIN_SUB
CONTAINS
SUBROUTINE INIBIN ( LUNHIT, WNOREQ, IFP, FAIL, ERRMSG )
!
! VERSION
!   19AUG24 AD Checked.
!   11AUG23 AD Add IFP argument, remove BINCOM.
!   22MAY23 AD Original. Taken from part of INIHFL.
!
! DESCRIPTION    
!   Initialise the HITRAN line data binary file
!   Called by INIHIT for each spectral range.
!   Set the forward pointers for all molecules in the file,
!   based on last record before WNOREQ
!
! VARIABLE KINDS
    USE KIND_DAT  
!
! SUBROUTINES
    USE SETIFP_SUB ! Set forward pointers for HITRAN binary file
!
  IMPLICIT NONE
!
! ARGUMENTS      
    INTEGER(I4),   INTENT(IN)  :: LUNHIT ! Logical Unit Number
    REAL(R8),      INTENT(IN)  :: WNOREQ ! Initial wavenumber
    INTEGER(I4),   INTENT(OUT) :: IFP(:) ! Forward Pointers for all molecs
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)  :: IDUM  ! Dummy read
    INTEGER(I4)  :: IOS   ! Value of IOSTAT for error messages
    INTEGER(I4)  :: IREC  ! Record#
    INTEGER(I4)  :: JREC  ! Lower Record# for search
    INTEGER(I4)  :: KREC  ! Upper Record# for search
    REAL(R8)     :: WNO   ! Wavenumber of HITRAN record
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FAIL = .FALSE.
!
! 1st record of binary file contains LSAT,IFMT then 1st/last data rec#
  READ ( LUNHIT, REC=1, IOSTAT=IOS, ERR=900 ) IDUM, IDUM, JREC, KREC
!
! Binary search for 1st record: follows code in the IBM SLUP routine:
! K is first record, L is last. Only need to read WNO for search
  DO WHILE ( JREC + 1 .LT. KREC )       
    IREC = (JREC + KREC) / 2
    READ ( LUNHIT, REC=IREC, IOSTAT=IOS, ERR=900 ) IDUM, IDUM, IDUM, WNO
    IF ( WNO .LT. WNOREQ ) THEN   ! Reqd WN higher, try half way up to L
      JREC = IREC 
    ELSE                     ! Reqd WN lower or =, try half way down to K
      KREC = IREC 
    ENDIF
  END DO
!
! If K & L differ by 1, have found required location where the K .LT. DWNLOW 
! and L .GE. DWNLOW, Choose record L (note K<L). If there is more than 1 
! record at exactly DWNLOW, will finish pointing to first. 
! Forward pointer block are labelled with wavenumber of next line to exploit 
! this. 
!
  IREC = KREC
! Now set the initial forward pointers for each gas
  CALL SETIFP ( LUNHIT, IREC, IFP, FAIL, ERRMSG ) 
  IF ( FAIL ) RETURN
!
900 CONTINUE
  FAIL = IOS .NE. 0
  IF ( FAIL ) WRITE ( ERRMSG, * ) &
    'F-INIBIN: Failed to read HITRAN file, rec#:', IREC, '. IOSTAT=', IOS
!
END SUBROUTINE INIBIN
END MODULE INIBIN_SUB
