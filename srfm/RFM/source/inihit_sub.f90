MODULE INIHIT_SUB
CONTAINS
SUBROUTINE INIHIT ( WNOREQ, WNOUPP, FAIL, ERRMSG )
!
! VERSION
!   20AUG24 AD Checked.
!   11AUG23 AD Restructured for multiple files.
!   22MAY23 AD Original. Replaces INIHFL.
!
! DESCRIPTION    
!   Initialise the HITRAN line data file(s)
!   Called by SPCWID and SPCFIN once for each spectral range.
!   Set the forward pointers for each required gas
!   and position ready to read required wavenumber as next record
!
! VARIABLE KINDS
    USE KIND_DAT  
!
! GLOBAL DATA
    USE HFLCOM_DAT ! HITRAN file data
!
! SUBROUTINES
    USE INIBIN_SUB ! Initialise the HITRAN line data binary file
    USE INIHDB_SUB ! Initialise pointer in HITRAN line database file
    USE INIPAR_SUB ! Initialise pointer in HITRAN line data .par file
    USE RECBIN_SUB ! Read record from HITRAN line data binary file
    USE RECHDB_SUB ! Read record from HITRAN database file
    USE RECPAR_SUB ! Read record from HITRAN line data .par file
!
  IMPLICIT NONE
!
! ARGUMENTS      
    REAL(R8),      INTENT(IN)  :: WNOREQ ! Initial wavenumber
    REAL(R8),      INTENT(IN)  :: WNOUPP ! Maximum wavenumber
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    LOGICAL     :: EOF    ! T=end-of-file (dummy)
    INTEGER(I4) :: IHFL   ! Counter for HITRAN data files
    INTEGER(I4) :: LUNHIT ! LUN for reading file data
    LOGICAL     :: USEIDM(MAXIDM) ! T=use molec# from file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Set upper limit on required wavenumber to avoid searching beyond this
  IF ( WNOUPP .GT. 0.0D0 ) WNOMAX = WNOUPP
!
  DO IHFL = 1, NHFL
    LUNHIT = HFL(IHFL)%LUN
    SELECT CASE ( HFL(IHFL)%TYP )
      CASE ( 'BIN' ) 
        CALL INIBIN ( LUNHIT, WNOREQ, HFL(IHFL)%IFP, FAIL, ERRMSG )
      CASE ( 'HDB' ) 
        CALL INIHDB ( LUNHIT, WNOREQ, HFL(IHFL)%HDB%WNOFMT, FAIL, ERRMSG )
      CASE ( 'PAR' ) 
        CALL INIPAR ( LUNHIT, WNOREQ, FAIL, ERRMSG )
      CASE DEFAULT; STOP 'F-INIHIT: Logical error'
    END SELECT
    IF ( FAIL ) RETURN
  END DO
!
! If more than one file, read next record from each file into HFL%HIT buffer
  IF ( NHFL .GT. 1 ) THEN
    DO IHFL = 1, NHFL
      LUNHIT = HFL(IHFL)%LUN
      USEIDM = HFL(IHFL)%IDM
      SELECT CASE ( HFL(IHFL)%TYP )
        CASE ( 'BIN' ) 
          CALL RECBIN ( LUNHIT, HFL(IHFL)%HIT, USEIDM, HFL(IHFL)%IFP, &
                        EOF, FAIL, ERRMSG )
        CASE ( 'HDB' ) 
          CALL RECHDB ( LUNHIT, HFL(IHFL)%HIT, USEIDM, WNOMAX, HFL(IHFL)%HDB, &
                        EOF, FAIL, ERRMSG )
        CASE ( 'PAR' ) 
          CALL RECPAR ( LUNHIT, HFL(IHFL)%HIT, USEIDM, WNOMAX, & 
                        EOF, FAIL, ERRMSG )
        CASE DEFAULT; STOP 'F-INIHIT: Logical error'
      END SELECT
      IF ( FAIL ) RETURN
    END DO
  END IF
!
END SUBROUTINE INIHIT
END MODULE INIHIT_SUB
