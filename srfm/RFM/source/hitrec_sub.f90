MODULE HITREC_SUB
CONTAINS
SUBROUTINE HITREC ( HIT, EOF, FAIL, ERRMSG ) 
!
! VERSION
!   14AUG24 Checked.
!   11AUG23 Rewritten
!
! DESCRIPTION
!   Read record from HITRAN line data file
!   Called by REAHIT.
!   If only a single HITRAN file is being accessed the line data are read
!   directly into the HIT structure. 
!   If multiple files, the HIT structure is filled from the file buffer with
!   the lowest wavenumber, and that buffer is reloaded with the next record
!   from the file.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE HFLCOM_DAT ! HITRAN line file general data
    USE HITCOM_DAT, ONLY: HITTYP ! HITRAN line data structure
!
! SUBROUTINES
    USE RECBIN_SUB ! Read record from HITRAN line data binary file
    USE RECHDB_SUB ! Read record from HITRAN database file
    USE RECPAR_SUB ! Read record from HITRAN line data .par file
!    
  IMPLICIT NONE
!
! ARGUMENTS      
    TYPE(HITTYP),  INTENT(OUT) :: HIT    ! HITRAN line parameters
    LOGICAL,       INTENT(OUT) :: EOF    ! Set TRUE if end-of-file reached
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IHFL           ! Counter for HITRAN data files
    INTEGER(I4) :: LUNHIT         ! LUN for HIT data file
    LOGICAL     :: USEIDM(MAXIDM) ! T=use molec# from file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  SELECT CASE ( NHFL ) 
    CASE ( 0 ) 
      EOF = .TRUE.
      RETURN
    CASE ( 1 )       ! Single HITRAN file, read directly into HIT
      LUNHIT = HFL(1)%LUN
      USEIDM = HFL(1)%IDM
      SELECT CASE ( HFL(1)%TYP )
        CASE ( 'BIN' ) 
          CALL RECBIN ( LUNHIT, HIT, USEIDM, HFL(1)%IFP, &
                        EOF, FAIL, ERRMSG )
        CASE ( 'HDB' ) 
          CALL RECHDB ( LUNHIT, HIT, USEIDM, WNOMAX, HFL(1)%HDB, &
                        EOF, FAIL, ERRMSG )
        CASE ( 'PAR' ) 
          CALL RECPAR ( LUNHIT, HIT, USEIDM, WNOMAX, &
                        EOF, FAIL, ERRMSG )
        CASE DEFAULT   ; STOP 'F-HITREC: Logical error'
      END SELECT
      IF ( FAIL ) RETURN
    CASE DEFAULT  ! Multiple files, read via buffers
      IHFL = MINLOC ( HFL(:)%HIT%WNO, 1 ) 
      EOF = HFL(IHFL)%HIT%WNO .GT. WNOMAX
      IF ( EOF ) RETURN
      HIT = HFL(IHFL)%HIT    ! Load from buffer
      LUNHIT = HFL(IHFL)%LUN
      USEIDM = HFL(IHFL)%IDM
      SELECT CASE ( HFL(IHFL)%TYP )  ! Load next into buffer
        CASE ( 'BIN' ) 
          CALL RECBIN ( LUNHIT, HFL(IHFL)%HIT, USEIDM, HFL(IHFL)%IFP, &
                        EOF, FAIL, ERRMSG )
        CASE ( 'HDB' ) 
          CALL RECHDB ( LUNHIT, HFL(IHFL)%HIT, USEIDM, WNOMAX, HFL(IHFL)%HDB, &
                        EOF, FAIL, ERRMSG )
        CASE ( 'PAR' ) 
          CALL RECPAR ( LUNHIT, HFL(IHFL)%HIT, USEIDM, WNOMAX, &
                        EOF, FAIL, ERRMSG )
        CASE DEFAULT   ; STOP 'F-HITREC: Logical error'
      END SELECT
      IF ( FAIL ) RETURN
      IF ( EOF ) HFL(IHFL)%HIT%WNO = 2 * WNOMAX
  END SELECT
!
END SUBROUTINE HITREC
END MODULE HITREC_SUB
