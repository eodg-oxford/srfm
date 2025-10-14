MODULE CIADEF_SUB
CONTAINS
SUBROUTINE CIADEF ( NAMDEF, FAIL, ERRMSG )
!
! VERSION
!   28JUL24 AD Checked.
!   24JUN19 AD Rewrite for through all molec.combinations, including Air.Checked
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Use default .cia filename to find any missing files
!   Called by DRVCIA if filename template found in *CIA section of driver table.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CIACOM_DAT ! Collision-induced absorption data
    USE GASCOM_DAT ! Molecule and isotope data
    USE LENREC_DAT ! Max length of input text record
    USE IDXCON_DAT, ONLY: IDXAIR ! RFM 'molecular' index for air
!
! SUBROUTINES
    USE CIAFIL_SUB ! Read Collision-Induced Abs. data file
    USE LEXIST_FNC ! Check if file exists
    USE UPCASE_FNC ! Convert text string to upper case
    USE USEMOL_FNC ! T = molecule is required
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: NAMDEF ! Default name of CIAfile
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4)       :: IGAS    ! Counter for first molecule in pair
    INTEGER(I4)       :: IPT     ! Location of '*' in NAMDEF
    INTEGER(I4)       :: JGAS    ! Counter for second molecule in pair
    CHARACTER(LENGAS) :: MOL1    ! Name of first molecule in pair
    CHARACTER(LENGAS) :: MOL2    ! Name of second molecule in pair
    CHARACTER(LENREC) :: NAMCIA  ! Name of CIA file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL WRTLOG ( 'I-CIADEF: Checking default CIA files' )
!
  IPT = INDEX ( NAMDEF, '*' )
  IF ( IPT .EQ. 0 ) STOP 'F-CIADEF: Logical error'
!
  DO IGAS = 1, NGAS
    DO JGAS = 1, NGAS
      IF ( NCIA .GT. 0 ) THEN 
        IF ( ANY ( CIA%IG1 .EQ. IGAS .AND. CIA%IG2 .EQ. JGAS ) .OR. &
             ANY ( CIA%IG2 .EQ. IGAS .AND. CIA%IG1 .EQ. JGAS ) ) CYCLE
      END IF
      MOL1 = UPCASE ( GAS(IGAS)%COD ) 
      MOL2 = UPCASE ( GAS(JGAS)%COD ) 
!
! Special cases for HITRAN CIA filenames
      IF ( MOL1 .EQ. 'HE'  ) MOL1 = 'He'
      IF ( MOL2 .EQ. 'HE'  ) MOL2 = 'He'
      IF ( MOL2 .EQ. 'AIR' ) MOL2 = 'Air'
!
      NAMCIA = NAMDEF(1:IPT-1) // TRIM(MOL1) // '-' // TRIM(MOL2) &
               // NAMDEF(IPT+1:)
      IF ( LEXIST ( NAMCIA ) ) THEN
        CALL CIAFIL ( NAMCIA, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END IF
    END DO
  END DO
!
! If 'air' not already specified in GASCOM, also try this as second absorber
  IF ( .NOT. USEMOL(IDXAIR) ) THEN 
    DO IGAS = 1, NGAS
      MOL1 = UPCASE ( GAS(IGAS)%COD ) 
      IF ( MOL1 .EQ. 'HE' ) MOL1 = 'He'
      NAMCIA = NAMDEF(1:IPT-1) // TRIM(MOL1) // '-Air' // NAMDEF(IPT+1:)
      IF ( LEXIST ( NAMCIA ) ) THEN
        CALL CIAFIL ( NAMCIA, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END IF
    END DO
  END IF      
!
END SUBROUTINE CIADEF
END MODULE CIADEF_SUB
