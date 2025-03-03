MODULE JACTGT_SUB
CONTAINS
SUBROUTINE JACTGT ( TARGET, JDX, FAIL, ERRMSG )
!
! VERSION
!   14JUL24 AD Checked.
!   01MAY17 AD F90 conversion of part of jacchk.for. Checked.
!
! DESCRIPTION
!   Check valid Jacobian target 
!   Called by DRVJAC for first field in each record of *JAC sec. of drv.table.
!   (corresponding to each target species for Jacobian).
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE JACCOM_DAT ! Jacobian data
!
! SUBROUTINES
    USE CHKGAS_SUB ! Check for valid molecule name, isotope, Vib.Level
    USE IDXGAS_FNC ! Index in GAS arrays of molecule,isotope
    USE JACISO_SUB ! Check valid parameters for isotopic Jacobian
    USE USEMOL_FNC ! T = molecule is required
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: TARGET ! Target for Jacobian
    INTEGER(I4),   INTENT(OUT) :: JDX    ! Jacobian target type identifier
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error is detected
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IDXISO ! Isotope#1 (1=most abundant, etc, 0=all isotopes)
    INTEGER(I4) :: IDXMOL ! HITRAN/RFM index of molecule
    INTEGER(I4) :: IDXVIB ! Global Quantum index of Vibrational level
!
! EXECUTABLE CODE ------------------------------------------------------------
!
  FAIL = .FALSE.
!
  SELECT CASE ( TARGET ) 
  CASE ( 'tem' ) 
    JDX = JDXTEM
  CASE ( 'pre' ) 
    JDX = JDXPRE
  CASE ( 'sfctem' ) 
    IF ( .NOT. SFCFLG ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-JACTGT: SFCTEM Jacobian requires SFC flag to be enabled'
      RETURN
    END IF
    JDX = JDXSFT
  CASE ( 'sfcems' ) 
    IF ( .NOT. SFCFLG ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-JACTGT: SFCEMS Jacobian requires SFC flag to be enabled'
      RETURN
    END IF
    JDX = JDXSFE
  CASE DEFAULT
    CALL CHKGAS ( TARGET, IDXMOL, IDXISO, IDXVIB, FAIL, ERRMSG )
    IF ( FAIL ) RETURN
    IF ( IDXMOL .EQ. 0 ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-JACTGT: unrecognised target for Jacobian calc,=' // TARGET
      RETURN
    ELSE IF ( .NOT. USEMOL(IDXMOL) ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-JACTGT: No ' // TARGET // ' profile loaded for Jacobian calc.'
      RETURN
    ELSE IF ( IDXISO .EQ. 0 ) THEN         ! Molecular Jacobian
      JDX = IDXGAS ( IDXMOL ) 
    ELSE IF ( IDXVIB .EQ. 0 ) THEN         ! Isotopic Jacobian
      CALL JACISO ( IDXMOL, IDXISO, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      JDX = IDXGAS ( IDXMOL, IDXISO ) 
    ELSE                                   ! Vib Tem Jacobian
      IF ( .NOT. NTEFLG ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-JACTGT: Vib Temp Jacobian requires NTE flag to be enabled'
        RETURN
      END IF
      JDX = IDXVIB
    END IF
  END SELECT
!
  IF ( NJAC .GT. 0 ) THEN
    IF ( ANY ( JAC%JDX .EQ. JDX ) ) THEN
      FAIL = .TRUE.
      ERRMSG = 'F-JACTGT: Repeated target for Jacobian calc.=' // TARGET
      RETURN
    END IF
  END IF
!
END SUBROUTINE JACTGT
END MODULE JACTGT_SUB
