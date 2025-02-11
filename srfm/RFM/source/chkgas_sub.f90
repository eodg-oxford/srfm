MODULE CHKGAS_SUB
CONTAINS
SUBROUTINE CHKGAS ( INPSTR, IDXMOL, IDXISO, IDXVIB, FAIL, ERRMSG, ADDNEW )
!
! VERSION
!   03SEP24 AD Checked.
!   01JUL19 AD Bug#22 Increase length of LSTR by +3. Checked.
!   05OCT17 AD F90 conversion.Checked.
!
! DESCRIPTION
!   Check for valid molecule name, isotope, Vib.Level
!   General purpose module.
!   Input string is assumed to be of form: molec(i)(v)
!   where brackets are optional, but (v) requires preceding (i).
!   i = isotope#, v = vib.level (for Vib Tem retrieval), 
!   If argument FAIL is absent, IDXMOL=0 used to indicate unidentified string
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! SUBROUTINES
    USE IDGSTR_SUB ! Decode string HITRAN/RFM index and isotope
    USE ISOLST_SUB ! List isotopomer weights of specific molecule
    USE LOCASE_FNC ! Convert text string to lower case
    USE MOLIDX_SUB ! Give molecule name for HITRAN/RFM index, or vice-versa
    USE REAQAL_SUB ! Decode Qualifier strings '(IQAL)' or '(IQAL:JQAL)'
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),  INTENT(IN)  :: INPSTR ! String to be tested 
    INTEGER(I4),   INTENT(OUT) :: IDXMOL ! HITRAN/RFM Index, 0 if unrecognised
    INTEGER(I4),   INTENT(OUT) :: IDXISO ! Isotope#n, 0 if not isotopomer
    INTEGER(I4), &
         OPTIONAL, INTENT(OUT) :: IDXVIB ! Encoded Vib.level, else 0
    LOGICAL, &
         OPTIONAL, INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error 
    CHARACTER(80), &
         OPTIONAL, INTENT(OUT) :: ERRMSG ! Error message if FAIL is TRUE
    LOGICAL, &
         OPTIONAL, INTENT(IN)  :: ADDNEW ! T=allow new molecs to be defined
!
! LOCAL VARIABLES
    LOGICAL      :: LFAIL ! T=Fatal error returned
    INTEGER(I4)  :: ILV   ! Lower state vibrational level index
    INTEGER(I4)  :: IPT   ! location in LSTR of first '(' character
    INTEGER(I4)  :: NISO  ! No. of isotopes for molecule
    CHARACTER(80)            :: LERMSG ! Local error message
    CHARACTER(LEN(INPSTR)+3) :: LSTR   ! lower case copy of INPSTR
    CHARACTER(LEN(INPSTR))   :: GASSTR ! part of LSTR containing molecule name
    CHARACTER(7)             :: MOLNAM ! Name of molecule
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( PRESENT ( FAIL ) ) FAIL = .FALSE.
!
  LSTR = LOCASE ( INPSTR )  ! convert to lower case
  IDXMOL = 0
  IDXISO = 0
  IF ( PRESENT ( IDXVIB ) ) IDXVIB = 0
! 
! Special cases expand as molecule + isotope (Length of LSTR increased here)
  IF ( INDEX ( LSTR, 'hdo' ) .EQ. 1 ) THEN
    IF ( LEN ( LSTR ) .EQ. 3 ) THEN  
      LSTR = 'h2o(4)' 
    ELSE
      LSTR = 'h2o(4)' // LSTR(4:)
    END IF
    CALL WRTLOG ( 'I-CHKGAS: Expanding HDO as H2O(4)' ) 
  ELSE IF ( INDEX ( LSTR, 'ch3d' ) .EQ. 1 ) THEN
    IF ( LEN ( LSTR ) .EQ. 4 ) THEN
      LSTR = 'ch4(3)' 
    ELSE
      LSTR = 'ch4(3)' // LSTR(5:)
    END IF
    CALL WRTLOG ( 'I-CHKGAS: Expanding CH3D as CH4(3)' ) 
  END IF
!
! Extract part of string up to any (optional) bracket
  GASSTR = LSTR
  IPT = INDEX ( LSTR, '(' ) 
  IF ( IPT .GT. 0 ) THEN 
    GASSTR = LSTR(1:IPT-1)
    LSTR = LSTR(IPT:)
  END IF
!
! See if GASSTR contains numerical code for molecule (and isotope)
  CALL IDGSTR ( GASSTR, IDXMOL, IDXISO ) 
  IF ( IDXMOL .GT. 0 ) THEN          ! interpreted molec index
    CALL MOLIDX ( IDXMOL, MOLNAM )   ! check if valid molec#
    IF ( MOLNAM .EQ. ' ' ) THEN      ! not valid so reset IDXMOL=0
      IDXMOL = 0
      IDXISO = 0
    END IF
  END IF
!  
! Set if GASSTR contains molecule name, in which case set IDXMOL to molec#
  IF ( IDXMOL .EQ. 0 ) THEN
    IF ( PRESENT ( ADDNEW ) ) THEN             ! allow new molec names 
      CALL MOLIDX ( IDXMOL, GASSTR, ADDNEW )
    ELSE
      CALL MOLIDX ( IDXMOL, GASSTR )
    END IF
  END IF
!
  IF ( IDXMOL .EQ. 0 ) THEN         ! No molecule identifiable
    LERMSG = 'F-CHKGAS: Unrecognised retrieval molecule: '//GASSTR
    GOTO 900
  END IF
!
! Check for isotope, first (remaining) '(' character
  IF ( IDXISO .EQ. 0 ) THEN        ! No isotope# yet 
    IF ( LSTR(1:1) .EQ. '(' ) THEN
      CALL REAQAL ( LSTR, IDXISO, LFAIL, LERMSG ) 
      IF ( LFAIL ) GOTO 900
    END IF
  END IF
!
! Normal exit with just molecule#
  IF ( IDXISO .EQ. 0 ) RETURN 
!
! Check valid isotope#
  CALL ISOLST ( IDXMOL, NISO )
  IF ( IDXISO .LT. 0 .OR. IDXISO .GT. NISO ) THEN
    WRITE ( LERMSG, * ) 'F-CHKGAS: molecule ' // GASSTR // &
                        ' has invalid Isotope#=', IDXISO
    GOTO 900
  END IF 
!
! Check for (optional) lower state quantum number
  IF ( LSTR(1:1) .EQ. ' ' ) RETURN
!
! Case where vib level present but not expected 
  IF ( .NOT. PRESENT ( IDXVIB ) ) THEN
    LERMSG = 'F-CHKGAS: Vib.level not expected for molec/iso'
    GOTO 900
  END IF
!
  CALL REAQAL ( LSTR, ILV, LFAIL, LERMSG ) 
  IF ( LFAIL ) GOTO 900
  IF ( ILV .LE. 1 ) THEN    ! 1=ground state, not valid for vib tem.
    WRITE ( LERMSG, * ) 'F-CHKGAS: retrieved Vib Temp for ' // GASSTR // &
                        ' has invalid Vib.level#=', ILV
    GOTO 900
  END IF
!
! Normal exit with IDXMOL, IDXISO and IDXVIB all set
  IDXVIB = IDXMOL * 1000000 + IDXISO * 1000 + ILV ! should match CHKVIB
  RETURN   
!
! Fatal error handling
900 CONTINUE
  IF ( PRESENT ( FAIL ) ) THEN
    FAIL = .TRUE.
    ERRMSG = LERMSG
  END IF
  IDXMOL = 0
!
END SUBROUTINE CHKGAS
END MODULE CHKGAS_SUB
