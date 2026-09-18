MODULE SPCOUT_SUB
CONTAINS
SUBROUTINE SPCOUT ( ISPC, FAIL, ERRMSG )
!
! VERSION
!   24JUL26 AD Checked.
!   01AUG25 AD Rewritten to redefine JTAN=0. Remove SPCPTR.
!   30MAR24 AD Checked.
!   04FEB19 AD Add SPCLOS to calculate LOS Jacobians
!   31JAN19 AD Fix Bug#14 - suppress off-diagonal outputs with JTP flag
!   08NOV17 AD F90 original. Checked.
!
! DESCRIPTION
!   Write spectral output data
!   Called by RFMSPC
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FLGCOM_DAT ! Option flags
    USE FULCOM_DAT ! Full grid data
    USE JACCOM_DAT ! Jacobian data
    USE NAMCOM_DAT ! RFM output filenames
    USE LEVCOM_DAT ! Intermediate output levels
    USE PHYCON_DAT, ONLY: C1,C2  ! Radiation constants
    USE TANCOM_DAT, ONLY: NTAN   ! No. of tangent paths for output
!
! SUBROUTINES
    USE BRIGHT_FNC ! Brightness Temperature calculation
    USE SPCLOS_FNC ! Calculate LOS Jacobian spectrum
    USE SPCWRT_SUB ! Write spectral data file
    USE WRTSTT_SUB ! Write widemesh statistics
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN)  :: ISPC   ! Spectral range number
    LOGICAL,       INTENT(OUT) :: FAIL   ! Set TRUE if a fatal error occurs
    CHARACTER(80), INTENT(OUT) :: ERRMSG ! Error message written if FAIL is TRUE
!
! LOCAL CONSTANTS
    LOGICAL, PARAMETER :: NOZERO = .FALSE. ! T=don't output zero Jacobian spectra
!
! LOCAL VARIABLES
    INTEGER(I4) :: IJAC ! Index of Jacobian
    INTEGER(I4) :: ILEV ! Index of intermediate output level
    INTEGER(I4) :: ISEC ! Counter for secondary output spectra
    INTEGER(I4) :: ITAN ! Counter for output tangent heights
    INTEGER(I4) :: JTAN ! Counter for tangent heights incl. Jacobian spectra
    INTEGER(I4) :: NSEC ! No. secondary output spectra per nominal spectrum
    REAL(R8)    :: SPCFUL(NFUL) ! Derived spectral outputs
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( NJAC .GT. 0 ) THEN 
    NSEC = NJAC
  ELSE IF ( MTXFLG ) THEN
    NSEC = NTAN 
  ELSE IF ( LEVFLG ) THEN
    NSEC = NLEV
  ELSE
    NSEC = 0
  END IF
!
  DO ITAN = 1, NTAN    
    DO ISEC = 0, NSEC     ! 0 = nominal spectra, 1:NSEC secondary spectra
      IF ( ISEC .EQ. 0 ) THEN        ! nominal spectrum
        IJAC = 0
        ILEV = 0
        JTAN = 0
      ELSE IF ( LEVFLG ) THEN                          ! secondary spectrum
        ILEV = ISEC
        JTAN = ITNLEV(ITAN,ILEV) 
      ELSE IF ( MTXFLG ) THEN
        JTAN = ITAN + ISEC*NTAN
      ELSE
        IJAC = ISEC
        JTAN = ITNJAC(ITAN,IJAC)
        IF ( JAC(IJAC)%COD .EQ. 'los' ) JTAN = -1 
      END IF
! No sec. output for this tan
      IF ( ISEC .GT. 0 .AND. JTAN .EQ. 0 .AND. NOZERO ) CYCLE  
      IF ( JTPFLG .AND. IJAC .GT. 0 .AND. JTAN .NE. -1 ) THEN
        IF ( JAC(IJAC)%ITN .NE. ITAN ) CYCLE ! Only tan.pt Jacobians
      END IF
!
      IF ( ABSFLG ) THEN
        IF ( ISEC .EQ. 0 ) THEN
          SPCFUL = 1.0D0 - TRAFUL(:,ITAN)
        ELSE IF ( JTAN .EQ. 0 ) THEN
          SPCFUL = 0.0D0
        ELSE IF ( JTAN .EQ. -1 ) THEN
          SPCFUL = 1.0D0 - SPCLOS ( ITAN, 'TRA' ) 
        ELSE
          SPCFUL = 1.0D0 - TRAFUL(:,JTAN)
        END IF
        CALL SPCWRT ( ABSNAM, 'ABS', NFUL, IRRFUL, WNOFUL, SPCFUL, &
                      FAIL, ERRMSG, &
                      IJAC=IJAC, ILEV=ILEV, ISPC=ISPC, ITAN=ITAN, JTAN=JTAN )
        IF ( FAIL ) RETURN
      END IF
!
      IF ( BBTFLG ) THEN
        IF ( ISEC .EQ. 0 ) THEN
          SPCFUL = BRIGHT ( RADFUL(:,ITAN), WNOFUL ) 
        ELSE IF ( JTAN .EQ. 0 ) THEN 
          SPCFUL = 0.0D0
        ELSE IF ( ILEV .GT. 0 ) THEN
          SPCFUL = BRIGHT ( RADFUL(:,JTAN), WNOFUL ) 
        ELSE IF ( JTAN .EQ. -1 ) THEN
          SPCFUL = BRIGHT ( SPCLOS(ITAN,'RAD') + RADFUL(:,ITAN), WNOFUL ) - &
                   BRIGHT ( RADFUL(:,ITAN), WNOFUL )  
        ELSE          ! other Jacobians
          SPCFUL = BRIGHT ( RADFUL(:,JTAN) + RADFUL(:,ITAN), WNOFUL ) - &
                   BRIGHT ( RADFUL(:,ITAN), WNOFUL )
        END IF
        CALL SPCWRT ( BBTNAM, 'BBT', NFUL, IRRFUL, WNOFUL, SPCFUL, &
                      FAIL, ERRMSG, &
                      IJAC=IJAC, ILEV=ILEV, ISPC=ISPC, ITAN=ITAN, JTAN=JTAN )
        IF ( FAIL ) RETURN
      END IF
!
      IF ( COOFLG ) THEN
        IF ( ISEC .EQ. 0 ) THEN
          SPCFUL = COOFUL(:,ITAN)
        ELSE
          SPCFUL = COOFUL(:,JTAN)
        END IF
        CALL SPCWRT ( COONAM, 'COO', NFUL, IRRFUL, WNOFUL, SPCFUL, &
                      FAIL, ERRMSG, &
                      IJAC=IJAC, ILEV=ILEV, ISPC=ISPC, ITAN=ITAN, JTAN=JTAN )
        IF ( FAIL ) RETURN
      END IF
!
      IF ( OPTFLG ) THEN
        IF ( ISEC .EQ. 0 ) THEN
          SPCFUL = OPTFUL(:,ITAN)
        ELSE IF ( JTAN .EQ. 0 ) THEN
          SPCFUL = 0.0D0
        ELSE IF ( JTAN .EQ. -1 ) THEN
          SPCFUL = SPCLOS ( ITAN, 'OPT' )
        ELSE
          SPCFUL = OPTFUL(:,JTAN)
        END IF
        CALL SPCWRT ( OPTNAM, 'OPT', NFUL, IRRFUL, WNOFUL, SPCFUL, &
                      FAIL, ERRMSG, &
                      IJAC=IJAC, ILEV=ILEV, ISPC=ISPC, ITAN=ITAN, JTAN=JTAN )
        IF ( FAIL ) RETURN
      END IF
!
      IF ( RADFLG ) THEN
        IF ( ISEC .EQ. 0 ) THEN
          SPCFUL = RADFUL(:,ITAN) 
        ELSE IF ( JTAN .EQ. 0 ) THEN
          SPCFUL = 0.0
        ELSE IF ( JTAN .EQ. -1 ) THEN
          SPCFUL = SPCLOS ( ITAN, 'RAD' )
        ELSE 
          SPCFUL = RADFUL(:,JTAN) 
        END IF
        IF ( FLXFLG .AND. .NOT. VRTFLG ) &
          SPCFUL = SPCFUL * 1.0E-5   ! Rad.flux: convert nW/cm2 to W/m2 
        CALL SPCWRT ( RADNAM, 'RAD', NFUL, IRRFUL, WNOFUL, SPCFUL, &
                      FAIL, ERRMSG, &
                      IJAC=IJAC, ILEV=ILEV, ISPC=ISPC, ITAN=ITAN, JTAN=JTAN )
        IF ( FAIL ) RETURN
      END IF
!
      IF ( RJTFLG ) THEN
        IF ( ISEC .EQ. 0 ) THEN
          SPCFUL = RADFUL(:,ITAN)
        ELSE IF ( JTAN .EQ. 0 ) THEN
          SPCFUL = 0.0
        ELSE IF ( JTAN .EQ. -1 ) THEN
          SPCFUL = SPCLOS ( ITAN, 'RAD' ) 
        ELSE
          SPCFUL = RADFUL(:,JTAN)
        END IF
        SPCFUL = C2 * SPCFUL / C1 / WNOFUL**2
        CALL SPCWRT ( RJTNAM, 'RJT', NFUL, IRRFUL, WNOFUL, SPCFUL, &
                      FAIL, ERRMSG, &
                      IJAC=IJAC, ILEV=ILEV, ISPC=ISPC, ITAN=ITAN, JTAN=JTAN )
        IF ( FAIL ) RETURN
      END IF
!
      IF ( TRAFLG ) THEN
        IF ( ISEC .EQ. 0 ) THEN
          SPCFUL = TRAFUL(:,ITAN)
        ELSE IF ( JTAN .EQ. 0 ) THEN
          SPCFUL = 0
        ELSE IF ( JTAN .EQ. -1 ) THEN
          SPCFUL = SPCLOS ( ITAN, 'TRA' ) 
        ELSE
          SPCFUL = TRAFUL(:,JTAN)
        END IF
        CALL SPCWRT ( TRANAM, 'TRA', NFUL, IRRFUL, WNOFUL, SPCFUL, &
                      FAIL, ERRMSG, &
                      IJAC=IJAC, ILEV=ILEV, ISPC=ISPC, ITAN=ITAN, JTAN=JTAN )
        IF ( FAIL ) RETURN
      END IF
!
    END DO
  END DO
!
  IF ( WIDFLG ) CALL WRTSTT ( ISPC, FAIL, ERRMSG )
!
END SUBROUTINE SPCOUT
END MODULE SPCOUT_SUB

