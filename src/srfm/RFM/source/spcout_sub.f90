MODULE SPCOUT_SUB
CONTAINS
SUBROUTINE SPCOUT ( ISPC, FAIL, ERRMSG )
!
! VERSION
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
    USE LEVCOM_DAT ! Intermediate output levels
    USE NAMCOM_DAT ! RFM output filenames
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
    REAL(R8), TARGET, ALLOCATABLE :: SPCFUL(:) ! Derived spectral outputs
    REAL(R8), POINTER             :: SPCPTR(:) ! Pointer for output spectrum
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NULLIFY ( SPCPTR ) 
  ALLOCATE ( SPCFUL(NFUL) )
!
  IF ( NJAC .GT. 0 ) THEN 
    NSEC = NJAC
  ELSE IF ( NLEV .GT. 0 ) THEN
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
        JTAN = ITAN
      ELSE                           ! secondary spectrum
        IF ( NJAC .GT. 0 ) THEN
          IJAC = ISEC
          JTAN = ITNJAC(ITAN,IJAC)
          IF ( JAC(IJAC)%COD .EQ. 'los' ) JTAN = -1 
        ELSE
          ILEV = ISEC
          JTAN = ITNLEV(ITAN,ILEV)
        END IF
        IF ( JTAN .EQ. 0 .AND. NOZERO ) CYCLE  ! No sec. output for this tan
        IF ( JTPFLG .AND. JTAN .NE. -1 ) THEN
          IF ( JAC(IJAC)%ITN .NE. ITAN ) CYCLE ! Only tan.pt Jacobians
        END IF
      END IF
!
      IF ( ABSFLG ) THEN
        IF ( JTAN .EQ. 0 ) THEN
          SPCFUL = 0.0D0
        ELSE IF ( JTAN .EQ. -1 ) THEN
          SPCFUL = 1.0D0 - SPCLOS ( ITAN, 'TRA' ) 
        ELSE
          SPCFUL = 1.0D0 - TRAFUL(:,JTAN)
        END IF
        SPCPTR => SPCFUL
        CALL SPCWRT ( ABSNAM, 'ABS', ISPC, ITAN, IJAC, ILEV, NFUL, IRRFUL, &
                      WNOFUL, SPCPTR, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END IF
!
      IF ( BBTFLG ) THEN
        IF ( IJAC .GT. 0 ) THEN
          IF ( JTAN .EQ. 0 ) THEN
            SPCFUL = 0.0D0
          ELSE IF ( JTAN .EQ. -1 ) THEN
            SPCFUL = BRIGHT ( SPCLOS(ITAN,'RAD') + RADFUL(:,ITAN), WNOFUL ) - &
                     BRIGHT ( RADFUL(:,ITAN), WNOFUL )  
          ELSE
            SPCFUL = BRIGHT ( RADFUL(:,JTAN) + RADFUL(:,ITAN), WNOFUL ) - &
                     BRIGHT ( RADFUL(:,ITAN), WNOFUL )
          END IF
        ELSE
          SPCFUL = BRIGHT ( RADFUL(:,JTAN), WNOFUL ) 
        END IF
        SPCPTR => SPCFUL
        CALL SPCWRT ( BBTNAM, 'BBT', ISPC, ITAN, IJAC, ILEV, NFUL, IRRFUL, &
                      WNOFUL, SPCPTR, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END IF
!
      IF ( COOFLG ) THEN
        SPCPTR => COOFUL(:,JTAN)
        CALL SPCWRT ( COONAM, 'COO', ISPC, ITAN, IJAC, ILEV, NFUL, IRRFUL, &
                      WNOFUL, SPCPTR, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END IF
!
      IF ( OPTFLG ) THEN
        IF ( JTAN .EQ. 0 ) THEN
          SPCFUL = 0.0D0
          SPCPTR => SPCFUL
        ELSE IF ( JTAN .EQ. -1 ) THEN
          SPCFUL = SPCLOS ( ITAN, 'OPT' )
          SPCPTR => SPCFUL
        ELSE
          SPCPTR => OPTFUL(:,JTAN)
        END IF
        CALL SPCWRT ( OPTNAM, 'OPT', ISPC, ITAN, IJAC, ILEV, NFUL, IRRFUL, &
                      WNOFUL, SPCPTR, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END IF
!
      IF ( RADFLG ) THEN
        IF ( JTAN .EQ. 0 ) THEN
          SPCFUL = 0.0
          SPCPTR => SPCFUL
        ELSE IF ( JTAN .EQ. -1 ) THEN
          SPCFUL = SPCLOS ( ITAN, 'RAD' )
          SPCPTR => SPCFUL
        ELSE IF ( FLXFLG .AND. .NOT. VRTFLG ) THEN
          SPCFUL = RADFUL(:,JTAN) * 1.0E-5   ! Rad.flux: convert nW/cm2 to W/m2 
          SPCPTR => SPCFUL
        ELSE
          SPCPTR => RADFUL(:,JTAN)
        END IF
        CALL SPCWRT ( RADNAM, 'RAD', ISPC, ITAN, IJAC, ILEV, NFUL, IRRFUL, &
                      WNOFUL, SPCPTR, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END IF
!
      IF ( RJTFLG ) THEN
        IF ( JTAN .EQ. 0 ) THEN
          SPCFUL = 0.0
        ELSE IF ( JTAN .EQ. -1 ) THEN
          SPCFUL = C2 * SPCLOS ( ITAN, 'RAD' ) / C1 / WNOFUL**2
        ELSE
          SPCFUL = C2 * RADFUL(:,JTAN) / C1 / WNOFUL**2
        END IF
        SPCPTR => SPCFUL
        CALL SPCWRT ( RJTNAM, 'RJT', ISPC, ITAN, IJAC, ILEV, NFUL, IRRFUL, &
                      WNOFUL, SPCPTR, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END IF
!
      IF ( TRAFLG ) THEN
        IF ( JTAN .EQ. 0 ) THEN
          SPCFUL = 0
          SPCPTR => SPCFUL
        ELSE IF ( JTAN .EQ. -1 ) THEN
          SPCFUL = SPCLOS ( ITAN, 'TRA' ) 
          SPCPTR => SPCFUL
        ELSE
          SPCPTR => TRAFUL(:,JTAN)
        END IF
        CALL SPCWRT ( TRANAM, 'TRA', ISPC, ITAN, IJAC, ILEV, NFUL, IRRFUL, &
                      WNOFUL, SPCPTR, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END IF
!
    END DO
  END DO
!
  IF ( WIDFLG ) CALL WRTSTT ( ISPC, FAIL, ERRMSG )
!
  NULLIFY ( SPCPTR ) 
!
END SUBROUTINE SPCOUT
END MODULE SPCOUT_SUB

