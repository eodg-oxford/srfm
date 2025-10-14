MODULE CHICO2_FNC
CONTAINS
FUNCTION CHICO2 ( TEM, PRE, PPA, WNOADJ, DWNO )
!
! VERSION
!   12FEB24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Chi-Factor for CO2 lines
!   Called by CHISHP for CO2 lines.
!   Note that this is a CPU-intensive routine.
!
! REFERENCES
!   C.Cousin,R.Le Doucen,C.Boulet and A.Henry, (1985)
!     Temperature dependence of the absorption in the region beyond the 4.3 
!     mic band head of CO2. 2: N2 and O2 broadening.
!     Appl.Opt. 24 (22) 3899-3907
!
!   Pure CO2 chi factor: Appl.Opt. 24 (6) 897-906 (1985)
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CHIDAT_DAT ! Chi-factor data
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: TEM     ! Temperature [K] 
    REAL(R4), INTENT(IN) :: PRE     ! Pressure [atm]
    REAL(R4), INTENT(IN) :: PPA     ! CO2 Partial Pressure [atm]
    REAL(R8), INTENT(IN) :: WNOADJ  ! Adjusted line wavenumber [cm-1]
    REAL(R8), INTENT(IN) :: DWNO(:) ! Array of Wavenumbers [cm-1]
!
! FUNCTION TYPE
    REAL(R4) :: CHICO2 ( SIZE(DWNO) )  ! Function returns array size of DWNO
!
! LOCAL VARIABLES
    INTEGER(I4) :: IFRQ        ! Index of tabulated data        
    INTEGER(I4) :: IWNO        ! Wavenumber array counter
    INTEGER(I4) :: JFRQ        ! Upper index for wavenumber interpolation
    INTEGER(I4) :: NWNO        ! No. Wno. points to be evaluted
    REAL(R4)    :: DFOR        ! Weights for foreign-broadening term
    REAL(R4)    :: DSLF        ! CO2 partial/total pressure [0-1]
    REAL(R4)    :: FFI,FFJ     ! Weights for lower,upper waveno interpolation
    REAL(R4)    :: RFRQ        ! Distance from line centre as a multiple of DFRQ
    REAL(R4)    :: TF1,TF2,TF3 ! Weights for lower,upper for.broad temperatures
    REAL(R4)    :: TS4, TS5    ! Weights for lower,upper self broad temperatures
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  NWNO = SIZE ( DWNO )
!
! Foreign broadening fraction
  DSLF = PPA / PRE
  DFOR = 1.0 - DSLF
!
! Temperature interpolation (foreign broadened)
  IF ( TEM .LT. TEMF1 ) THEN
    TF1 = DFOR
    TF2 = 0.0
    TF3 = 0.0
  ELSE IF ( TEM .LT. TEMF2 ) THEN
    TF1 = ( TEMF2 - TEM ) / ( TEMF2 - TEMF1 ) 
    TF2 = ( 1.0 - TF1 ) * DFOR
    TF1 = TF1 * DFOR
    TF3 = 0.0
  ELSE IF ( TEM .LT. TEMF3 ) THEN
    TF1 = 0.0
    TF2 = ( TEMF3 - TEM ) / ( TEMF3 - TEMF2 )
    TF3 = ( 1.0 - TF2 ) * DFOR
    TF2 = TF2 * DFOR
  ELSE 
    TF1 = 0.0
    TF2 = 0.0
    TF3 = DFOR
  END IF
!
! Temperature interpolation (self broadened)
  IF ( TEM .LT. TEMS4 ) THEN
    TS4 = DSLF
    TS5 = 0.0
  ELSE IF ( TEM .LT. TEMS5 ) THEN
    TS4 = ( TEMS5 - TEM ) / ( TEMS5 - TEMS4 )
    TS5 = ( 1.0 - TS4 ) * DSLF
    TS4 = TS4 * DSLF
  ELSE
    TS4 = 0.0
    TS5 = DSLF
  END IF
!
  DO IWNO = 1, NWNO 
!
! Frequency interpolation
    RFRQ = 1.0 + ABS ( SNGL ( DWNO(IWNO) - WNOADJ ) ) / DFRQ 
    IF ( RFRQ .LT. NFRQ ) THEN 
      IFRQ = INT ( RFRQ ) 
      FFJ = ( RFRQ - IFRQ ) / DFRQ
    ELSE
      IFRQ = NFRQ - 1
      FFJ = 1.0
    END IF
    JFRQ = 1 + IFRQ
    FFI = ( 1.0 - FFJ ) 
!     
    CHICO2(IWNO) = TF1 * ( FFI * CHI(IFRQ,1) + FFJ * CHI(JFRQ,1) ) &
                 + TF2 * ( FFI * CHI(IFRQ,2) + FFJ * CHI(JFRQ,2) ) &
                 + TF3 * ( FFI * CHI(IFRQ,3) + FFJ * CHI(JFRQ,3) ) &
                 + TS4 * ( FFI * CHI(IFRQ,4) + FFJ * CHI(JFRQ,4) ) &
                 + TS5 * ( FFI * CHI(IFRQ,5) + FFJ * CHI(JFRQ,5) )
  END DO
!
END FUNCTION CHICO2
END MODULE CHICO2_FNC
