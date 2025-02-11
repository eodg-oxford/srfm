MODULE CIAINT_FNC
CONTAINS
FUNCTION CIAINT ( ID1, ID2, TEM, NWNO, WNOLST )
!
! VERSION
!   24JUL24 AD Checked.
!   01MAY17 AD F90 conversion of intcia.for. Checked.
!
! DESCRIPTION
!   Interpolate CIA spectrum
!   Called by SPCCIA.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CIACOM_DAT ! Collision-induced absorption data
!
! SUBROUTINES
    USE INTERP_GEN ! Interpolate array
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: ID1       ! RFM/HITRAN index of 1st collis molec.
    INTEGER(I4), INTENT(IN) :: ID2       ! RFM/HITRAN index of 2nd collis molec.
    REAL(R4),    INTENT(IN) :: TEM       ! Temperature [K]
    INTEGER(I4), INTENT(IN) :: NWNO      ! No.spectral points
    REAL(R8),    INTENT(IN) :: WNOLST(:) ! Wavenumber axis [cm-1]
!
! FUNCTION TYPE
    REAL(R4) :: CIAINT(NWNO)   ! Fnc returns array size NWNO
!
! LOCAL VARIABLES
    LOGICAL     :: FIRST  ! T=start new search for tabulations
    INTEGER(I4) :: ICIA   ! Counter for CIA tabulations
    INTEGER(I4) :: ICL    ! Index of lower temperature CIA tabulation
    INTEGER(I4) :: ICU    ! Index of upper temperature CIA tabulation
    REAL(R4)    :: A      ! Interpolation fraction
    REAL(R4)    :: TEMLOW ! Lower tabulated temperature [K]
    REAL(R4)    :: TEMUPP ! Upper tabulated temperature [K]
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  FIRST = .TRUE.
  DO ICIA = 1, NCIA
    IF ( CIA(ICIA)%ID1 .NE. ID1 ) CYCLE
    IF ( CIA(ICIA)%ID2 .NE. ID2 ) CYCLE
    IF ( CIA(ICIA)%WNU .LT. WNOLST(1) ) CYCLE
    IF ( CIA(ICIA)%WNL .GT. WNOLST(NWNO) ) CYCLE
! Found pair of molecules and wavenumber range, so check temperatures
    IF ( FIRST ) THEN      ! If first, use for both lower,upper TEM
      IF ( CIA(ICIA)%TEM .GE. TEM ) THEN
        TEMLOW = -1.0
        TEMUPP = CIA(ICIA)%TEM
      ELSE
        TEMLOW = CIA(ICIA)%TEM
        TEMUPP = 1.0E6
      END IF
      ICL = ICIA
      ICU = ICIA 
      FIRST = .FALSE.
    ELSE IF ( ( CIA(ICIA)%TEM .GT. TEMLOW .AND. CIA(ICIA)%TEM .LE. TEM ) .OR. &
              ( CIA(ICIA)%TEM .LT. TEMLOW .AND. TEMLOW .GT. TEM )        ) THEN
      TEMLOW = CIA(ICIA)%TEM
      ICL = ICIA
    ELSE IF ( ( CIA(ICIA)%TEM .LT. TEMUPP .AND. CIA(ICIA)%TEM .GE. TEM ) .OR. &
              ( CIA(ICIA)%TEM .GT. TEMUPP .AND. TEMUPP .LT. TEM )        ) THEN
      TEMUPP = CIA(ICIA)%TEM
      ICU = ICIA
    END IF
  END DO
!
  IF ( FIRST ) THEN                ! No CIA data within spectral range
    CIAINT = 0.0
  ELSE IF ( ICU .EQ. ICL ) THEN    ! Temperature beyond tabulated values
    CIAINT = INTERP ( CIA(ICL)%WNO, WNOLST, CIA(ICL)%ABS )
  ELSE                             ! Interpolate between tabulated values
    A = ( TEM - TEMLOW ) / ( TEMUPP - TEMLOW ) 
    CIAINT = (1.0 - A ) * INTERP ( CIA(ICL)%WNO, WNOLST, CIA(ICL)%ABS ) + &
                      A * INTERP ( CIA(ICU)%WNO, WNOLST, CIA(ICU)%ABS ) 
  END IF
!
END FUNCTION CIAINT
END MODULE CIAINT_FNC
           
