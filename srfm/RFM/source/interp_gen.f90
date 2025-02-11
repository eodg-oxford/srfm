MODULE INTERP_GEN
!
! VERSION
!   05SEP24 AD Checked.
!   06SEP21 AD Add INTERP_RRS
!   01MAY17 AD Original. Checked.
!
! DESCRIPTION
!   Interpolate array
!   General purpose module.
!   Returns an array or scalar of interpolated values.
!   If EXTRAP is set FALSE, then where XTAB<XINT or XTAB>XINT the end data 
!   values of the YTAB array are duplicated.
!   If EXTRAP is set TRUE, then outside points are extrapolated (if NXTAB>1)
!   If LOGINT is set TRUE, interpolation is linear in ln(YTAB). 
!   If XTAB is a single element, output is replicated YTAB value.
!
! VARIABLE KINDS
    USE KIND_DAT
!
INTERFACE INTERP
  MODULE PROCEDURE INTERP_RR, INTERP_RRS, INTERP_DD, INTERP_DR
END INTERFACE

CONTAINS

PURE FUNCTION INTERP_RR ( XTAB, XINT, YTAB, LOGINT, EXTRAP )
!
! SUBROUTINES 
    USE VAL1DI_GEN ! Interpolate value from 1D array
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: XTAB(:) ! List of tabulated coordinates
    REAL(R4), INTENT(IN) :: XINT(:) ! List of interpolation coordinates
    REAL(R4), INTENT(IN) :: YTAB(:) ! List of tabulated data values at XTAB
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: LOGINT  ! TRUE=interpolate linearly in Log(YTAB)
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: EXTRAP  ! TRUE=extrapolate beyond ends of XTAB
!
! FUNCTION TYPE
    REAL(R4) :: INTERP_RR ( SIZE(XINT) ) ! Function returns array size of XINT 
!
! LOCAL VARIABLES
    LOGICAL     :: LINT ! T=Log interpolation, F=linear interpolation
    LOGICAL     :: LEXT ! T=extrapolation, F=no extrapolation
    INTEGER(I4) :: I    ! Counter for interpolated points
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( PRESENT ( LOGINT ) ) THEN
    LINT = LOGINT
  ELSE
    LINT = .FALSE.
  END IF
!
  IF ( PRESENT ( EXTRAP ) ) THEN
    LEXT = EXTRAP
  ELSE
    LEXT = .FALSE.
  END IF
!
  DO I = 1, SIZE ( XINT ) 
    INTERP_RR(I) = VAL1DI ( XTAB, XINT(I), YTAB, LINT, LEXT ) 
  END DO
!
END FUNCTION INTERP_RR

PURE FUNCTION INTERP_RRS ( XTAB, XINT, YTAB, LOGINT, EXTRAP )
!
! SUBROUTINES 
    USE VAL1DI_GEN ! Interpolate value from 1D array
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: XTAB(:) ! List of tabulated coordinates
    REAL(R4), INTENT(IN) :: XINT    ! Interpolation coordinate
    REAL(R4), INTENT(IN) :: YTAB(:) ! List of tabulated data values at XTAB
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: LOGINT  ! TRUE=interpolate linearly in Log(YTAB)
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: EXTRAP  ! TRUE=extrapolate beyond ends of XTAB
!
! FUNCTION TYPE
    REAL(R4) :: INTERP_RRS ! Function returns scalar
!
! LOCAL VARIABLES
    LOGICAL     :: LINT ! T=Log interpolation, F=linear interpolation
    LOGICAL     :: LEXT ! T=extrapolation, F=no extrapolation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( PRESENT ( LOGINT ) ) THEN
    LINT = LOGINT
  ELSE
    LINT = .FALSE.
  END IF
!
  IF ( PRESENT ( EXTRAP ) ) THEN
    LEXT = EXTRAP
  ELSE
    LEXT = .FALSE.
  END IF
!
  INTERP_RRS = VAL1DI ( XTAB, XINT, YTAB, LINT, LEXT ) 
!
END FUNCTION INTERP_RRS

PURE FUNCTION INTERP_DD ( XTAB, XINT, YTAB, LOGINT, EXTRAP )
!
! SUBROUTINES 
    USE VAL1DI_GEN ! Interpolate value from 1D array
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(IN) :: XTAB(:) ! List of tabulated coordinates
    REAL(R8), INTENT(IN) :: XINT(:) ! List of interpolation coordinates
    REAL(R8), INTENT(IN) :: YTAB(:) ! List of tabulated data values at XTAB
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: LOGINT  ! TRUE=interpolate linearly in Log(YTAB)
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: EXTRAP  ! TRUE=extrapolate beyond ends of XTAB
!
! FUNCTION TYPE
    REAL(R8) :: INTERP_DD ( SIZE(XINT) ) ! Function returns array size of XINT 
!
! LOCAL VARIABLES
    LOGICAL     :: LINT   ! T=Log interpolation, F=linear interpolation
    LOGICAL     :: LEXT   ! T=extrapolation, F=no extrapolation
    INTEGER(I4) :: I    ! Counter for interpolated points
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( PRESENT ( LOGINT ) ) THEN
    LINT = LOGINT
  ELSE
    LINT = .FALSE.
  END IF
!
  IF ( PRESENT ( EXTRAP ) ) THEN
    LEXT = EXTRAP
  ELSE
    LEXT = .FALSE.
  END IF
!
  DO I = 1, SIZE ( XINT ) 
    INTERP_DD(I) = VAL1DI ( XTAB, XINT(I), YTAB, LINT, LEXT ) 
  END DO
!
END FUNCTION INTERP_DD

PURE FUNCTION INTERP_DR ( XTAB, XINT, YTAB, LOGINT, EXTRAP )
!
! SUBROUTINES 
    USE VAL1DI_GEN ! Interpolate value from 1D array
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(IN) :: XTAB(:) ! List of tabulated coordinates
    REAL(R8), INTENT(IN) :: XINT(:) ! List of interpolation coordinates
    REAL(R4), INTENT(IN) :: YTAB(:) ! List of tabulated data values at XTAB
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: LOGINT  ! TRUE=interpolate linearly in Log(YTAB)
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: EXTRAP  ! TRUE=extrapolate beyond ends of XTAB
!
! FUNCTION TYPE
    REAL(R4) :: INTERP_DR ( SIZE(XINT) ) ! Function returns array size of XINT 
!
! LOCAL VARIABLES
    LOGICAL     :: LINT ! T=Log interpolation, F=linear interpolation
    LOGICAL     :: LEXT ! T=extrapolation, F=no extrapolation
    INTEGER(I4) :: I    ! Counter for interpolated points
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( PRESENT ( LOGINT ) ) THEN
    LINT = LOGINT
  ELSE
    LINT = .FALSE.
  END IF
!
  IF ( PRESENT ( EXTRAP ) ) THEN
    LEXT = EXTRAP
  ELSE
    LEXT = .FALSE.
  END IF
!
  DO I = 1, SIZE ( XINT ) 
    INTERP_DR(I) = VAL1DI ( XTAB, XINT(I), YTAB, LINT, LEXT ) 
  END DO
!
END FUNCTION INTERP_DR

END MODULE INTERP_GEN

