MODULE VAL1DI_GEN
!
! VERSION
!   29MAY24 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Interpolate value from 1D array
!   General purpose module
!   If EXTRAP is set TRUE, then outside points are extrapolated (if NXTAB>1)
!   otherwise the end points are duplicated.
!   If LOGINT is set TRUE, interpolation is linear in ln(YTAB). 
!   If XTAB is a single element, output is replicated YTAB value.
!
! VARIABLE KINDS
    USE KIND_DAT
!
INTERFACE VAL1DI
  MODULE PROCEDURE VAL1DI_RR, VAL1DI_DD, VAL1DI_DR
END INTERFACE

CONTAINS

REAL(R4) PURE FUNCTION VAL1DI_RR ( XTAB, XINT, YTAB, LOGINT, EXTRAP )
!
! SUBROUTINES
    USE IBRAKT_GEN ! Lower index of array interpolation
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: XTAB(:) ! List of tabulated coordinates
    REAL(R4), INTENT(IN) :: XINT    ! interpolation coordinate
    REAL(R4), INTENT(IN) :: YTAB(:) ! List of tabulated data values at XTAB
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: LOGINT  ! TRUE=interpolate linearly in Log(YTAB)
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: EXTRAP  ! TRUE=extrapolate beyond ends of XTAB
!
! LOCAL VARIABLES
    LOGICAL     :: LEXT   ! T=extrapolation, F=no extrapolation
    LOGICAL     :: LINT   ! T=Log interpolation, F=linear interpolation
    INTEGER(I4) :: IX     ! Index of lower coordinate in XTAB
    INTEGER(I4) :: NTAB   ! Size of XTAB,YTAB arrays
    REAL(R4)    :: DX     ! XINT as fraction of XTAB interval
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NTAB = SIZE ( XTAB )
!
! Special case of just a single element in tabulated values
  IF ( NTAB .EQ. 1 ) THEN
    VAL1DI_RR = YTAB(1)
    RETURN
  END IF
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
  IF ( LEXT ) THEN
    IX = IBRAKT ( XTAB, XINT, LIMIT=.TRUE. )  ! T=limit IX to 1:NTAB-1
  ELSE
    IX = IBRAKT ( XTAB, XINT, LIMIT=.FALSE. ) 
  END IF
!
  IF ( IX .EQ. 0 ) THEN
    IX = 1
    DX = 0.0
  ELSE IF ( IX .EQ. NTAB ) THEN
    IX = NTAB - 1
    DX = 1.0
  ELSE 
    DX = ( XINT - XTAB(IX) ) / ( XTAB(IX+1) - XTAB(IX) )
  END IF
!
  IF ( LINT ) THEN
    VAL1DI_RR = EXP ( (1.0-DX) * LOG ( MAX ( TINY(1.0), YTAB(IX)   ) ) + &
                        DX  * LOG ( MAX ( TINY(1.0), YTAB(IX+1) ) )    ) 
  ELSE
    VAL1DI_RR = (1.0-DX) * YTAB(IX) + DX * YTAB(IX+1) 
  END IF
!
END FUNCTION VAL1DI_RR

REAL(R8) PURE FUNCTION VAL1DI_DD ( XTAB, XINT, YTAB, LOGINT, EXTRAP )
!
! SUBROUTINES
    USE IBRAKT_GEN ! Lower index of array interpolation
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(IN) :: XTAB(:) ! List of tabulated coordinates
    REAL(R8), INTENT(IN) :: XINT    ! interpolation coordinate
    REAL(R8), INTENT(IN) :: YTAB(:) ! List of tabulated data values at XTAB
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: LOGINT  ! TRUE=interpolate linearly in Log(YTAB)
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: EXTRAP  ! TRUE=extrapolate beyond ends of XTAB
!
! LOCAL VARIABLES
    LOGICAL     :: LEXT   ! T=extrapolation, F=no extrapolation
    LOGICAL     :: LINT   ! T=Log interpolation, F=linear interpolation
    INTEGER(I4) :: IX     ! Index of lower coordinate in XTAB
    INTEGER(I4) :: NTAB   ! Size of XTAB,YTAB arrays
    REAL(R8)    :: DX     ! XINT as fraction of XTAB interval
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NTAB = SIZE ( XTAB )
!
! Special case of just a single element in tabulated values
  IF ( NTAB .EQ. 1 ) THEN
    VAL1DI_DD = YTAB(1)
    RETURN
  END IF
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
  IF ( LEXT ) THEN
    IX = IBRAKT ( XTAB, XINT, LIMIT=.TRUE. )  ! T=limit IX to 1:NTAB-1
  ELSE
    IX = IBRAKT ( XTAB, XINT, LIMIT=.FALSE. ) 
  END IF
!
  IF ( IX .EQ. 0 ) THEN
    IX = 1
    DX = 0.0D0
  ELSE IF ( IX .EQ. NTAB ) THEN
    IX = NTAB - 1
    DX = 1.0D0
  ELSE 
    DX = ( XINT - XTAB(IX) ) / ( XTAB(IX+1) - XTAB(IX) )
  END IF
!
  IF ( LINT ) THEN
    VAL1DI_DD = EXP ( (1.0D0-DX) * LOG ( MAX ( TINY(1.0D0), YTAB(IX)   ) ) + &
                        DX  * LOG ( MAX ( TINY(1.0D0), YTAB(IX+1) ) )    ) 
  ELSE
    VAL1DI_DD = (1.0D0-DX) * YTAB(IX) + DX * YTAB(IX+1) 
  END IF
!
END FUNCTION VAL1DI_DD

REAL(R4) PURE FUNCTION VAL1DI_DR ( XTAB, XINT, YTAB, LOGINT, EXTRAP )
!
! SUBROUTINES
    USE IBRAKT_GEN ! Lower index of array interpolation
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8), INTENT(IN) :: XTAB(:) ! List of tabulated coordinates
    REAL(R8), INTENT(IN) :: XINT    ! interpolation coordinate
    REAL(R4), INTENT(IN) :: YTAB(:) ! List of tabulated data values at XTAB
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: LOGINT  ! TRUE=interpolate linearly in Log(YTAB)
    LOGICAL, OPTIONAL, &
              INTENT(IN) :: EXTRAP  ! TRUE=extrapolate beyond ends of XTAB
!
! LOCAL VARIABLES
    LOGICAL     :: LEXT   ! T=extrapolation, F=no extrapolation
    LOGICAL     :: LINT   ! T=Log interpolation, F=linear interpolation
    INTEGER(I4) :: IX     ! Index of lower coordinate in XTAB
    INTEGER(I4) :: NTAB   ! Size of XTAB,YTAB arrays
    REAL(R4)    :: DX     ! XINT as fraction of XTAB interval
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NTAB = SIZE ( XTAB )
!
! Special case of just a single element in tabulated values
  IF ( NTAB .EQ. 1 ) THEN
    VAL1DI_DR = YTAB(1)
    RETURN
  END IF
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
  IF ( LEXT ) THEN
    IX = IBRAKT ( XTAB, XINT, LIMIT=.TRUE. )  ! T=limit IX to 1:NTAB-1
  ELSE
    IX = IBRAKT ( XTAB, XINT, LIMIT=.FALSE. ) 
  END IF
!
  IF ( IX .EQ. 0 ) THEN
    IX = 1
    DX = 0.0
  ELSE IF ( IX .EQ. NTAB ) THEN
    IX = NTAB - 1
    DX = 1.0
  ELSE 
    DX = SNGL ( ( XINT - XTAB(IX) ) / ( XTAB(IX+1) - XTAB(IX) ) ) 
  END IF
!
  IF ( LINT ) THEN
    VAL1DI_DR = EXP ( (1.0-DX) * LOG ( MAX ( TINY(1.0), YTAB(IX)   ) ) + &
                        DX  * LOG ( MAX ( TINY(1.0), YTAB(IX+1) ) )    ) 
  ELSE
    VAL1DI_DR = (1.0-DX) * YTAB(IX) + DX * YTAB(IX+1) 
  END IF
!
END FUNCTION VAL1DI_DR

END MODULE VAL1DI_GEN

