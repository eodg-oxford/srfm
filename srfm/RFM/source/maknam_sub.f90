MODULE MAKNAM_SUB
CONTAINS
SUBROUTINE MAKNAM ( INPNAM, OUTNAM, IGAS, IJAC, ILEV, IPSI, ISPC, ITAN )
!
! VERSION
!   09SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Construct filename for RFM output files.
!   General purpose module.
!   CMINUS can be changed to some other character if operating system does 
!   not allow '-' sign in filenames.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE JACCOM_DAT ! Jacobian data
    USE LEVCOM_DAT ! Intermediate output levels
    USE RUN_ID_DAT ! RFM Run ID string
    USE SPCCOM_DAT ! Spectral range data
    USE TANCOM_DAT ! Tangent path data
    USE GRACOM_DAT, ONLY:PSIGRA ! Psi values [deg] of each 2D atm location
!
! SUBROUTINES
    USE HGTSTR_FNC ! Convert altitude to C*5 string
    USE NAMGAS_FNC ! Return molecule name + (iso) associated with GASCOM index
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*),          INTENT(IN)  :: INPNAM ! Input (template) filename
    CHARACTER(*),          INTENT(OUT) :: OUTNAM ! Output filename
    INTEGER(I4), OPTIONAL, INTENT(IN)  :: IGAS ! Absorbing species index
    INTEGER(I4), OPTIONAL, INTENT(IN)  :: IJAC ! Jacobian index
    INTEGER(I4), OPTIONAL, INTENT(IN)  :: ILEV ! Level index
    INTEGER(I4), OPTIONAL, INTENT(IN)  :: IPSI ! Horiz.location index
    INTEGER(I4), OPTIONAL, INTENT(IN)  :: ISPC ! Spectral range index
    INTEGER(I4), OPTIONAL, INTENT(IN)  :: ITAN ! Tangent path index
!
! LOCAL VARIABLES
    INTEGER(I4)   :: L      ! Location for start of inserted strings
    CHARACTER(6)  :: PSISTR ! Horiz.Angle [deg] * 1000 written as string
    CHARACTER(LENJAC) :: JACSTR ! String containing Jacobian info
    CHARACTER(10) :: LEVSTR ! String containing Level dir & height info
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Find if there is a wildcard in the input name. 
  L = INDEX ( INPNAM, '*') 
  IF ( L .EQ. 0 ) THEN          ! No wildcard, so same filename
    OUTNAM = INPNAM
  ELSE                         ! Wildcard, so insert strings (if any)
    OUTNAM = INPNAM(1:L-1) // INPNAM(L+1:)    ! remove asterisk
!
! Insert any spectral range label
    IF ( PRESENT ( ISPC ) ) THEN
      OUTNAM = OUTNAM(1:L-1) // TRIM ( SPC(ISPC)%LAB ) // OUTNAM(L:)
      L = L + LEN_TRIM ( SPC(ISPC)%LAB )
    END IF
!
! Insert any absorbing species name
    IF ( PRESENT ( IGAS ) ) THEN
      IF ( IGAS .NE. 0 ) THEN
        OUTNAM = OUTNAM(1:L-1) // TRIM(NAMGAS(IGAS)) // OUTNAM(L:)
        L = L + LEN_TRIM ( NAMGAS(IGAS) ) 
      END IF
    END IF
!
! Insert any tangent path value
    IF ( PRESENT ( ITAN ) ) THEN
      OUTNAM = OUTNAM(1:L-1) // TRIM(TAN(ITAN)%STR) // OUTNAM(L:)
      L = L + LEN_TRIM ( TAN(ITAN)%STR ) 
    END IF
!
! Insert any horizontal angular coordinate value
    IF ( PRESENT ( IPSI ) ) THEN
      WRITE ( PSISTR, '(SP,I6.5)' ) NINT ( PSIGRA(IPSI) * 1000.0 ) 
      OUTNAM = OUTNAM(1:L-1) // TRIM(PSISTR) // OUTNAM(L:)
      L = L + LEN_TRIM ( PSISTR ) 
    END IF
!
! Insert any Jacobian identification
    IF ( PRESENT ( IJAC ) ) THEN
      IF ( IJAC .NE. 0 ) THEN
        JACSTR = JAC(IJAC)%COD 
        IF ( .NOT. JAC(IJAC)%COL ) &
          JACSTR = TRIM ( JACSTR ) // HGTSTR ( JAC(IJAC)%HGT ) 
        OUTNAM = OUTNAM(1:L-1) // '_' // TRIM ( JACSTR ) // OUTNAM(L:)
        L = L + 1 + LEN_TRIM ( JACSTR ) 
      END IF
    END IF
!
! Insert any Intermediate level identification
! Note that IDR refers to direction of integration from remote to observer
! so zenith view = 'up' is IDR=-1
    IF ( PRESENT ( ILEV ) ) THEN
      IF ( ILEV .NE. 0 ) THEN
        SELECT CASE ( LEV(ILEV)%IDR ) 
          CASE ( 1 )   ; LEVSTR = 'down' // HGTSTR ( LEV(ILEV)%HGT ) 
          CASE (-1 )   ; LEVSTR = 'up' // HGTSTR ( LEV(ILEV)%HGT ) 
          CASE DEFAULT ; LEVSTR = HGTSTR ( LEV(ILEV)%HGT ) 
        END SELECT
        OUTNAM = OUTNAM(1:L-1) // '_' // TRIM ( LEVSTR ) // OUTNAM(L:)
        L = L + 1 + LEN_TRIM ( LEVSTR ) 
      END IF
    END IF
!
  END IF
!
! Append Run ID (if any) 
  OUTNAM = TRIM ( OUTNAM ) // RUN_ID
!
END SUBROUTINE MAKNAM
END MODULE MAKNAM_SUB

