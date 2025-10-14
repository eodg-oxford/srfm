MODULE HGTSTR_FNC
CONTAINS
CHARACTER(LENTAN) PURE FUNCTION HGTSTR ( HGT ) 
!
! VERSION
!   22JAN24 AD Checked.
!   12OCT18 AD Bug#10: fix problem with overflow for large -ve angles
!   01MAY17 AD F90 orginal. Checked.
!
! DESCRIPTION
!   Convert altitude to C*6 string
!   General purpose module.
!   Used for constructing altitude component of RFM output filenames.
!   This converts any altitude in range -999:999km to metres, altitudes outside
!   this range are kept as km.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE TANCOM_DAT, ONLY:LENTAN ! Length of tan.ht. info part of filename
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: HGT ! Altitude [km]
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  IF ( HGT .GE. 0.0 ) THEN
    IF ( NINT ( HGT * 1000.0 ) .LE. 999999 ) THEN
      WRITE ( HGTSTR, '(I6.6)' ) NINT ( HGT * 1000.0 ) 
    ELSE
      WRITE ( HGTSTR, '(I6.6)' ) NINT ( HGT ) 
    ENDIF
  ELSE
    IF ( NINT ( HGT * 1000.0 ) .GE. -999999 ) THEN
      WRITE ( HGTSTR, '(I7.6)' ) NINT ( HGT * 1000.0 ) 
! For small negative values, there's a chance that rounding may lead to 0
! in which case ensure that there is a minus sign otherwise there will be a
! space in the output filename
      IF ( HGTSTR .EQ. ' 000000' ) HGTSTR = '-000000'
    ELSE
      WRITE ( HGTSTR, '(I7.6)' ) NINT ( HGT )
    ENDIF
  END IF
!
END FUNCTION HGTSTR
END MODULE HGTSTR_FNC

