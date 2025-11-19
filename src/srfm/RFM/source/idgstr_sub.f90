MODULE IDGSTR_SUB
CONTAINS
SUBROUTINE IDGSTR ( GASSTR, IDXMOL, IDXISO ) 
!
! VERSION
!   20FEB23 AD Checked.
!   02OCT17 AD Original. Checked.
!
! DESCRIPTION
!   Decode string HITRAN/RFM index and isotope
!   Called by CHKGAS.
!   Returns IDXMOL = IDXISO = 0 if unable to decode
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE 
!
! ARGUMENTS
    CHARACTER(*), INTENT(IN)  :: GASSTR ! String to be decoded
    INTEGER(I4),  INTENT(OUT) :: IDXMOL ! HITRAN/RFM Molec#ID
    INTEGER(I4),  INTENT(OUT) :: IDXISO ! Isotope# (or 0)
!
! LOCAL VARIABLES
    INTEGER(I4)  :: IDOT ! Location of '.' or 'i' within GASSTR
    INTEGER(I4)  :: IOS  ! Value of IOSTAT after READ
    INTEGER(I4)  :: L    ! Length of string
! 
! EXECUTABLE CODE -------------------------------------------------------------
!
  L = LEN ( GASSTR ) 
!
  IDOT = INDEX ( GASSTR, '.' ) 
  IF ( IDOT .EQ. 0 ) IDOT = INDEX ( GASSTR, 'i' ) 
!
  IF ( IDOT .GT. 1 .AND. IDOT .LT. L ) THEN        ! Read molec#, isotope#
    READ ( GASSTR(1:IDOT-1), *, IOSTAT=IOS ) IDXMOL
    READ ( GASSTR(IDOT+1:L), *, IOSTAT=IOS ) IDXISO
  ELSE
    READ ( GASSTR, *, IOSTAT=IOS ) IDXMOL             ! Read molec#
    IDXISO = 0
  END IF
!
! If any read error, or IDXMOL not +ve, or IDXISO -ve, assume not valid string
  IF ( IOS .NE. 0 .OR. IDXMOL .LE. 0 .OR. IDXISO .LT. 0 ) THEN
    IDXMOL = 0
    IDXISO = 0
  END IF
!
END SUBROUTINE IDGSTR
END MODULE IDGSTR_SUB
