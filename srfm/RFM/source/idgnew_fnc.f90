MODULE IDGNEW_FNC
CONTAINS
INTEGER(I4) PURE FUNCTION IDGNEW ( IDGOLD )
!
! VERSION
!   05MAR23 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Convert old RFM index for .xsc data to new value
!   General purpose module.
!   This is the inverse of function IDGOLD
!   Returns IDGNEW = 0 if IDGOLD not listed.
!   Returns IDGNEW = IDGOLD if IDGOLD GE 100 (assumes already new index)
!   or if IDGOLD LE 49
!   See MOLIDX for definitive list of new index assignments
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE 
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IDGOLD ! Old index
! 
! EXECUTABLE CODE -------------------------------------------------------------
!
  SELECT CASE ( IDGOLD ) 
  CASE ( 50 ) ; IDGNEW = 100        ! aerosol
  CASE ( 51 ) ; IDGNEW = 111        ! f11
  CASE ( 52 ) ; IDGNEW = 112        ! f12
  CASE ( 53 ) ; IDGNEW = 117        ! f13
  CASE ( 54 ) ; IDGNEW = 118        ! f14
  CASE ( 55 ) ; IDGNEW = 121        ! f21
  CASE ( 56 ) ; IDGNEW = 122        ! f22
  CASE ( 57 ) ; IDGNEW = 113        ! f113
  CASE ( 58 ) ; IDGNEW = 114        ! f114
  CASE ( 59 ) ; IDGNEW = 115        ! f115
  CASE ( 60 ) ; IDGNEW = 104        ! ccl4
  CASE ( 61 ) ; IDGNEW = 101        ! clono2
  CASE ( 62 ) ; IDGNEW = 102        ! n2o5
  CASE ( 63 ) ; IDGNEW = 105        ! hno4
  CASE ( 64 ) ; IDGNEW = 103        ! sf6
  CASE ( 70 ) ; IDGNEW = 123        ! f123
  CASE ( 71 ) ; IDGNEW = 124        ! f124
  CASE ( 72 ) ; IDGNEW = 125        ! f141b
  CASE ( 73 ) ; IDGNEW = 126        ! f142b
  CASE ( 74 ) ; IDGNEW = 127        ! f225ca
  CASE ( 75 ) ; IDGNEW = 128        ! f225cb
  CASE ( 76 ) ; IDGNEW = 132        ! f32
  CASE ( 77 ) ; IDGNEW = 131        ! f125
  CASE ( 78 ) ; IDGNEW = 133        ! f134
  CASE ( 79 ) ; IDGNEW = 134        ! f134a
  CASE ( 80 ) ; IDGNEW = 135        ! f143a
  CASE ( 81 ) ; IDGNEW = 136        ! f152a
  CASE ( 82 ) ; IDGNEW = 116        ! f116
  CASE ( 83 ) ; IDGNEW = 106        ! sf5cf3
  CASE ( 84 ) ; IDGNEW = 145        ! pan
  CASE ( 85 ) ; IDGNEW = 142        ! ch3cn
  CASE ( 86 ) ; IDGNEW = 153        ! c6h6
  CASE ( 87 ) ; IDGNEW = 151        ! c2h6
  CASE ( 88 ) ; IDGNEW = 152        ! c3h8
  CASE ( 89 ) ; IDGNEW = 144        ! acetone
  CASE DEFAULT ; IDGNEW = IDGOLD
  END SELECT
!
END FUNCTION IDGNEW
END MODULE IDGNEW_FNC
