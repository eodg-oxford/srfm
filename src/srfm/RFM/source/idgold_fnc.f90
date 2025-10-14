MODULE IDGOLD_FNC
CONTAINS
INTEGER(I4) PURE FUNCTION IDGOLD ( IDGNEW )
!
! VERSION
!   08MAR23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Convert new RFM index for .xsc data to old value
!   General purpose module.
!   This is the inverse of function IDGNEW
!   Returns IDGOLD = 0 if IDGNEW not listed
!   Returns IDGOLD = IDGNEW if IDGNEW LE 49 - assumes unchanged
!   See molidx_sub.f90 for definitive list of new index assignments
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE 
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IDGNEW ! Old index
! 
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( IDGNEW .LE. 49 ) THEN
    IDGOLD = IDGNEW    ! no change, already old index
    RETURN
  END IF
!
  IDGOLD = 0
  SELECT CASE ( IDGNEW ) 
    CASE ( 100 ) ; IDGOLD = 50        ! aerosol
    CASE ( 111 ) ; IDGOLD = 51        ! f11
    CASE ( 112 ) ; IDGOLD = 52        ! f12
    CASE ( 117 ) ; IDGOLD = 53        ! f13
    CASE ( 118 ) ; IDGOLD = 54        ! f14
    CASE ( 121 ) ; IDGOLD = 55        ! f21
    CASE ( 122 ) ; IDGOLD = 56        ! f22
    CASE ( 113 ) ; IDGOLD = 57        ! f113
    CASE ( 114 ) ; IDGOLD = 58        ! f114
    CASE ( 115 ) ; IDGOLD = 59        ! f115
    CASE ( 104 ) ; IDGOLD = 60        ! ccl4
    CASE ( 101 ) ; IDGOLD = 61        ! clono2
    CASE ( 102 ) ; IDGOLD = 62        ! n2o5
    CASE ( 105 ) ; IDGOLD = 63        ! hno4
    CASE ( 103 ) ; IDGOLD = 64        ! sf6
!
    CASE ( 123 ) ; IDGOLD = 70        ! f123
    CASE ( 124 ) ; IDGOLD = 71        ! f124
    CASE ( 125 ) ; IDGOLD = 72        ! f141b
    CASE ( 126 ) ; IDGOLD = 73        ! f142b
    CASE ( 127 ) ; IDGOLD = 74        ! f225ca
    CASE ( 128 ) ; IDGOLD = 75        ! f225cb
    CASE ( 132 ) ; IDGOLD = 76        ! f32
    CASE ( 131 ) ; IDGOLD = 77        ! f125
    CASE ( 133 ) ; IDGOLD = 78        ! f134
    CASE ( 134 ) ; IDGOLD = 79        ! f134a
    CASE ( 135 ) ; IDGOLD = 80        ! f143a
    CASE ( 136 ) ; IDGOLD = 81        ! f152a
    CASE ( 116 ) ; IDGOLD = 82        ! f116
    CASE ( 106 ) ; IDGOLD = 83        ! sf5cf3
    CASE ( 145 ) ; IDGOLD = 84        ! pan
    CASE ( 142 ) ; IDGOLD = 85        ! ch3cn
    CASE ( 153 ) ; IDGOLD = 86        ! c6h6
    CASE ( 151 ) ; IDGOLD = 87        ! c2h6
    CASE ( 152 ) ; IDGOLD = 88        ! c3h8
    CASE ( 144 ) ; IDGOLD = 89        ! acetone
  END SELECT
!
END FUNCTION IDGOLD
END MODULE IDGOLD_FNC
