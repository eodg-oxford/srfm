MODULE ISOLST_SUB
CONTAINS
SUBROUTINE ISOLST ( IDXMOL, NISO, WGTISO )
!
! VERSION
!   27SEP22 AD Add HITRAN2020 Molecules/isotopes. Checked.
!   03JUN19 AD Add GEISA molecules/isotopes. Checked.
!   04MAY18 AD Add HITRAN2016 isotopes/molecules: add H2O-262, CO2-737, 
!              change#48 C3H8>C2N2, #50 BrO>SO, #51 GeH4>C3H4, #52 C3H8>CH3
!              #53 C2N2>CS2,                
!   01JUN17 AD F90 version.
!
! DESCRIPTION
!   List isotopomer weights of specific molecule
!   General purpose module.
!   If WGTISO argument is missing, just returns no.of defined isotopomers.
!   Returns NISO=0 if no isotopes defined for molecule
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
  SAVE
!
! ARGUMENTS
    INTEGER(I4),        INTENT(IN)  :: IDXMOL    ! HITRAN/RFM index of molecule
    INTEGER(I4),        INTENT(OUT) :: NISO      ! No. of isotopes
    REAL(R4), OPTIONAL, INTENT(OUT) :: WGTISO(:) ! List of isotopic weights
!
! LOCAL VARIABLES 
! Effectively these are constants but have to be stored as variables in order to
! use WGT as a pointer
!
! #1  H2O                             161  181  171  162  182  172  262
    REAL(R4), TARGET :: WGT01(7) = (/ 18., 20., 19., 19., 21., 20., 20./)
! #2  CO2                             626  636  628  627  638  637  828  827  727  838  837  737 
    REAL(R4), TARGET :: WGT02(12) = (/ 44., 45., 46., 45., 47., 46., 48., 47., 46., 49., 48., 47. /)
! #3  O3                             666  668  686  667  676
    REAL(R4), TARGET :: WGT03(5) = (/ 48., 50., 50., 49., 49. /)
! #4  N2O                            446  456  546  448  447  458* 548* 556*  *GEISA
    REAL(R4), TARGET :: WGT04(8) = (/ 44., 45., 45., 46., 45., 47., 47., 46. /)
! #5  CO                              26   36   28   27   38   37
    REAL(R4), TARGET :: WGT05(6) = (/ 28., 29., 30., 29., 31., 30. /)
! #6  CH4                            211  311  212  312
    REAL(R4), TARGET :: WGT06(4) = (/ 16., 17., 17., 18. /)
! #7  O2                              66   68   67
    REAL(R4), TARGET :: WGT07(3) = (/ 32., 34., 33. /)
! #8  NO                              46   56   48
    REAL(R4), TARGET :: WGT08(3) = (/ 30., 31., 32. /)
! #9  SO2                            626  646  628* 636*    *old GEISA, but keep in case required
    REAL(R4), TARGET :: WGT09(4) = (/ 64., 66., 66., 65. /)
! #10 NO2                            646  656*   *HITRAN2020
    REAL(R4), TARGET :: WGT10(2) = (/ 46., 47. /)
! #11 NH3                           4111 5111
    REAL(R4), TARGET :: WGT11(2) = (/ 17., 18. /)
! #12 HNO3                           146  156
    REAL(R4), TARGET :: WGT12(2) = (/ 63., 64. /)
! #13 OH                              61   81   62
    REAL(R4), TARGET :: WGT13(3) = (/ 17., 19., 18. /)
! #14 HF                              19   29
    REAL(R4), TARGET :: WGT14(2) = (/ 20., 21. /)
! #15 HCl                             15   17   25   27
    REAL(R4), TARGET :: WGT15(4) = (/ 36., 38., 37., 39. /)
! #16 HBr                             19   11   29   21
    REAL(R4), TARGET :: WGT16(4) = (/ 80., 82., 81., 83. /)
! #17 HI                                17    27
    REAL(R4), TARGET :: WGT17(2) =  (/ 128., 129. /)
! #18 ClO                             56   76
    REAL(R4), TARGET :: WGT18(2) = (/ 51., 53. /)
! #19 OCS                            622  624  632  623  822  634* *Geisa
    REAL(R4), TARGET :: WGT19(6) = (/ 60., 62., 61., 61., 62., 63. /)
! #20 H2CO                           126  136  128
    REAL(R4), TARGET :: WGT20(3) = (/ 30., 31., 32. /) 
! #21 HOCl                           165  167 
    REAL(R4), TARGET :: WGT21(2) = (/ 52., 54. /)
! #22 N2                              44   45
    REAL(R4), TARGET :: WGT22(2) = (/ 28., 29. /)
! #23 HCN                            124  134  125  224*    *Geisa
    REAL(R4), TARGET :: WGT23(4) = (/ 27., 28., 28., 28. /)
! #24 CH3Cl                          215  217
    REAL(R4), TARGET :: WGT24(2) = (/ 50., 52. /)
! #25 H2O2                          1661
    REAL(R4), TARGET :: WGT25(1) = (/ 34. /)
! #26 C2H2                          1221 1231 1222
    REAL(R4), TARGET :: WGT26(3) = (/ 26., 27., 27. /)
! #27 C2H6                          1221 1231
    REAL(R4), TARGET :: WGT27(2) = (/ 30., 31. /)
! #28 PH3                           1111
    REAL(R4), TARGET :: WGT28(1) = (/ 34. /)
! #29 COF2                           269  369
    REAL(R4), TARGET :: WGT29(2) = (/ 66., 67. /)
! #30 SF6                             29
    REAL(R4), TARGET :: WGT30(1) = (/ 146. /)
! #31 H2S                            121  141  131
    REAL(R4), TARGET :: WGT31(3) = (/ 34., 36., 35. /)
! #32 HCOOH                          126
    REAL(R4), TARGET :: WGT32(1) = (/ 46. /)               
! #33 HO2                            166
    REAL(R4), TARGET :: WGT33(1) = (/ 33. /)               
! #34 O                                6
    REAL(R4), TARGET :: WGT34(1) = (/ 16. /)               
! #35 ClONO2                        5646 7646
    REAL(R4), TARGET :: WGT35(2) = (/ 97., 99. /)          
! #36 NO+                             46
    REAL(R4), TARGET :: WGT36(1) = (/ 30. /)               
! #37 HOBr                           169  161
    REAL(R4), TARGET :: WGT37(2) = (/ 96., 98. /)               
! #38 C2H4                           221  231  
    REAL(R4), TARGET :: WGT38(2) = (/ 28., 29. /)               
! #39 CH3OH                         2161
    REAL(R4), TARGET :: WGT39(1) = (/ 32. /)
! #40 CH3Br                          219  211 
    REAL(R4), TARGET :: WGT40(2) = (/ 94., 96. /)
! #41 CH3CN                         2124
    REAL(R4), TARGET :: WGT41(1) = (/ 41. /)               
! #42 CF4                             29
    REAL(R4), TARGET :: WGT42(1) = (/ 88. /)               
! #43 C4H2                          2211
    REAL(R4), TARGET :: WGT43(1) = (/ 50. /)               
! #44 HC3N                          1224
    REAL(R4), TARGET :: WGT44(1) = (/ 51. /)
! #45 H2                             11   12 
    REAL(R4), TARGET :: WGT45(2) = (/ 2., 3. /)
! #46 CS                              22   24   32   23 
    REAL(R4), TARGET :: WGT46(4) = (/ 44., 46., 45., 45. /)               
! #47 SO3                             26
    REAL(R4), TARGET :: WGT47(1) = (/ 80. /)               
! #48 C2N2                          4224
    REAL(R4), TARGET :: WGT48(1) = (/ 52. /)
! #49 COCl2                         2655 2657 2677
    REAL(R4), TARGET :: WGT49(3) = (/ 98., 102., 104. /)
! #50 SO                              26   46   28 
    REAL(R4), TARGET :: WGT50(3) = (/ 48., 50., 50. /)
! #51 CH3F                           219 (HITRAN 2020)
    REAL(R4), TARGET :: WGT51(1) = (/ 34. /)
! #52 GeH4                           411  211  011  311  611 (HITRAN 2020)
    REAL(R4), TARGET :: WGT52(5) = (/ 78., 76., 74., 77., 80. /)
! #53 CS2                            222  224  223  232 (HITRAN 2020)
    REAL(R4), TARGET :: WGT53(4) = (/ 76., 78., 77., 77. /)
! #54 CH3I                            217 (HITRAN 2020)
    REAL(R4), TARGET :: WGT54(1) = (/ 142. /) 
! #55 NF3                             49 (HITRAN 2020)
    REAL(R4), TARGET :: WGT55(1) = (/ 71. /) 
! #56 C3H4                          1221 (TIPS 2020)
    REAL(R4), TARGET :: WGT56(1) = (/ 40. /)
! #57 CH3                           2111 (TIPS 2020)
    REAL(R4), TARGET :: WGT57(1) = (/ 15. /)
!
! GEISA molecules
! #60 GeH4                           411*
!    REAL(R4), TARGET :: WGT60(1) = (/ 77.,  /)           
! #61 C3H8                           221*
    REAL(R4), TARGET :: WGT61(1) = (/ 44. /)               
! #62 HNC                            142*
    REAL(R4), TARGET :: WGT62(1) = (/ 27. /)               
! #63 C6H6                           266*
    REAL(R4), TARGET :: WGT63(1) = (/ 78. /)
!
    REAL(R4), POINTER :: WGT(:)   ! used to point to appropriate WGTnn array
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  NULLIFY ( WGT ) 
!
  SELECT CASE ( IDXMOL )
    CASE ( 1 )  ; WGT => WGT01
    CASE ( 2 )  ; WGT => WGT02
    CASE ( 3 )  ; WGT => WGT03
    CASE ( 4 )  ; WGT => WGT04
    CASE ( 5 )  ; WGT => WGT05
    CASE ( 6 )  ; WGT => WGT06
    CASE ( 7 )  ; WGT => WGT07
    CASE ( 8 )  ; WGT => WGT08
    CASE ( 9 )  ; WGT => WGT09
    CASE ( 10 ) ; WGT => WGT10
    CASE ( 11 ) ; WGT => WGT11
    CASE ( 12 ) ; WGT => WGT12
    CASE ( 13 ) ; WGT => WGT13
    CASE ( 14 ) ; WGT => WGT14
    CASE ( 15 ) ; WGT => WGT15
    CASE ( 16 ) ; WGT => WGT16
    CASE ( 17 ) ; WGT => WGT17
    CASE ( 18 ) ; WGT => WGT18
    CASE ( 19 ) ; WGT => WGT19
    CASE ( 20 ) ; WGT => WGT20
    CASE ( 21 ) ; WGT => WGT21
    CASE ( 22 ) ; WGT => WGT22
    CASE ( 23 ) ; WGT => WGT23
    CASE ( 24 ) ; WGT => WGT24
    CASE ( 25 ) ; WGT => WGT25
    CASE ( 26 ) ; WGT => WGT26
    CASE ( 27 ) ; WGT => WGT27
    CASE ( 28 ) ; WGT => WGT28
    CASE ( 29 ) ; WGT => WGT29
    CASE ( 30 ) ; WGT => WGT30
    CASE ( 31 ) ; WGT => WGT31
    CASE ( 32 ) ; WGT => WGT32
    CASE ( 33 ) ; WGT => WGT33
    CASE ( 34 ) ; WGT => WGT34
    CASE ( 35 ) ; WGT => WGT35
    CASE ( 36 ) ; WGT => WGT36
    CASE ( 37 ) ; WGT => WGT37
    CASE ( 38 ) ; WGT => WGT38
    CASE ( 39 ) ; WGT => WGT39
    CASE ( 40 ) ; WGT => WGT40
    CASE ( 41 ) ; WGT => WGT41
    CASE ( 42 ) ; WGT => WGT42
    CASE ( 43 ) ; WGT => WGT43
    CASE ( 44 ) ; WGT => WGT44
    CASE ( 45 ) ; WGT => WGT45
    CASE ( 46 ) ; WGT => WGT46
    CASE ( 47 ) ; WGT => WGT47
    CASE ( 48 ) ; WGT => WGT48
    CASE ( 49 ) ; WGT => WGT49
    CASE ( 50 ) ; WGT => WGT50
    CASE ( 51 ) ; WGT => WGT51
    CASE ( 52 ) ; WGT => WGT52
    CASE ( 53 ) ; WGT => WGT53
    CASE ( 54 ) ; WGT => WGT54
    CASE ( 55 ) ; WGT => WGT55
    CASE ( 56 ) ; WGT => WGT56
    CASE ( 57 ) ; WGT => WGT57
!
!    CASE ( 60 ) ; WGT => WGT60
    CASE ( 61 ) ; WGT => WGT61
    CASE ( 62 ) ; WGT => WGT62
    CASE ( 63 ) ; WGT => WGT63
    CASE DEFAULT ; CONTINUE        ! Any other molecule index
  END SELECT
!
  IF ( ASSOCIATED ( WGT ) ) THEN
    IF ( PRESENT ( WGTISO ) ) THEN
      NISO = MIN ( SIZE ( WGT ), SIZE ( WGTISO ) ) ! limit to size of WGTISO
      WGTISO(1:NISO) = WGT(1:NISO)
    ELSE
      NISO = SIZE ( WGT ) 
    END IF
    NULLIFY ( WGT )
  ELSE
    NISO = 0
  END IF
!
END SUBROUTINE ISOLST
END MODULE ISOLST_SUB
