MODULE MOLIDX_SUB
CONTAINS
SUBROUTINE MOLIDX ( IDX, MOL, ADDNEW )
!
! VERSION
!   27SEP22 AD HITRAN2020 modifications: #51 CH3F (was C3H4),
!              #52 GeH4 (was CH3), #53 CS2 (unchanged), #54 CH3I (new),
!              #55 NF3 (new), #56 C3H4 (reassigned), #57 CH3 (reassigned) 
!              #60 deassigned (was GeH4) 
!   09MAR20 AD Bug#26 check allocation status of MOLNEW
!              Also rewrite to remove fixed array size MAXMOL
!   12JUL18 AD Assign #60 GeH4, #61 C3H8, #62 HNC, #63 C6H6
!   09MAY18 AD Redefine #48,50-53, deassign #54-56.
!   23JUN17 AD Add ADDNEW argument to allow new molecules to be added
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Give molecule name for HITRAN/RFM index, or vice-versa
!   General purpose module.
!   If IDX <= 0 returns IDX for given MOL, or unchanged if MOL unrecognised 
!   If IDX > 0 returns MOL for given IDX, or ' ' if IDX is unrecognised
!   Also locally stores list of assignments for molecules which can be 
!   represented by both line and cross-section data, ie choice of 2 indices.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL CONSTANTS
    USE GASCOM_DAT, ONLY: LENGAS ! Maximum length of molecule name
!
! SUBROUTINES
    USE LOCASE_FNC ! Convert text string to lower case
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
  SAVE
!
! ARGUMENTS
    INTEGER(I4),  INTENT(INOUT) :: IDX ! Molecule index
    CHARACTER(*), INTENT(INOUT) :: MOL ! Molecule name
    LOGICAL, OPTIONAL, INTENT(IN) :: ADDNEW
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: IOFNEW = 184 ! Offset index for new molecules
    INTEGER(I4), PARAMETER :: MAXMOL = 30  ! No. different molecules to store
!
! LOCAL VARIABLES
    INTEGER(I4)         :: ICHK      ! 1=check source, 2=check & warn
    INTEGER(I4)         :: IMOL      ! Counter for listed molecules
    INTEGER(I4)         :: INEW      ! Counter for new molecules
    INTEGER(I4)         :: NMOL = 0  ! Number of listed molecules
    INTEGER(I4)         :: NNEW = 0  ! Number of new molecules added
    CHARACTER(LEN(MOL)) :: LMOL      ! Lower case version of MOL
    CHARACTER(20), ALLOCATABLE :: MOLNEW(:) ! New molecules
    CHARACTER(20), ALLOCATABLE :: NEWSAV(:) ! Saved NEWMOL during reallocation
!
    TYPE :: MOLTYP
      INTEGER(I4)       :: IDX     ! Stored HITRAN/RFM indices of mols.
      CHARACTER(LENGAS) :: NAM     ! List of molecules stored so far
    END TYPE MOLTYP
    TYPE(MOLTYP), ALLOCATABLE :: MOLLST(:) ! List of re/assigned molecules
    TYPE(MOLTYP), ALLOCATABLE :: MOLSAV(:) ! Saved MOLLST during reallocation
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  ICHK = 0
! HITRAN line molecules
  IF ( IDX .GT. 0 ) THEN
    SELECT CASE ( IDX ) 
    CASE ( 1 )   ; MOL = 'h2o' 
    CASE ( 2 )   ; MOL = 'co2' 
    CASE ( 3 )   ; MOL = 'o3'  ; ICHK = 1
    CASE ( 4 )   ; MOL = 'n2o' ; ICHK = 1
    CASE ( 5 )   ; MOL = 'co' 
    CASE ( 6 )   ; MOL = 'ch4'
    CASE ( 7 )   ; MOL = 'o2' 
    CASE ( 8 )   ; MOL = 'no' 
    CASE ( 9 )   ; MOL = 'so2' ; ICHK = 1
    CASE ( 10 )  ; MOL = 'no2' ; ICHK = 1
    CASE ( 11 )  ; MOL = 'nh3'
    CASE ( 12 )  ; MOL = 'hno3'
    CASE ( 13 )  ; MOL = 'oh' 
    CASE ( 14 )  ; MOL = 'hf' 
    CASE ( 15 )  ; MOL = 'hcl'
    CASE ( 16 )  ; MOL = 'hbr'
    CASE ( 17 )  ; MOL = 'hi' 
    CASE ( 18 )  ; MOL = 'clo'
    CASE ( 19 )  ; MOL = 'ocs'
    CASE ( 20 )  ; MOL = 'h2co' ; ICHK = 1
    CASE ( 21 )  ; MOL = 'hocl'
    CASE ( 22 )  ; MOL = 'n2' 
    CASE ( 23 )  ; MOL = 'hcn'
    CASE ( 24 )  ; MOL = 'ch3cl'
    CASE ( 25 )  ; MOL = 'h2o2'
    CASE ( 26 )  ; MOL = 'c2h2'
    CASE ( 27 )  ; MOL = 'c2h6' ; ICHK = 1
    CASE ( 28 )  ; MOL = 'ph3' 
    CASE ( 29 )  ; MOL = 'cof2'
    CASE ( 30 )  ; MOL = 'sf6'  ; ICHK = 2
    CASE ( 31 )  ; MOL = 'h2s' 
    CASE ( 32 )  ; MOL = 'hcooh'
    CASE ( 33 )  ; MOL = 'ho2'
    CASE ( 34 )  ; MOL = 'o' 
    CASE ( 35 )  ; MOL = 'clono2' ; ICHK = 2
    CASE ( 36 )  ; MOL = 'no+'
    CASE ( 37 )  ; MOL = 'hobr'
    CASE ( 38 )  ; MOL = 'c2h4'   ; ICHK = 1
    CASE ( 39 )  ; MOL = 'ch3oh'  ; ICHK = 2
    CASE ( 40 )  ; MOL = 'ch3br'  
    CASE ( 41 )  ; MOL = 'ch3cn'  ; ICHK = 2
    CASE ( 42 )  ; MOL = 'cf4'    ; ICHK = 2
    CASE ( 43 )  ; MOL = 'c4h2'
    CASE ( 44 )  ; MOL = 'hc3n'
    CASE ( 45 )  ; MOL = 'h2' 
    CASE ( 46 )  ; MOL = 'cs' 
    CASE ( 47 )  ; MOL = 'so3'
    CASE ( 48 )  ; MOL = 'c2n2'
    CASE ( 49 )  ; MOL = 'cocl2'
    CASE ( 50 )  ; MOL = 'so'
    CASE ( 51 )  ; MOL = 'ch3f'
    CASE ( 52 )  ; MOL = 'geh4'
    CASE ( 53 )  ; MOL = 'cs2'
    CASE ( 54 )  ; MOL = 'ch3i'
    CASE ( 55 )  ; MOL = 'nf3'
! TIPS molecules
    CASE ( 56 )  ; MOL = 'c3h4'
    CASE ( 57 )  ; MOL = 'ch3'
! GEISA molecules
!    CASE ( 60 )  ; MOL = 'geh4'
    CASE ( 61 )  ; MOL = 'c3h8'  ; ICHK = 2  ! Propane
    CASE ( 62 )  ; MOL = 'hnc'
    CASE ( 63 )  ; MOL = 'c6h6'  ; ICHK = 2  ! Benzene
!
    CASE ( 99 )  ; MOL = 'air'
!
! Heavy  molecules
    CASE ( 100 ) ; MOL = 'aerosol'
    CASE ( 101 ) ; MOL = 'clono2' ; ICHK = 1  ! Chlorine Nitrate
    CASE ( 102 ) ; MOL = 'n2o5'               ! DiNitrogen Pentoxide
    CASE ( 103 ) ; MOL = 'sf6'    ; ICHK = 1  ! Sulphur Hexafluoride
    CASE ( 104 ) ; MOL = 'ccl4'    ! Carbon Tetrachloride
    CASE ( 105 ) ; MOL = 'hno4'    ! Peroxynitric Acid
    CASE ( 106 ) ; MOL = 'sf5cf3'  !Trifluoromethyl Sulphur Pentaflouride 
    CASE ( 107 ) ; MOL = 'brono2'  ! Bromine Nitrate
    CASE ( 108 ) ; MOL = 'cloocl'  ! Chlorine Peroxide
    CASE ( 109 ) ; MOL = 'x109'    ! spare
    CASE ( 110 ) ; MOL = 'x110'    ! spare
! CFCs
    CASE ( 111 ) ; MOL = 'f11'     ! ccl3f
    CASE ( 112 ) ; MOL = 'f12'     ! ccl2f2
    CASE ( 113 ) ; MOL = 'f113'    ! c2cl3f3
    CASE ( 114 ) ; MOL = 'f114'    ! c2cl2f4
    CASE ( 115 ) ; MOL = 'f115'    ! c2clf5
    CASE ( 116 ) ; MOL = 'f116'    ! c2f6
    CASE ( 117 ) ; MOL = 'f13'     ! cclf3
    CASE ( 118 ) ; MOL = 'f14'     ; ICHK = 1  ! cf4
! HCFCs
    CASE ( 121 ) ; MOL = 'f21'     ! chcl2f
    CASE ( 122 ) ; MOL = 'f22'     ! chclf2
    CASE ( 123 ) ; MOL = 'f123'    ! chcl2cf3
    CASE ( 124 ) ; MOL = 'f124'    ! chclfcf3
    CASE ( 125 ) ; MOL = 'f141b'   ! ch3ccl2f
    CASE ( 126 ) ; MOL = 'f142b '  ! ch3cclf2
    CASE ( 127 ) ; MOL = 'f225ca'  ! chcl2cf2cf3
    CASE ( 128 ) ; MOL = 'f225cb'  ! cclf2cf2chclf
! HFCs
    CASE ( 131 ) ; MOL = 'f125'    ! chf2cf3
    CASE ( 132 ) ; MOL = 'f32'     ! ch2f2
    CASE ( 133 ) ; MOL = 'f134'    ! chf2chf2
    CASE ( 134 ) ; MOL = 'f134a'   ! cfh2cf3
! 01OCT13: Insert 135, f143, from GEISA and renumber f143a,f152a
    CASE ( 135 ) ; MOL = 'f143'    ! ch2fchf2 1,1,2-Trifluoroethane
    CASE ( 136 ) ; MOL = 'f143a'   ! cf3ch3   1,1,1-Trifluoroethane
    CASE ( 137 ) ; MOL = 'f152a'   ! ch3chf2
! New GEISA molecule
    CASE ( 138 ) ; MOL = 'f365mfc' ! CF3CH2CF2CH3, pentaflourobutane
! MeHCs
    CASE ( 141 ) ; MOL = 'ch3oh'   ; ICHK = 1 ! Methanol
    CASE ( 142 ) ; MOL = 'ch3cn'   ; ICHK = 1 ! Acetonitrile
    CASE ( 143 ) ; MOL = 'ch3cho'  ! Acetaldehyde
    CASE ( 144 ) ; MOL = 'acetone' ! ch3coch3
    CASE ( 145 ) ; MOL = 'pan'     ! ch3c(o)oono2
! NMCs
    CASE ( 151 ) ; MOL = 'c2h6'   ; ICHK = 2   ! Ethane
    CASE ( 152 ) ; MOL = 'c3h8'   ; ICHK = 1   ! Propane
    CASE ( 153 ) ; MOL = 'c6h6'   ; ICHK = 1   ! Benzene
! New GEISA molecules
    CASE ( 154 ) ; MOL = 'c2h2'  ; ICHK = 2 ! Acetylene
    CASE ( 155 ) ; MOL = 'c2h4'  ; ICHK = 2 ! Ethylene
! Halocarbons
    CASE ( 161 ) ; MOL = 'c4f8'    ! Octafluorocyclobutane
! UV Cross-sections
    CASE ( 171 ) ; MOL = 'o3'   ; ICHK = 2
    CASE ( 172 ) ; MOL = 'n2o'  ; ICHK = 2
    CASE ( 173 ) ; MOL = 'so2'  ; ICHK = 2
    CASE ( 174 ) ; MOL = 'no2'  ; ICHK = 2
    CASE ( 175 ) ; MOL = 'h2co' ; ICHK = 2
    CASE ( 176 ) ; MOL = 'bro'  ; ICHK = 1
    CASE ( 177 ) ; MOL = 'no3'
    CASE ( 178 ) ; MOL = 'oclo'
    CASE ( 181 ) ; MOL = 'c7h8'
    CASE ( 182 ) ; MOL = 'oxylene'               ! o-C8H10
    CASE ( 183 ) ; MOL = 'mxylene'               ! m-C8H10
    CASE ( 184 ) ; MOL = 'pxylene'               ! p-C8H10
    CASE ( IOFNEW+1: ) 
      IF ( IDX-IOFNEW .LE. NNEW ) THEN
        MOL = MOLNEW(IDX-IOFNEW)
      ELSE
        MOL = ''
      END IF
    CASE DEFAULT ; MOL = ''
    END SELECT
    LMOL = MOL       ! MOL set locally to lower case already
!
! Return index for given molecule name
  ELSE
    LMOL = LOCASE ( MOL )     ! MOL on input could be upper case
    SELECT CASE ( LMOL ) 
    CASE ( 'h2o'     ) ; IDX = 1
    CASE ( 'co2'     ) ; IDX = 2
    CASE ( 'o3'      ) ; IDX = 3  ; ICHK = 1
    CASE ( 'n2o'     ) ; IDX = 4  ; ICHK = 1
    CASE ( 'co'      ) ; IDX = 5
    CASE ( 'ch4'     ) ; IDX = 6
    CASE ( 'o2'      ) ; IDX = 7
    CASE ( 'no'      ) ; IDX = 8
    CASE ( 'so2'     ) ; IDX = 9  ; ICHK = 1
    CASE ( 'no2'     ) ; IDX = 10 ; ICHK = 1
    CASE ( 'nh3'     ) ; IDX = 11
    CASE ( 'hno3'    ) ; IDX = 12
    CASE ( 'oh'      ) ; IDX = 13
    CASE ( 'hf'      ) ; IDX = 14
    CASE ( 'hcl'     ) ; IDX = 15
    CASE ( 'hbr'     ) ; IDX = 16
    CASE ( 'hi'      ) ; IDX = 17
    CASE ( 'clo'     ) ; IDX = 18
    CASE ( 'ocs'     ) ; IDX = 19
    CASE ( 'h2co'    ) ; IDX = 20 ; ICHK = 1
    CASE ( 'hocl'    ) ; IDX = 21
    CASE ( 'n2'      ) ; IDX = 22
    CASE ( 'hcn'     ) ; IDX = 23
    CASE ( 'ch3cl'   ) ; IDX = 24
    CASE ( 'h2o2'    ) ; IDX = 25
    CASE ( 'c2h2'    ) ; IDX = 26
    CASE ( 'c2h6'    ) ; IDX = 27 ; ICHK = 1
    CASE ( 'ph3'     ) ; IDX = 28
    CASE ( 'cof2'    ) ; IDX = 29
    CASE ( 'sf6q'    ) ; IDX = 30 ; MOL = 'sf6' ; ICHK = 2
    CASE ( 'h2s'     ) ; IDX = 31
    CASE ( 'hcooh'   ) ; IDX = 32
    CASE ( 'ho2'     ) ; IDX = 33
    CASE ( 'o'       ) ; IDX = 34
    CASE ( 'clono2q' ) ; IDX = 35 ; MOL = 'clono2' ; ICHK =2
    CASE ( 'no+'     ) ; IDX = 36
    CASE ( 'hobr'    ) ; IDX = 37
    CASE ( 'c2h4'    ) ; IDX = 38 ; ICHK = 1
    CASE ( 'ch3ohq'  ) ; IDX = 39 ; MOL = 'ch3oh' ; ICHK = 2
    CASE ( 'ch3br'   ) ; IDX = 40  
    CASE ( 'ch3cnq'  ) ; IDX = 41 ; MOL = 'ch3cn' ; ICHK = 2
    CASE ( 'cf4q', 'f14q' ) ; IDX = 42 ; MOL = 'f14' ; ICHK = 2
    CASE ( 'c4h2'    ) ; IDX = 43
    CASE ( 'hc3n'    ) ; IDX = 44
    CASE ( 'h2'      ) ; IDX = 45
    CASE ( 'cs'      ) ; IDX = 46
    CASE ( 'so3'     ) ; IDX = 47
    CASE ( 'c2n2'    ) ; IDX = 48 
    CASE ( 'cocl2'   ) ; IDX = 49
! TIPS molecules
    CASE ( 'so'      ) ; IDX = 50
    CASE ( 'ch3f'    ) ; IDX = 51
    CASE ( 'geh4'    ) ; IDX = 52
    CASE ( 'cs2'     ) ; IDX = 53
    CASE ( 'ch3i'    ) ; IDX = 54
    CASE ( 'nf3'     ) ; IDX = 55
    CASE ( 'c3h4'    ) ; IDX = 56
    CASE ( 'ch3'     ) ; IDX = 57
! GEISA molecules
!    CASE ( 'geh4'    ) ; IDX = 60
    CASE ( 'c3h8q'   ) ; IDX = 61 ; MOL = 'c3h8' ; ICHK = 2
    CASE ( 'hnc'     ) ; IDX = 62
    CASE ( 'c6h6q'   ) ; IDX = 63 ; MOL = 'c6h6' ; ICHK = 2
!
    CASE ( 'air' )     ; IDX = 99
!
! Heavy molecules
    CASE ( 'aerosol' ) ; IDX = 100
    CASE ( 'clono2'  ) ; IDX = 101 ; ICHK = 1     ! Chlorine nitrate
    CASE ( 'n2o5'    ) ; IDX = 102                ! DiNitrogen Pentoxide
    CASE ( 'sf6'     ) ; IDX = 103 ; ICHK = 1     ! Sulphur Hexafluoride
    CASE ( 'ccl4'    ) ; IDX = 104     ! Carbon Tetrachloride
    CASE ( 'hno4'    ) ; IDX = 105     ! Peroxynitric Acid
    CASE ( 'sf5cf3'  ) ; IDX = 106     ! Trifluoromethyl Sulphur Pentaflouride 
    CASE ( 'brono2'  ) ; IDX = 107     ! Bromine Nitrate
    CASE ( 'cloocl'  ) ; IDX = 108     ! Chlorine Peroxide
    CASE ( 'x109'    ) ; IDX = 109     ! spare
    CASE ( 'x110'    ) ; IDX = 110     ! spare
!
! CFCs
    CASE ( 'ccl3f',   'cfcl3',  'f11' ) ; IDX = 111
    CASE ( 'ccl2f2',  'cf2cl2', 'f12' ) ; IDX = 112
    CASE ( 'c2cl3f3', 'f113' )          ; IDX = 113
    CASE ( 'c2cl2f4', 'f114' )          ; IDX = 114
    CASE ( 'c2clf5',  'f115' )          ; IDX = 115
    CASE ( 'c2f6',    'f116' )          ; IDX = 116
    CASE ( 'cclf3',   'f13'  )          ; IDX = 117
    CASE ( 'cf4',     'f14'  )          ; IDX = 118 ; MOL = 'f14' ; ICHK = 1
!
! HCFCs
    CASE ( 'chcl2f',   'f21'   ) ; IDX = 121
    CASE ( 'chclf2',   'f22'   ) ; IDX = 122
    CASE ( 'chcl2cf3', 'f123'  ) ; IDX = 123
    CASE ( 'chclfcf3', 'f124'  ) ; IDX = 124
    CASE ( 'ch3ccl2f', 'f141b' ) ; IDX = 125
    CASE ( 'ch3cclf2', 'f142b' ) ; IDX = 126

    CASE ( 'chcl2cf2cf3',   'f225ca' ) ; IDX = 127
    CASE ( 'cclf2cf2chclf', 'f225cb' ) ; IDX = 128
!
! HFCs
    CASE ( 'chf2cf3',  'f125'  ) ; IDX = 131
    CASE ( 'ch2f2',    'f32'   ) ; IDX = 132
    CASE ( 'chf2chf2', 'f134'  ) ; IDX = 133
    CASE ( 'cfh2cf3',  'f134a' ) ; IDX = 134
    CASE (             'f143'  ) ; IDX = 135
    CASE ( 'cf3ch3',   'f143a' ) ; IDX = 136
    CASE ( 'ch3chf2',  'f152a' ) ; IDX = 137
    CASE (           'f365mfc' ) ; IDX = 138
!
! MeHCs
    CASE ( 'ch3oh'    ) ; IDX = 141 ; ICHK = 1    ! Methanol
    CASE ( 'ch3cn'    ) ; IDX = 142 ; ICHK = 1    ! Acetonitrile
    CASE ( 'ch3cho'   ) ; IDX = 143               ! Acetaldehyde
    CASE ( 'ch3coch3', 'acetone' ) ; IDX = 144    ! Acetone
    CASE ( 'ch3c(o)oono2', 'pan' ) ; IDX = 145    ! PAN
!
! NMCs
    CASE ( 'c2h6x'  ) ; IDX = 151 ; MOL = 'c2h6' ; ICHK = 2 ! Ethane
    CASE ( 'c3h8'   ) ; IDX = 152 ; MOL = 'c3h8' ; ICHK = 1 ! Propane
    CASE ( 'c6h6'   ) ; IDX = 153 ; ICHK = 1                ! Benzene
    CASE ( 'c2h2x'  ) ; IDX = 154 ; MOL = 'c2h2' ; ICHK = 2 ! Acetylene
    CASE ( 'c2h4x'  ) ; IDX = 155 ; MOL = 'c2h5' ; ICHK = 2 ! Ethylene
!
! Halocarbons
    CASE ( 'c4f8' ) ; IDX = 161         ! Octafluorocyclobutane
!
! UV Cross-sections
    CASE ( 'o3x'   ) ; IDX = 171 ; MOL = 'o3'   ; ICHK = 2
    CASE ( 'n2ox'  ) ; IDX = 172 ; MOL = 'n2o'  ; ICHK = 2
    CASE ( 'so2x'  ) ; IDX = 173 ; MOL = 'so2'  ; ICHK = 2
    CASE ( 'no2x'  ) ; IDX = 174 ; MOL = 'no2'  ; ICHK = 2
    CASE ( 'h2cox' ) ; IDX = 175 ; MOL = 'h2co' ; ICHK = 2
    CASE ( 'bro'   ) ; IDX = 176 ; ICHK = 1
    CASE ( 'no3'   ) ; IDX = 177
    CASE ( 'oclo'  ) ; IDX = 178
    CASE ( 'c7h8'  ) ; IDX = 181           ! Toluene
    CASE ( 'oxylene' ) ; IDX = 182        ! o-C8H10
    CASE ( 'mxylene' ) ; IDX = 183        ! m-C8H10
    CASE ( 'pxylene' ) ; IDX = 184        ! p-C8H10
!
    CASE DEFAULT 
      IDX = 0
      DO INEW = 1, NNEW
        IF ( LMOL .EQ. MOLNEW(INEW) ) THEN
          IDX = INEW + IOFNEW 
          EXIT
        END IF
      END DO
      IF ( IDX .EQ. 0 .AND. PRESENT ( ADDNEW ) ) THEN  ! Assign new molecule
        IF ( ADDNEW ) THEN
          IF ( ALLOCATED ( MOLNEW ) ) CALL MOVE_ALLOC ( MOLNEW, NEWSAV )
          NNEW = NNEW + 1
          ALLOCATE ( MOLNEW(NNEW) ) 
          IF ( ALLOCATED ( NEWSAV ) ) MOLNEW(1:NNEW-1) = NEWSAV
          MOLNEW(NNEW) = LMOL
          IDX = NNEW + IOFNEW
          CALL WRTLOG ( 'W-MOLIDX: Assuming user-defined .xsc molecule: ' & 
                        // LMOL ) 
        END IF
      END IF
    END SELECT
  END IF
!
  IF ( ICHK .GT. 0 ) THEN 
! See if source already assigned for this molecule
    DO IMOL = 1, NMOL 
      IF ( LMOL .EQ. MOLLST(IMOL)%NAM ) THEN
        IDX = MOLLST(IMOL)%IDX
        RETURN
      END IF
    END DO
!
! New molecule so assign index
    IF ( ALLOCATED ( MOLLST ) ) CALL MOVE_ALLOC ( MOLLST, MOLSAV ) 
    NMOL = NMOL + 1
    ALLOCATE ( MOLLST(NMOL) ) 
    IF ( ALLOCATED ( MOLSAV ) ) MOLLST(1:NMOL-1) = MOLSAV
    MOLLST(NMOL)%NAM = LMOL
    MOLLST(NMOL)%IDX = IDX
!
    IF ( ICHK .EQ. 2 ) THEN
      IF ( IDX .LT. 100 ) THEN
        CALL WRTLOG ( 'W-MOLIDX: changing to use line data for gas=' // MOL )
      ELSE
        CALL WRTLOG ( 'W-MOLIDX: changing to use x/s data for gas=' // MOL )
      END IF
    END IF
  END IF
!
END SUBROUTINE MOLIDX
END MODULE MOLIDX_SUB

