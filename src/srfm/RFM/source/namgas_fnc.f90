MODULE NAMGAS_FNC
CONTAINS
CHARACTER(LENGAS+2) PURE FUNCTION NAMGAS ( IGAS )
!
! VERSION
!   04MAR23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Return molecule name + i[#iso] associated with GASCOM index
!   General purpose module.
!
! GLOBAL DATA
    USE GASCOM_DAT ! Molecule and isotope data

  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER, INTENT(IN) :: IGAS ! Index of GAS in GASCOM 
!
! LOCAL VARIABLES
    CHARACTER(2)        :: ISO ! String containing isotope ID
    CHARACTER(LENGAS+2) :: NAM ! Local name of molecule (+iso)
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  NAM = GAS(IGAS)%COD
  IF ( ISOMOL ( GAS(IGAS)%IDM ) ) THEN
    WRITE ( ISO, '(A,I1)' )  'i', GAS(IGAS)%IDI
    NAM = TRIM ( NAM ) // ISO
  END IF
!
  NAMGAS = NAM      
!
END FUNCTION NAMGAS
END MODULE NAMGAS_FNC

