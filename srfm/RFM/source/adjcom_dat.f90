MODULE ADJCOM_DAT
!
! VERSION
!   26DEC23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Path-adjusted line data 
!   Set by module ADJUST. 
!
! VARIABLE KINDS
    USE KIND_DAT 
!
  IMPLICIT NONE
  SAVE
!
! GLOBAL VARIABLES
    REAL(R4) :: PPAADJ ! Partial pressure [atm]
    REAL(R4) :: PREADJ ! Pressure [atm]
    REAL(R4) :: TEMADJ ! Temperature [K]
    REAL(R4) :: WIDADJ ! half width [cm-1]
    REAL(R4) :: STRADJ ! line strength  [cm-1]
    REAL(R4) :: DOPADJ ! Doppler half width [cm-1]
    REAL(R4) :: YMXADJ ! Y mixing coeffient
!
    REAL(R8) :: WNOADJ ! Wavenumber [cm-1]
!
END MODULE ADJCOM_DAT

