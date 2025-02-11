MODULE ATMGRA_SUB
CONTAINS
SUBROUTINE ATMGRA
!
! VERSION
!   05DEC23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Interpolate to fill 2-D atmospheric profile field
!   Called once by DRVATM (if GRA flag) and DRVNTE (if also NTE flag)
!   For each species (including pressure and temperature) the procedure is
!   to establish the psi locations at which profile values have been supplied
!   and interpolate values for any missing internal locations, and duplicate
!   end values for any missing external locations.
!   Horizontal interpolation is linear for T, log for p (linear if p small),
!   and depends upon LIN flag for VMR as in ATMLAY. 
!   The end result should be that full profiles for p, T and VMR are set at
!   all specified psi values, and any TVib profiles loaded via .atm files
!   Also loads reference profile (psi=0) into atmcom.inc
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GRACOM_DAT ! Atmospheric 2-D field
    USE ATMCOM_DAT ! Atmospheric profile data
!
! SUBROUTINES
    USE INTGRA_SUB ! Interpolate to fill 2-D atmospheric profile field
    USE MOVGRA_SUB ! Move profiles within GRACOM arrays
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IVIB ! Counter for vibrational temperatures
    INTEGER(I4) :: IVMR ! Counter for absorbing species
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL INTGRA ( PSIGRA, FLTGRA, .TRUE.,  TEMGRA ) ! Temperature
  CALL INTGRA ( PSIGRA, FLPGRA, .FALSE., PREGRA ) ! Pressure
!
  DO IVMR = 1, NVMR                               ! VMRs
    CALL INTGRA ( PSIGRA, FLVGRA(:,IVMR), LINVMR(IVMR), VMRGRA(:,:,IVMR) ) 
  END DO
!
  DO IVIB = 1, NVIB                               ! Vib Temps
    CALL INTGRA ( PSIGRA, FLNGRA(:,IVIB), .TRUE., VIBGRA(:,:,IVIB) ) 
  END DO
!
  CALL MOVGRA ( NOMGRA, 0 )    ! is this necessary?
!
END SUBROUTINE ATMGRA
END MODULE ATMGRA_SUB

