MODULE SPCCTM_SUB
CONTAINS
SUBROUTINE SPCCTM 
!
! VERSION
!   15MAR23 AD Add new H2O continuum v4.1 and C32 flag
!   02NOV18 AD New H2O continuum MT_CKD v3.2. Checked.
!   16NOV17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   Calculate continuum absorption
!   Called by RFMSPC if CTM flag enabled.
!   Calculates continuum at wide mesh points.
!   Continuum data available for absorbers H2O, O2, N2, CO2.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE GASCOM_DAT ! Molecule and isotope data
    USE IDXCON_DAT ! RFM/HITRAN indices for specific molecules
    USE FLGCOM_DAT, ONLY: C32FLG  ! T=use v3.2 H2O continuum
    USE WIDCOM_DAT, ONLY: NLBL, IDXLBL ! Indices for LBL path segments
!
! SUBROUTINES
    USE CTMC32_SUB ! MT_CKD v3.2 H2O continuum
    USE CTMC25_SUB ! H2O continuum MT_CKD v2.5
    USE CTMCKD_SUB ! CKD H2O continuum
    USE CTMCO2_SUB ! CO2 continuum
    USE CTMH2O_SUB ! H2O continuum
    USE CTMN2_SUB  ! N2 continuum
    USE CTMO2_SUB  ! O2 continuum
!
  IMPLICIT NONE
!
! LOCAL CONSTANTS
    LOGICAL, PARAMETER :: USECKD = .FALSE. ! T = use old CKD H2O continuum
    LOGICAL, PARAMETER :: USEC25 = .FALSE. ! T = use old MT_CKD v2.5 H2O ctm.
!                                          ! F = use new MT_CKD v3.2 H2O ctm.
! LOCAL VARIABLES
    INTEGER(I4) :: IGAS   ! Gas counter
    INTEGER(I4) :: ICLC   ! Gas/segment Path counter
    INTEGER(I4) :: ILBL   ! Counter for line-by-line path segments
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO ILBL = 1, NLBL
    ICLC = IDXLBL(ILBL)
    IGAS = CLC(ICLC)%IGS
    IF ( GAS(IGAS)%CTM ) THEN
      SELECT CASE ( GAS(IGAS)%IDM )
      CASE ( IDXH2O ) 
        IF ( USECKD ) THEN 
          CALL CTMCKD ( ILBL )
        ELSE IF ( USEC25 ) THEN
          CALL CTMC25 ( ILBL ) 
        ELSE IF ( C32FLG ) THEN 
          CALL CTMC32 ( ILBL )
        ELSE
          CALL CTMH2O ( ILBL )
        END IF
      CASE ( IDXCO2 ) 
        CALL CTMCO2 ( ILBL )
      CASE ( IDXO2 ) 
        CALL CTMO2 ( ILBL )
      CASE ( IDXN2 ) 
        CALL CTMN2 ( ILBL )
      END SELECT
    END IF
  END DO
!
END SUBROUTINE SPCCTM
END MODULE SPCCTM_SUB

