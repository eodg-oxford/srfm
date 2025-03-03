MODULE DIMPCG_SUB
CONTAINS
SUBROUTINE DIMPCG 
!
! VERSION
!   19OCT23 AD Checked.
!   01MAY17 AD F90 conversion of part of dimpre.for. Tested
!
! DESCRIPTION
!   Set TAB p-axis to CG values of internal ATM profile
!   Called by DRVDIM if TAB flag and PCG present in *DIM section of driver file.
!   The CG pressure calculation is extracted from subroutine VRTPTH
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE TABCOM_DAT, ONLY: NPTAB,PAXTAB ! Pressure Axis for tab. abs. coeff.
!
! SUBROUTINES
    USE C9REAL_GEN ! Write real number as C*9 string
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! LOCAL CONSTANTS
! Note: thin layers can cause problems in evaluating DPV etc, so avoid
! 1E-3 = 0.1%, about 7m min layer thickness
    REAL(R4), PARAMETER :: DLPMIN = 1.0E-3 ! Min val for dlnp for explicit diff
!
! LOCAL VARIABLES
    INTEGER(I4) :: IATM       ! Atmos. layer# for path
    INTEGER(I4) :: JATM       ! Index of profile point at top of segment
    REAL(R4)    :: B          ! VMR profile parameter: VMR = VMR1*(P/P1)^X
    REAL(R4)    :: DLNP       ! ln(p_l/p_u)
    REAL(R4)    :: DP         ! p_l - p_u
    REAL(R4)    :: DP2        ! p_l^2 - p_u^2
    REAL(R4)    :: DP2V       ! p_l^2 v_l - p_u^2 v_u
    REAL(R4)    :: DP2V2      ! p_l^2 v_l^2 - p_u^2 v_u^2
    REAL(R4)    :: DPV        ! p_l v_l - p_u v_l
    REAL(R4)    :: DV         ! v_l - v_u
    REAL(R4)    :: PCG        ! Curtis-Godson pressure [mb]
    REAL(R4)    :: PRE1, PRE2 ! Pressure [mb] at upp,low segment bounds
    REAL(R4)    :: USUM       ! Absorber Amount in layer [10^5*kmole/m^2]
    REAL(R4)    :: VMR1, VMR2 ! Vmr [ppv] at upp,low seg.bounds
!
! EXECUTABLE CODE --------------------------------------------------------------
!
! Determine no.of p-axis grid points that will be set
  CALL WRTLOG ( 'I-DIMPCG: Setting TAB p-axis to CG pressure ' // &
                'values of internal .atm profile' )
  NPTAB = NATM + 1       ! NATM-1 layer values plus max+min level values
  ALLOCATE ( PAXTAB(NPTAB) )
!
  PAXTAB(1) = PREATM(1)    ! Set first value at max pressure level
  VMR2 = VMRATM(1,1)
  PRE2 = PREATM(1)
  DO IATM = 1, NATM-1
    JATM = IATM + 1
! Copy previous segment upper value to current segment lower values
    VMR1 = VMR2
    PRE1 = PRE2
! Set current segment upper values from profile layer values
    VMR2 = VMRATM(JATM,1) 
    PRE2 = PREATM(JATM)
! Calculate summation terms for CG integrals (allow for VMR1=0 and/or VMR2=0)
    DLNP  = LOG(PRE1/PRE2)
    DP2V  = PRE1*PRE1*VMR1 - PRE2*PRE2*VMR2
    DP2V2 = PRE1*PRE1*VMR1*VMR1 - PRE2*PRE2*VMR2*VMR2
    DPV   = PRE1*VMR1 - PRE2*VMR2
    IF ( LINVMR(1) .OR. VMR1 .LT. TINY(1.0) .OR. VMR2 .LT. TINY(1.0) ) THEN
      DP    = PRE1 - PRE2
      DP2   = PRE1*PRE1 - PRE2*PRE2
      DV    = VMR1 - VMR2
      USUM  = DPV - DP * DV / DLNP
      IF ( USUM .NE. 0.0 .AND. DLNP .GE. DLPMIN ) THEN
        PCG = 0.5 * ( DP2V - 0.5*DP2*DV/DLNP ) / USUM
      ELSE 
        PCG = SQRT ( PRE1 * PRE2 )
      END IF
    ELSE
      B = LOG(VMR1/VMR2) / DLNP
      IF ( ABS ( 1.0 + B ) .LE. 1.0E-3 .OR. DLNP .LT. DLPMIN ) THEN
        USUM = PRE1*VMR1*DLNP
        PCG  = DP2V / ( 2.0 + B ) / USUM
      ELSE 
        USUM = DPV / ( 1.0 + B )
        IF ( ABS ( 2.0 + B ) .LE. 1.0E-3 ) THEN
          PCG = PRE1*PRE1*VMR1*DLNP
        ELSE
          PCG = DP2V / ( 2.0 + B ) / USUM
        END IF
      END IF
    END IF
    PAXTAB(JATM) = PCG
  END DO
  PAXTAB(NATM+1) = PREATM(NATM)   ! Set last point at min pressure level
!
! Print values to log file
  CALL WRTLOG ( 'p-axis values: ', .TRUE. )
  DO IATM = 1, NPTAB
    CALL WRTLOG ( C9REAL ( PAXTAB(IATM) ), .TRUE. )
  END DO
  CALL WRTLOG ( '', .FALSE. ) 
!
END SUBROUTINE DIMPCG
END MODULE DIMPCG_SUB
