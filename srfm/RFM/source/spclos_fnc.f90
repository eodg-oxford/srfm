MODULE SPCLOS_FNC
CONTAINS
FUNCTION SPCLOS ( ITAN, TYP ) 
!
! VERSION
!   29MAR24 AD Checked.
!   04FEB19 AD Original. Adapted from earlier subroutine. Checked.
!
! DESCRIPTION
!   Calculate LOS Jacobian spectrum
!   Called by SPCOUT if LOS flag enabled
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FULCOM_DAT ! Full grid data    
    USE TANCOM_DAT ! Tangent path data
    USE PTBCON_DAT, ONLY: PTBLOS ! LOS perturbation [km]
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),  INTENT(IN) :: ITAN ! Nominal limb tangent path
    CHARACTER(3), INTENT(IN) :: TYP  ! 'OPT', 'RAD' or 'TRA' spectra type
!
! FUNCTION TYPE
    REAL(R8) :: SPCLOS(NFUL)   ! Return spectrum size NFUL
!
! LOCAL VARIABLES
    INTEGER(I4) :: I1,I2,I3    ! Indices of tan.paths used for interpolation
    INTEGER(I4) :: IDX(3)      ! Indices for 3-point fit
    REAL(R8)    :: W           ! Weight for linear interpolation
    REAL(R8)    :: WGT(3)      ! Weights for quadratic interpolation 
    REAL(R8)    :: Z21,Z32,Z31 ! Differences between tangent heights
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  SELECT CASE ( NTAN ) 
  CASE ( 1 )                                ! NTAN=1 not allowed for LOS flag
    STOP 'F-SPCLOS: Logical error#1'
  CASE ( 2 ) 
    I1 = MAX ( ITAN - 1, 1 ) 
    I3 = MIN ( ITAN + 1, NTAN ) 
    W = PTBLOS / ( TAN(I3)%USR - TAN(I1)%USR ) 
    SELECT CASE ( TYP )
      CASE ( 'OPT' ) ; SPCLOS = ( OPTFUL(:,I3) - OPTFUL(:,I1) ) * W
      CASE ( 'RAD' ) ; SPCLOS = ( RADFUL(:,I3) - RADFUL(:,I1) ) * W
      CASE ( 'TRA' ) ; SPCLOS = ( TRAFUL(:,I3) - TRAFUL(:,I1) ) * W
      CASE DEFAULT   ; STOP 'F-SPCLOS: Logical error#2'
    END SELECT
  CASE ( 3: ) 
    I1 = ITAN-1
    I2 = ITAN
    I3 = ITAN+1
    IF ( ITAN .EQ. 1 ) I1 = 3
    IF ( ITAN .EQ. NTAN ) I3 = NTAN-2
    Z21 = ( TAN(I2)%USR - TAN(I1)%USR ) / PTBLOS
    Z32 = ( TAN(I3)%USR - TAN(I2)%USR ) / PTBLOS
    Z31 = Z32 + Z21
    WGT(1) = -Z32 / Z21 / Z31 
    WGT(2) = (Z32-Z21) / Z21 / Z32
    WGT(3) = Z21 / Z31 / Z32
    IDX = (/ I1, I2, I3 /)
    SELECT CASE ( TYP )
      CASE ( 'OPT' ) ; SPCLOS = MATMUL ( OPTFUL(:,IDX),  WGT )
      CASE ( 'RAD' ) ; SPCLOS = MATMUL ( RADFUL(:,IDX),  WGT )
      CASE ( 'TRA' ) ; SPCLOS = MATMUL ( TRAFUL(:,IDX),  WGT )
      CASE DEFAULT   ; STOP 'F-SPCLOS: Logical error#3'
    END SELECT
  END SELECT
!
END FUNCTION SPCLOS
END MODULE SPCLOS_FNC
