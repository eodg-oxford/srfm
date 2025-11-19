MODULE DSHVAL_FNC
CONTAINS
REAL(R8) FUNCTION DSHVAL ( HGT, PSI )
!
! VERSION
!   16JAN24 AD Checked.
!   01MAY17 AD F90 version. Checked.
!
! DESCRIPTION
!   Interpolate Density Scale Height
!   General purpose function
!   Density Scale Height applies to each layer rather than level, so no interp.
!   in height.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE ATMCOM_DAT ! Atmospheric profile data
    USE GRACOM_DAT ! Atmospheric 2-D field
!
! SUBROUTINES
    USE INTRVL_SUB ! Calculate interpolation interval
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4),           INTENT(IN) :: HGT  ! Altitude [km]
    REAL(R4), OPTIONAL, INTENT(IN) :: PSI ! Horiz.angle [deg]
!
! LOCAL VARIABLES
    INTEGER(I4) :: IATM ! Lower vertical profile index
    INTEGER(I4) :: IPSI ! Lower horizontal profile index
    REAL(R4)    :: DATM ! (Dummy) Fraction of vert intrvl of HGT above IATM
    REAL(R4)    :: DPSI ! Fraction of horiz interval of PSI above IPSI
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  CALL INTRVL ( HGTATM, HGT, IATM, DATM )  
!
  IF ( PRESENT ( PSI ) .AND. NPSI .GT. 1 ) THEN
    CALL INTRVL ( PSIGRA, PSI, IPSI, DPSI )  
    DSHVAL = DBLE ( (1.0-DPSI)*DSHGRA(IATM,IPSI) + DPSI*DSHGRA(IATM,IPSI+1) )
  ELSE 
    DSHVAL = DBLE ( DSHATM(IATM) )
  END IF
!
END FUNCTION DSHVAL
END MODULE DSHVAL_FNC


