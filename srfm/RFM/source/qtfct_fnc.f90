MODULE QTFCT_FNC
CONTAINS
REAL(R4) FUNCTION QTFCT ( IDGAS, ISO, TEM )
!
! VERSION
!   17JUL24 AD Checked.
!   22JAN18 AD Adapted for new (HITRAN2016) structure of TIPS data. 
!   01MAY17 AD F90 conversion, originally a subroutine. Checked.
!
! DESCRIPTION
!   Calculate total internal partition sum
!   Called by ADJUST, NTECLC, QTNTE.
!   The ratio of total internal partition sums, Q(296K)/Q(Path Temp) is 
!   calculated for a given molecule and isotopic species. The line strength 
!   can then be adjusted by this factor.
!
!   The TIPS values Q(T)/Q(296) is stored in one large block QTEM
!   Array ISPOFF stores the *first* index in QTEM for each molecule/isotope 
!   combination.
!   Array IDGOFF stores the *offset* index in ISPOFF for each molecule.
!
!   The new (2016) Gamache TIPS data is tabulated at 1 or 2K intervals over 
!   temperature ranges 1K to an upper limit which varies between 
!   species/isotope. Since this is huge, the data used here in TPSDAT_DAT
!   represents a subsampled version (20K intervals).
!
!   This uses 4th order Lagrangian interpolation.     
!   In general, Lagrangian interpolation of y(x) is given by
!       y(x) = sum_ijkl  yi(x-xj)(x-xk)(x-xl)/[(xi-xj)(xi-xk)(xi-xl)]
!   where ijkl are the four different points used for interpolation.
!   In this code these points are chosen such that x lies between j and k, 
!   except at the extreme ends of the tabulated temperature where the points 
!   used are the edge points of the tabulation. In this code it it assumed that
!   the xi are equally spaced which allows for significant simplification.
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE TPSDAT_DAT ! Total Internal Partition Sum data
!
! SUBROUTINES
    USE QTWARN_SUB ! Warning messages from QTFCT
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IDGAS ! HITRAN gas ID
    INTEGER(I4), INTENT(IN) :: ISO   ! HITRAN isotope ID 
    REAL(R4),    INTENT(IN) :: TEM   ! Path temperature [K]
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: DTEM3 = DTEM**3 ! Cube of temp.increment
!
! LOCAL VARIABLES
    INTEGER(I4) :: IQTOFF   ! Offset for tabulation if molec/iso data
    INTEGER(I4) :: ISPE     ! Counter for isotope species (all gases)
    INTEGER(I4) :: JT       ! Indices of 4 consecutive tabulation points
    INTEGER(I4) :: NT       ! No of tabulated T points
    REAL(R4)    :: DTI,DTJ,DTK,DTL ! Position of TEM relative to tab. points
    REAL(R4)    :: TEM2     ! Uppper tabulated temperature
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( IDGAS .GE. MAXIDG ) THEN      ! Assume extra molecule
    CALL QTWARN ( IDGAS, -1, TEM )   ! -1 flags this as no TIPS for molec
    QTFCT = ( TEMREF / TEM )**1.5
    RETURN
  END IF
!
  ISPE = IDGOFF(IDGAS) + ISO     ! Position in ISPOFF of molecule/isotope
! Note that IDGOFF array extends to max IDGAS plus 1 to allow for this check
  IF ( ISPE .GT. IDGOFF(IDGAS+1) ) THEN
    CALL QTWARN ( IDGAS, ISO, TEM ) 
    ISPE = IDGOFF(IDGAS) + 1           ! Just use main isotope value instead 
  END IF 
!
  IQTOFF = ISPOFF(ISPE) - 1
  NT = ISPOFF(ISPE+1) - ISPOFF(ISPE)   ! No.tabulated QTEM for mol/iso.
!
  IF ( NT .EQ. 1 ) THEN                ! Special case for monatomic oxygen
    QTFCT = 1.0
    RETURN
  END IF
! Should be a mininum of 4 tabulated temperatures for every other molec/iso.
  IF ( NT .LT. 4 ) STOP 'F-QTFCT: Logical error'
!
! Find indices of 4 consecutive points for interpolation, IT,JT,KT,LT, aiming
! to bracket TEM between JT LE TEM LT KT, but ensuring that the 4 chosen points
! remain within the tabulation if TEM is near the edges
  JT = INT ( ( TEM - TEM1 ) / DTEM ) + 1
  IF ( JT .LT. 2 ) THEN 
    IF ( TEM .LT. TEM1 ) CALL QTWARN ( IDGAS, ISO, TEM, TEM1=TEM1 )
    JT = 2
  ELSE IF ( JT .GT. NT-2 ) THEN
    TEM2 = TEM1 + ( NT - 1 ) * DTEM
    IF ( TEM .GT. TEM2 ) CALL QTWARN ( IDGAS, ISO, TEM, TEM2=TEM2 )
    JT = NT - 2
  END IF
!
  DTI = TEM - ( TEM1 + (JT-2) * DTEM ) 
  DTJ = DTI - DTEM
  DTK = DTJ - DTEM
  DTL = DTK - DTEM
!
  QTFCT = DTEM3 / ( -1.0/6.0*DTJ*DTK*DTL * QTEM(IQTOFF+JT-1)  &
                    +1.0/2.0*DTI*DTK*DTL * QTEM(IQTOFF+JT)  &
                    -1.0/2.0*DTI*DTJ*DTL * QTEM(IQTOFF+JT+1)    &
                    +1.0/6.0*DTI*DTJ*DTK * QTEM(IQTOFF+JT+2)    )
!
END FUNCTION QTFCT
END MODULE QTFCT_FNC
