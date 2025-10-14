MODULE YMIX_FNC
CONTAINS
REAL(R4) FUNCTION YMIX ( TEM, PRE, PPA )
!
! VERSION
!   08SEP24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
! 
! DESCRIPTION
!   Calculate line-mixing y-coefficient
!   Called by ADJUST if MIXFLG enabled.
!                      CO2 Q-branch mixing coefficient data for:
!                  1   618 cm-1  3  2, 626 (Q50 - Q 2,  615 - 619 cm-1)
!                  2   648 cm-1  2  1, 636 (Q 2 - Q50,  648 - 652 cm-1)
! (currently dummy)3   662 cm-1  2  1, 628 (Q 1 - Q50,  662 - 665 cm-1)
!                  4   667 cm-1  2  1, 626 (Q 2 - Q50,  667 - 671 cm-1)
!                  5   668 cm-1  4  2, 626 (Q 2 - Q81,  667 - 675 cm-1)
!                  6   720 cm-1  5  2, 626 (Q50 - Q 2,  718 - 721 cm-1)
!                  7   741 cm-1  8  4, 626 (Q81 - Q 2,  733 - 742 cm-1)
!                  8   791 cm-1  8  3, 626 (Q 2 - Q50,  791 - 794 cm-1)
!                  9  1932 cm-1  6  1, 626 (Q 2 - Q50, 1932 -1937 cm-1)
!                 10  2076 cm-1  8  1, 626 (Q 2 - Q50, 2076 -2080 cm-1)
!                 11  2093 cm-1 14  2, 626 (Q 2 - Q70, 2093 -2095 cm-1)
!                 12  2129 cm-1 15  2, 626 (Q50 - Q 2, 2128 -2130 cm-1)
!   Refer to Tobin & Strow for details of quantum numbers.
!   To be used with CO2MIX lines with LABEL '101293'
!
!   1st order line mixing = Y(T)*(200/T)^0.75, with
!   Y(T) = COMIX(1) + COMIX(2)*(T - 200) +
!   COMIX(3)*(T - 200)**2 + COMIX(4)*(T - 200)**3
!
!   Data relating to the ranges of sets is incorporated in WIDMIX
!
! REFERENCES
!   TOBIN, D. and L.L. STROW
!   A Compilation of First-order Line-mixing Coefficients for CO2 Q-branches
!   J.Q.S.R.T., 52, 281 (1994).
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE HITCOM_DAT ! HITRAN line data
    USE MIXDAT_DAT ! CO2 line-mixing data
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R4), INTENT(IN) :: TEM ! CG temperature [K]
    REAL(R4), INTENT(IN) :: PRE ! CG Pressure [atm]
    REAL(R4), INTENT(IN) :: PPA ! CG Partial Pressure [atm]
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: TEMREF = 200.0 ! Reference temp. for calculations
!  
! LOCAL VARIABLES
    INTEGER(I4) :: ILIN   ! Counter for lines in each set [IQ]
    INTEGER(I4) :: ISET   ! Counter for number of sets [IS]
    REAL(R8)    :: TEMDIF ! Difference from calc.ref. temperature
    REAL(R8)    :: YMIXF  ! Foreign-broadened mixing
    REAL(R8)    :: YMIXS  ! Self-broadened mixing
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO ISET = 1, NSET
    IF ( IMMOL(ISET) .EQ. HIT%IDM .AND. & 
         IMISO(ISET) .EQ. HIT%IDI .AND. &
          MGQU(ISET) .EQ. HIT%IUS .AND. &
          MGQL(ISET) .EQ. HIT%ILS        ) THEN
      DO ILIN = 1, NMLN(ISET)
        IF ( HIT%BLQ(5:9) .EQ. MLQL(ILIN,ISET)(5:9) ) THEN
          TEMDIF = DBLE ( TEM - TEMREF )
          YMIXS =   COMIXS(1,ILIN,ISET) + TEMDIF * &
                  ( COMIXS(2,ILIN,ISET) + TEMDIF * &
                  ( COMIXS(3,ILIN,ISET) + TEMDIF * &
                    COMIXS(4,ILIN,ISET)            )) 
          YMIXF =   COMIXF(1,ILIN,ISET) + TEMDIF * &
                  ( COMIXF(2,ILIN,ISET) + TEMDIF * &
                  ( COMIXF(3,ILIN,ISET) + TEMDIF * &
                    COMIXF(4,ILIN,ISET)            ))
          YMIX = (TEMREF/TEM)**0.75 * &
                 SNGL ( DBLE(PPA)*YMIXS + DBLE(PRE-PPA)*YMIXF )
          RETURN
        END IF
      END DO
      EXIT
    END IF 
  END DO
!
  YMIX = 0.0
!
END FUNCTION YMIX
END MODULE YMIX_FNC
