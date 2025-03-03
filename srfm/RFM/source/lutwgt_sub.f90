MODULE LUTWGT_SUB
CONTAINS
SUBROUTINE LUTWGT ( FINISH, LINVMR, IGAS, PRE, TEM, PPA, IDXLKP, WGTLKP )
!
! VERSION
!   16JUN24 AD Checked.
!   01JUL17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   Calculate LUT p,T,q-axis interpolation factors for path
!   Called by LUTPTH for each LUT file.
!
! VARIABLE KINDS
    USE KIND_DAT
!                  
! GLOBAL DATA
    USE TABCOM_DAT ! Axes for tabulated absorption coefficients
    USE PHYCON_DAT, ONLY: ATMB ! Std Atmos. Pressure [mb]
!
! SUBROUTINES
    USE IBRAKT_GEN ! Lower index of array interpolation
    USE LUTWRN_SUB ! Warn if TAB-LUT p,T,q axis limits exceeded
!
  IMPLICIT NONE
!
! ARGUMENTS
    LOGICAL,     INTENT(IN)  :: FINISH    ! T=last path, F=continue
    LOGICAL,     INTENT(IN)  :: LINVMR    ! T=linear interp, F=log interp VMR
    INTEGER(I4), INTENT(IN)  :: IGAS      ! Index of absorber
    REAL(R4),    INTENT(IN)  :: PRE       ! Pressure [mb]
    REAL(R4),    INTENT(IN)  :: TEM       ! Temperature [K]
    REAL(R4),    INTENT(IN)  :: PPA       ! Partial pressure [mb]
    INTEGER(I4), INTENT(OUT) :: IDXLKP(:) ! Interpolation points 
    REAL(I4),    INTENT(OUT) :: WGTLKP(:) ! Interpolation weights
!
! LOCAL VARIABLES
    INTEGER(I4) :: IP = 1     ! Index of lower point in p axis
    INTEGER(I4) :: IQ = 1     ! Index of lower point in q axis
    INTEGER(I4) :: IT = 1     ! Index of lower point in T axis
    INTEGER(I4) :: JP,JQ,JT   ! Indices of upper points in p,q,T-axes
    INTEGER(I4) :: KP,KQ,KT   ! Counter for lower/upper points in p,q,T-axes
    INTEGER(I4) :: ILKP       ! Counter for interpolation points
    INTEGER(I4) :: NP1LIM = 0 ! Counter: no.times p outside PAXTAB lower index
    INTEGER(I4) :: NP2LIM = 0 ! Counter: no.times p outside PAXTAB upper index
    INTEGER(I4) :: NQ1LIM = 0 ! Counter: no.times q outside QAXTAB lower index
    INTEGER(I4) :: NQ2LIM = 0 ! Counter: no.times q outside QAXTAB upper index
    INTEGER(I4) :: NT1LIM = 0 ! Counter: no.times T outside TAXTAB lower index
    INTEGER(I4) :: NT2LIM = 0 ! Counter: no.times T outside TAXTAB upper index
    INTEGER(I4) :: NTOT   = 0 ! Total number of calculated paths for IGAS
    INTEGER(I4) :: NTQ        ! No.of T*q points for each p-axis value
    REAL(R4)    :: DP,DQ,DT   ! Fraction of p,q,T-axis interval above IP,IQ,IT
    REAL(R4)    :: PAX        ! Pressure [mb] of path segment 
    REAL(R4)    :: PFACT      ! Interpolation weight for p-axis (DP or 1-DP)
    REAL(R4)    :: QAX        ! % of tab.VMR of path segment 
    REAL(R4)    :: QFACT      ! Interpolation weight for q-axis (DQ or 1-DQ)
    REAL(R4)    :: TAX        ! (Relative) Temperature [K] of path segment
    REAL(R4)    :: TFACT      ! Interpolation weight for T-axis (DT or 1-DT)
    REAL(R4)    :: VMR        ! Volume mixing ratio [ppmv] of path segment
    REAL(R4)    :: XFACT      ! Product of interpolation weights P,Q,TFACT
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  NTOT = NTOT + 1
!
  PAX = PRE * ATMB                ! convert [atm] to [mb]
! Arrays PAXTAB, TAXTAB etc may have been sub-setted by SUBAXS so need to be
! explicit with number of elements currently filled
  IP = IBRAKT ( PAXTAB(1:NPTAB), PAX, IP ) 
  IF ( IP .LT. 1 ) THEN
    IP = 1
    DP = 0.0
    NP1LIM = NP1LIM + 1
  ELSE IF ( IP .EQ. NPTAB ) THEN
    IP = NPTAB - 1
    DP = 1.0
    NP2LIM = NP2LIM + 1
  ELSE
    DP = LOG ( PAX/PAXTAB(IP) ) / LOG (PAXTAB(IP+1)/PAXTAB(IP)) 
  END IF
! Note that number of interpolation points is fixed for each LUT, so need to
! ensure that all points are filled, eg JP=IP+1 unless NPTAB=1
  JP = MIN ( NPTAB, IP + 1 )     
!
  TAX = TEM
  IF ( OFFTAB ) TAX = TAX - ( TEMTAB(IP) * (1.0-DP) + TEMTAB(JP) * DP ) 
  IT = IBRAKT ( TAXTAB(1:NTTAB), TAX, IT ) 
  IF ( IT .LT. 1 ) THEN
    IT = 1
    DT = 0.0
    NT1LIM = NT1LIM + 1
  ELSE IF ( IT .EQ. NTTAB ) THEN
    IT = NTTAB - 1
    DT = 1.0
    NT2LIM = NT2LIM + 1
  ELSE
    DT = ( TAX - TAXTAB(IT) ) / ( TAXTAB(IT+1) - TAXTAB(IT) )
  END IF
  JT = MIN ( NTTAB, IT + 1 ) 
!           
  VMR = PPA / PRE * 1.0E6
  IF ( LINVMR ) THEN
    QAX = 100.0 * VMR / ( VMRTAB(IP,1) * (1.0-DP) + VMRTAB(JP,1) * DP ) 
  ELSE
    QAX = 100.0 * VMR / &
          EXP ( LOG ( VMRTAB(IP,1) ) * (1.0-DP) + LOG ( VMRTAB(JP,1) ) * DP ) 
  END IF  
  IQ = IBRAKT ( QAXTAB(1:NQTAB), QAX, IQ )
  IF ( IQ .LT. 1 ) THEN
    IQ = 1
    DQ = 0.0
    NQ1LIM = NQ1LIM + 1
  ELSE IF ( IQ .EQ. NQTAB ) THEN
    IQ = NQTAB - 1
    DQ = 1.0
    NQ2LIM = NQ2LIM + 1
  ELSE
    IF ( LINVMR ) THEN 
      DQ = ( QAX - QAXTAB(IQ) ) / ( QAXTAB(IQ+1) - QAXTAB(IQ) ) 
    ELSE
      DQ = LOG (QAX/QAXTAB(IQ)) / LOG (QAXTAB(IQ+1)/QAXTAB(IQ)) 
    END IF
  END IF
  JQ = MIN ( NQTAB, IQ + 1 ) 
!
! Calculate weights for lower and upper points in p,T,q axes
  ILKP = 0
  NTQ = NTTAB * NQTAB
  PFACT = 1.0 - DP
  DO KP = IP, JP
    TFACT = 1.0 - DT
    DO KT = IT, JT
      QFACT = 1.0 - DQ
      DO KQ = IQ, JQ
        XFACT = PFACT * TFACT * QFACT 
        ILKP = ILKP + 1
        WGTLKP(ILKP) = XFACT
        IDXLKP(ILKP) = (KP-1)*NTQ + (KT-1)*NQTAB + KQ
        QFACT = DQ
      END DO
      TFACT = DT
    END DO
    PFACT = DP
  END DO                                
!
  IF ( ALLOCATED(IXORG) ) IDXLKP = IXORG(IDXLKP) 
!
  IF ( .NOT. FINISH ) RETURN
!
  CALL LUTWRN ( IGAS, NTOT, NP1LIM, NP2LIM, NT1LIM, NT2LIM, NQ1LIM, NQ2LIM ) 

! Reset counters
  NP1LIM = 0
  NP2LIM = 0
  NT1LIM = 0
  NT2LIM = 0
  NQ1LIM = 0
  NQ2LIM = 0
  NTOT = 0
!
! Ensure initial guess values are defined for IBRAKT
  IP = 1
  IT = 1
  IQ = 1
!
END SUBROUTINE LUTWGT
END MODULE LUTWGT_SUB
