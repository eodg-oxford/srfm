MODULE SUBAXS_SUB
CONTAINS
SUBROUTINE SUBAXS ( NDP, NDT )
!
! VERSION
!   02JUL24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   Subsample (p,T) axes from LUT TAB input file
!   Called by LUTTAB for each TAB file if qualifiers used to resample p,T axes.
!   Notes: interpretation of resampling increments NDP, NDT
!      +ve sample from first point on axis
!      -ve sample from last point on axis (ie keep last point on axis 
!          rather than first point)
!      0   Take single point nearest middle of axes (changes P1,T1).
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE TABCOM_DAT ! Axes for tabulated absorption coefficients
!
! SUBROUTINES
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE WRTLOG_SUB ! Write text message to log file
!                  
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: NDP ! p-axis subsample factor
    INTEGER(I4), INTENT(IN) :: NDT ! T-axis subsample factor
!
! LOCAL VARIABLES
    INTEGER(I4) :: IP       ! DO loop counter for original -lnp axis
    INTEGER(I4) :: IP1, IP2 ! DO loop limits for subsampling p-axis
    INTEGER(I4) :: IPO      ! Original p-axis index
    INTEGER(I4) :: IPSPAN   ! Range of Npts of orig p-axis that will be sampled
    INTEGER(I4) :: IPSTEP   ! Step size for subsampling p-axis
    INTEGER(I4) :: IQ       ! Counter for q-axis points
    INTEGER(I4) :: IT       ! DO loop counter for original Temp axis
    INTEGER(I4) :: IT1, IT2 ! DO loop limits for subsampling T-axis
    INTEGER(I4) :: ITO      ! Original T-axis index
    INTEGER(I4) :: ITSPAN   ! Range of Npts of Orig T-axis that will be sampled
    INTEGER(I4) :: ITSTEP   ! Step size for subsampling T-axis
    INTEGER(I4) :: IX       ! Counter for original axis
    INTEGER(I4) :: NPORG    ! Original no. of p-axis points
    INTEGER(I4) :: NTORG    ! Original no. of T-axis points
    INTEGER(I4), ALLOCATABLE :: IPORG(:) ! Original p-axis indices of new pts
    INTEGER(I4), ALLOCATABLE :: ITORG(:) ! Original T-axis indices of new pts
!
! EXECUTABLE CODE --------------------------------------------------------------
!
   CALL WRTLOG ( 'I-LUTTAB: Subsampling p,T axes by factors ' // &
                  TRIM ( C11INT(NDP) ) // TRIM ( C11INT(NDT) ) )
!
  ALLOCATE ( IPORG(NPTAB) )
  ALLOCATE ( ITORG(NTTAB) ) 
!
  IF ( NDP .EQ. 0 ) THEN       
    IP1 = 1 + INT ( ( NPTAB - 1 ) / 2 ) ! NP=1,2,3,4,5 :IP1=1,1,2,2,3
    IP2 = IP1
    IPSTEP = 1
  ELSE
    IPSTEP = ABS ( NDP ) 
    IPSPAN = ( ( NPTAB - 1 ) / IPSTEP ) * IPSTEP
    IF ( NDP .GT. 0 ) THEN
      IP1 = 1
    ELSE
      IP1 = NPTAB - IPSPAN
    END IF
    IP2 = IP1 + IPSPAN
  END IF
  IX = 0
  DO IP = IP1, IP2, IPSTEP
    IX = IX + 1
    PAXTAB(IX) = PAXTAB(IP)
    TEMTAB(IX) = TEMTAB(IP)
    VMRTAB(IX,1) = VMRTAB(IP,1)
    IPORG(IX) = IP
  END DO
  NPORG = NPTAB
  NPTAB = IX
!
  IF ( NDT .EQ. 0 ) THEN       
    IT1 = 1 + INT ( ( NTTAB - 1 ) / 2 ) 
    IT2 = IT1
    ITSTEP = 1
  ELSE
    ITSTEP = ABS ( NDT ) 
    ITSPAN = ( ( NTTAB - 1 ) / ITSTEP ) * ITSTEP
    IF ( NDT .GT. 0 ) THEN
      IT1 = 1
    ELSE
      IT1 = NTTAB - ITSPAN
    END IF
    IT2 = IT1 + ITSPAN
  END IF
  IX = 0
  DO IT = IT1, IT2, ITSTEP
    IX = IX + 1
    TAXTAB(IX) = TAXTAB(IT)
    ITORG(IX) = IT
  END DO
  NTORG = NTTAB
  NTTAB = IX
!
  IX = 0
  ALLOCATE ( IXORG(NPTAB*NTTAB*NQTAB) )   ! IXORG in TABCOM
  DO IP = 1, NPTAB
    IPO = IPORG(IP)
    DO IT = 1, NTTAB
      ITO = ITORG(IT) 
      DO IQ = 1, NQTAB
        IX = IX + 1
        IXORG(IX) = (IPO-1) * NTORG * NQTAB + (ITO-1) * NQTAB + IQ
      END DO
    END DO
  END DO 
!
END SUBROUTINE SUBAXS
END MODULE SUBAXS_SUB
