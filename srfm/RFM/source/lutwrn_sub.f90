MODULE LUTWRN_SUB
CONTAINS
SUBROUTINE LUTWRN ( IGAS, NTOT, NP1LIM, NP2LIM, NT1LIM, NT2LIM, NQ1LIM, NQ2LIM )
!
! VERSION
!   06JUN24 AD Checked.
!   01MAY17 AD F90 conversion of part of lutpth.for. Checked.
! 
! DESCRIPTION    
!   Warn if TAB-LUT p,T,q axis limits exceeded
!   Called by LUTWGT
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE TABCOM_DAT ! Axes for tabulated absorption coefficients
!
! SUBROUTINES
    USE C11INT_FNC ! Write integer as left-adjusted string
    USE C9REAL_GEN ! Write real number as C*9 string
    USE NAMGAS_FNC ! Return molecule name + (iso) associated with GASCOM index
    USE WRTLOG_SUB ! Write text message to log file
!                  
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IGAS   ! Index of absorber 
    INTEGER(I4), INTENT(IN) :: NTOT   ! Total number of calculated paths
    INTEGER(I4), INTENT(IN) :: NP1LIM ! No.times p outside lower index
    INTEGER(I4), INTENT(IN) :: NP2LIM ! No.times p outside upper index
    INTEGER(I4), INTENT(IN) :: NT1LIM ! No.times T outside lower index
    INTEGER(I4), INTENT(IN) :: NT2LIM ! No.times T outside upper index
    INTEGER(I4), INTENT(IN) :: NQ1LIM ! No.times q outside lower index
    INTEGER(I4), INTENT(IN) :: NQ2LIM ! No.times q outside upper index
!
! LOCAL VARIABLES
    INTEGER(I4)   :: NPMAX  ! No. times axis max pressure exceeded
    INTEGER(I4)   :: NQMAX  ! No. times axis max VMR factor exceeded
    INTEGER(I4)   :: NTMAX  ! No. times axis max temperature exceeded
    INTEGER(I4)   :: NTMIN  ! No. times axis min temperature exceeded
    CHARACTER(9)  :: PCTSTR ! Percentage written as a string
    CHARACTER(80) :: WRNMSG ! Warning message
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  IF ( NPTAB .GT. 1 ) THEN
    IF ( PAXTAB(NPTAB) .GT. PAXTAB(1) ) THEN
      NPMAX = NP2LIM
    ELSE 
      NPMAX = NP1LIM
    END IF
    IF ( NPMAX .GT. 0 ) THEN
      PCTSTR = C9REAL ( 100.0 * FLOAT ( NPMAX ) / FLOAT ( NTOT ) ) 
      WRNMSG = 'W-LUTWRN: ' // TRIM ( NAMGAS(IGAS) ) // &
               ' Path Pre > LUT Max Pre for ' // TRIM ( C11INT(NPMAX) ) // &
               ' calc paths (' // TRIM ( PCTSTR ) // '%)'
      CALL WRTLOG ( WRNMSG ) 
    END IF
  END IF
!
  IF ( NTTAB .GT. 1 ) THEN
    IF ( TAXTAB(NTTAB) .GT. TAXTAB(1) ) THEN
      NTMAX = NT2LIM
      NTMIN = NT1LIM
    ELSE
      NTMAX = NT1LIM
      NTMIN = NT2LIM
    END IF
    IF ( NTMAX .GT. 0 ) THEN
      PCTSTR = C9REAL ( 100.0 * FLOAT ( NTMAX ) / FLOAT ( NTOT ) ) 
      WRNMSG = 'W-LUTWRN: ' // TRIM ( NAMGAS(IGAS) ) // &
               ' Path Tem > LUT Max Tem for ' // TRIM ( C11INT(NTMAX) ) // &
               ' calc paths (' // TRIM ( PCTSTR ) // '%)'
      CALL WRTLOG ( WRNMSG )
    END IF
    IF ( NTMIN .GT. 0 ) THEN
      PCTSTR = C9REAL ( 100.0 * FLOAT ( NTMIN ) / FLOAT ( NTOT ) )
      WRNMSG = 'W-LUTWRN: ' // TRIM ( NAMGAS(IGAS) ) // &
               ' Path Tem < LUT Min Tem for ' // TRIM ( C11INT(NTMIN) ) // &
               ' calc paths (' // TRIM ( PCTSTR ) // '%)'
      CALL WRTLOG ( WRNMSG )
    END IF
  END IF
!
  IF ( NQTAB .GT. 1 ) THEN
    IF ( QAXTAB(NQTAB) .GT. QAXTAB(1) ) THEN
      NQMAX = NQ2LIM
    ELSE 
      NQMAX = NQ1LIM
    END IF
    IF ( NQMAX .GT. 0 ) THEN
      PCTSTR = C9REAL ( 100.0 * FLOAT ( NQMAX ) / FLOAT ( NTOT ) ) 
      WRNMSG = 'W-LUTWRN: ' // TRIM ( NAMGAS(IGAS) ) // &
               ' Path VMR > LUT Max VMR for ' // TRIM ( C11INT(NQMAX) ) // &
               ' calc paths (' // TRIM ( PCTSTR ) // '%)'
      CALL WRTLOG ( WRNMSG ) 
    END IF
  END IF
!
END SUBROUTINE LUTWRN
END MODULE LUTWRN_SUB

