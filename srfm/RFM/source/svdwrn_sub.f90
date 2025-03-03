MODULE SVDWRN_SUB
CONTAINS
SUBROUTINE SVDWRN ( IGAS, NP, DP, NT, DT, &
                    NTOT, NP1LIM, NP2LIM, NT1LIM, NT2LIM )
! VERSION
!   01JUL24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   Warn if SVD-LUT p,T axis limits exceeded
!   Called by SVDWGT
!
! VARIABLE KINDS
    USE KIND_DAT
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
    INTEGER(I4), INTENT(IN) :: NP     ! No. -lnp axis points
    REAL(R4),    INTENT(IN) :: DP     ! -lnp axis interval
    INTEGER(I4), INTENT(IN) :: NT     ! No. Tem axis points
    REAL(R4),    INTENT(IN) :: DT     ! Tem axis interval
    INTEGER(I4), INTENT(IN) :: NTOT   ! Total number of calculated paths
    INTEGER(I4), INTENT(IN) :: NP1LIM ! No.times p outside lower index
    INTEGER(I4), INTENT(IN) :: NP2LIM ! No.times p outside upper index
    INTEGER(I4), INTENT(IN) :: NT1LIM ! No.times T outside lower index
    INTEGER(I4), INTENT(IN) :: NT2LIM ! No.times T outside upper index
!
! LOCAL VARIABLES
    INTEGER(I4)   :: NPMAX  ! No. times axis max pressure exceeded
    INTEGER(I4)   :: NTMAX  ! No. times axis max temperature exceeded
    INTEGER(I4)   :: NTMIN  ! No. times axis min temperature exceeded
    CHARACTER(9)  :: PCTSTR ! Percentage written as a string
    CHARACTER(80) :: WRNMSG ! Warning message
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  IF ( NP .GT. 1 ) THEN
    IF ( DP .GT. 0.0 ) THEN 
      NPMAX = NP1LIM                    ! DP +ve is max(p) as P1
    ELSE 
      NPMAX = NP2LIM
    END IF
    IF ( NPMAX .GT. 0 ) THEN
      PCTSTR = C9REAL ( 100.0 * FLOAT ( NPMAX ) / FLOAT ( NTOT ) ) 
      WRNMSG = 'W-SVDWRN: ' // TRIM ( NAMGAS(IGAS) ) // &
               ' Path Pre > LUT Max Pre for' // TRIM ( C11INT(NPMAX) ) // &
               ' calc paths (' // TRIM ( PCTSTR ) // '%)'
      CALL WRTLOG ( WRNMSG ) 
    END IF
  END IF
!
  IF ( NT .GT. 1 ) THEN
    IF ( DT .LT. 0.0 ) THEN  ! DT -ve is max(T) as T1
      NTMAX = NT1LIM
      NTMIN = NT2LIM
    ELSE
      NTMAX = NT2LIM
      NTMIN = NT1LIM
    END IF
    IF ( NTMAX .GT. 0 ) THEN
      PCTSTR = C9REAL ( 100.0 * FLOAT ( NTMAX ) / FLOAT ( NTOT ) ) 
      WRNMSG = 'W-SVDWRN: ' // TRIM ( NAMGAS(IGAS) ) // &
               ' Path Tem > LUT Max Tem for' // TRIM ( C11INT(NTMAX) ) // &
               ' calc paths (' // TRIM ( PCTSTR ) // '%)'
      CALL WRTLOG ( WRNMSG )
    END IF
    IF ( NTMIN .GT. 0 ) THEN
      PCTSTR = C9REAL ( 100.0 * FLOAT ( NTMIN ) / FLOAT ( NTOT ) )
      WRNMSG = 'W-SVDWRN: ' // TRIM ( NAMGAS(IGAS) ) // &
               ' Path Tem < LUT Min Tem for' // TRIM ( C11INT(NTMIN) ) // &
               ' calc paths (' // TRIM ( PCTSTR ) // '%)'
      CALL WRTLOG ( WRNMSG )
    END IF
  END IF
!
END SUBROUTINE SVDWRN
END MODULE SVDWRN_SUB

