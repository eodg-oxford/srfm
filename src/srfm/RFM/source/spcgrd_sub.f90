MODULE SPCGRD_SUB
CONTAINS
SUBROUTINE SPCGRD ( IWID, NOGRID )
!
! VERSION
!   14APR24 AD Checked.
!   03MAR19 AD Bug#17: Add NOGRID argument. Checked.
!   01MAY17 AD F90 conversion of F77 module RFMGRD. Checked.
!
! DESCRIPTION    
!   Sets up fine-resolution grid for each widemesh interval
!   Called by RFMSPC for each widemesh interval.
!   If the GRD option is enabled, sets NFIN and wavenumber grid points
!   If not enabled, just sets uniform grid.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FINCOM_DAT ! Finemesh data
    USE FULCOM_DAT ! Full grid data
    USE WIDCOM_DAT ! Widemesh data
    USE CLCCOM_DAT, ONLY: NCLC ! No. calculated path segments
    USE RFMCON_DAT, ONLY: FEXC ! Window for fine mesh calcs
!
! SUBROUTINES 
    USE IBRAKT_GEN ! Lower index of array interpolation
!
  IMPLICIT NONE
!
! ARGUMENTS      
    INTEGER(I4), INTENT(IN)  :: IWID   ! Index of current wide mesh intvl [IPW]
    LOGICAL,     INTENT(OUT) :: NOGRID ! T=No grid pts within interval
!
! LOCAL VARIABLES
    REAL(R8)    :: DW       ! Negligible wno [/cm] to allow for rounding
    INTEGER(I4) :: LFIN = 0 ! Last value of NFIN
!
! EXECUTABLE CODE -------------------------------------------------------------
! 
  DW    = WNRFUL * 0.1D0
  IFUL1 = IFUL2 + 1
  IFUL2 = IBRAKT ( WNOFUL, WNOWID(IWID) - DW, IFUL2 )
! If there are no irreg grid pts in current widemesh interval IFUL2 will be
! returned as the same highest grid pt as was used in previous interval
  NOGRID = IFUL2 .LT. IFUL1 
  IF ( NOGRID ) RETURN
!
  IF1FIN = IFUL1
  IF2FIN = IFUL2
  NFIN = IFUL2 - IFUL1 + 1
!
  IF ( NFIN .NE. LFIN ) THEN  
    IF ( ALLOCATED(WNOFIN) ) DEALLOCATE ( WNOFIN, ABSFIN )
    ALLOCATE ( WNOFIN(NFIN), ABSFIN(NFIN,NCLC) ) 
    IF ( USECNT ) THEN
      IF ( ALLOCATED ( CNTFIN ) ) DEALLOCATE ( CNTFIN ) 
      ALLOCATE ( CNTFIN(NFIN,NCLC) )
    END IF
  END IF
  LFIN = NFIN
  ABSFIN = 0.0
!
  WNOFIN = WNOFUL(IFUL1:IFUL2)
  WN1FIN = WNOFIN(1)
  WN2FIN = WNOFIN(NFIN)
  WNLFIN = WNOWID(IWID-1) - FEXC
  WNUFIN = WNOWID(IWID) + FEXC
!
END SUBROUTINE SPCGRD
END MODULE SPCGRD_SUB
