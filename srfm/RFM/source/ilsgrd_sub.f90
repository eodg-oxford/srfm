MODULE ILSGRD_SUB
CONTAINS
SUBROUTINE ILSGRD ( IFNOUT, IP1ILS, IP2ILS, FNCILS, IGRD1, IGRD2, FNCIRR ) 
!
! VERSION
!   19JAN23 AD Bug#37: Revise Logical Error#5. Checked.
!   19APR22 AD Bug#34: Allow for all contrib.grid points above irreg grid.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Calculate ILS fn for irregular grid
!   Called by SPCILS if ILS or AVG flags enabled with irregular grid
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE GRDCOM_DAT ! Irregular grid
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN)    :: IFNOUT    ! Output pt on fine grid
    INTEGER(I4), INTENT(IN)    :: IP1ILS    ! Lower bound of ILS function
    INTEGER(I4), INTENT(IN)    :: IP2ILS    ! Upper bound of ILS function
    REAL(R8),    INTENT(IN)    :: FNCILS(IP1ILS:IP2ILS) ! ILS function
    INTEGER(I4), INTENT(INOUT) :: IGRD1     ! Lower irreg grid pt for convol.
    INTEGER(I4), INTENT(OUT)   :: IGRD2     ! Upper irreg grid pt for convol. 
    REAL(R8), ALLOCATABLE, &
                 INTENT(OUT)   :: FNCIRR(:) ! Convolution weights for irr.grid
!
! LOCAL VARIABLES
    INTEGER(I4) :: IFIN    ! Counter for fine grid points
    INTEGER(I4) :: IFINGL  ! Fine grid index of lower irreg grid point
    INTEGER(I4) :: IFINGU  ! Fine grid index of upper irreg grid point
    INTEGER(I4) :: IFINL   ! Index of lower fine grid point
    INTEGER(I4) :: IFINU   ! Index of upper fine grid point
    INTEGER(I4) :: IGRD    ! Counter for irreg.grid points
    INTEGER(I4) :: IPILS   ! Counter for ILS function fine grid points
    REAL(R8)    :: A       ! Interpolation weight 
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Establish range of full ILS convolution on fine grid IFINL:IFINU
  IFINL = IFNOUT + IP1ILS
  IFINU = IFNOUT + IP2ILS
!
! In case this subroutine is called with no overlap between ILS and
! irregular grid, return equivalent of zero contribution 
  IF ( IFINU .LE. GRD(1)%IFN .OR. IFINL .GE. GRD(NGRD)%IFN ) THEN
    IGRD1 = 1
    IGRD2 = 1
    IF ( ALLOCATED ( FNCIRR ) ) DEALLOCATE ( FNCIRR ) 
    ALLOCATE ( FNCIRR(1) )
    FNCIRR = 0.0D0
    RETURN
  END IF   
!
! Also possible that called with ILS overlapping end of grid, so adjust
  IFINL = MAX ( IFINL, GRD(1)%IFN )
  IFINU = MIN ( IFINU, GRD(NGRD)%IFN )
  IF ( IFINU .LE. IFINL ) STOP 'F-ILSGRD: Logical error#1'
!
! Ensure that IGRD1 can be accessed
  IF ( IGRD1 .LT. 1 .OR. IGRD1 .GE. NGRD ) IGRD1 = 1 
!
! Establish lower irreg.grid point IGRD1 .le. IFINL.
! Assume called with increasing output grid points so if IGRD1 above IFINL reset
  IF ( GRD(IGRD1)%IFN .GT. IFINL ) IGRD1 = 1     
! Increase IGRD1 until next grid point lies above IFINL (so IGRD1 .le. IFINL)
  DO WHILE ( GRD(IGRD1+1)%IFN .LE. IFINL ) 
    IGRD1 = IGRD1 + 1
    IF ( IGRD1 .EQ. NGRD ) STOP 'F-ILSGRD: Logical error#2'
  END DO 
! 
! Establish upper irreg.grid point IGRD2 .ge. IFINU
! Since IFINU limited to GRD(NGRD)%IFN IGRD2 is limited to NGRD as upper limit
  IF ( IGRD1 .GE. NGRD ) STOP 'F-ILSGRD: Logical error#3'
  IGRD2 = IGRD1 + 1
  DO WHILE ( GRD(IGRD2)%IFN .LT. IFINU )
    IGRD2 = IGRD2 + 1                     ! end with IGRD2 .ge. IFINU
    IF ( IGRD2 .GT. NGRD ) STOP 'F-ILSGRD: Logical error#4'
  END DO
!
  IF ( ALLOCATED ( FNCIRR ) ) DEALLOCATE ( FNCIRR ) 
  ALLOCATE ( FNCIRR(IGRD1:IGRD2) ) 
  FNCIRR = 0.0D0
!
! Indices on full fine grid of lower,upper irreg grid points
  IGRD = IGRD1
  IFINGL = GRD(IGRD)%IFN
  IFINGU = GRD(IGRD+1)%IFN
  DO IFIN = IFINL, IFINU    ! Loop over ILS width on fine grid
    IPILS = IFIN - IFNOUT   ! Index of pt within ILS fn.
    A = DBLE ( ( IFIN - IFINGL ) / ( IFINGU - IFINGL ) )
    FNCIRR(IGRD)   = FNCIRR(IGRD)   + (1.0-A) * FNCILS(IPILS)
    FNCIRR(IGRD+1) = FNCIRR(IGRD+1) + A * FNCILS(IPILS)
    IF ( IFIN .EQ. IFINGU ) THEN ! next irreg grid pt
      IGRD = IGRD + 1
      IFINGL = IFINGU
      IFINGU = GRD(IGRD+1)%IFN
      IF ( IGRD .EQ. NGRD .AND. IFIN .NE. IFINU ) &
        STOP 'F-ILSGRD: Logical error#5'
    END IF
  END DO
!
END SUBROUTINE ILSGRD
END MODULE ILSGRD_SUB

