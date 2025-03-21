MODULE LSTABS_SUB
CONTAINS
SUBROUTINE LSTABS ( WNOL, WNOH, IMIN, NREQ, MOLLST, OPTLST )
!
! VERSION
!   09JUL24 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Construct ordered list of absorbers for given spec.range
!   Called by GASALL
!   Molecules are returned in order of decreasing absorption
!   (largest value of OPTLST is first element, corresponding to maximum
!   absorption).      
!
! VARIABLE KINDS
    USE KIND_DAT 
!
! GLOBAL DATA
    USE OPTDAT_DAT ! Optical depth data for each molecule
!
  IMPLICIT NONE
!
! ARGUMENTS
    REAL(R8),    INTENT(IN)  :: WNOL ! Lower Wno [cm-1] of required range
    REAL(R8),    INTENT(IN)  :: WNOH ! Upper Wno [cm-1] of required range
    INTEGER(I4), INTENT(IN)  :: IMIN ! Threshold for min. optical depth
    INTEGER(I4), INTENT(OUT) :: NREQ ! No. molecules found
    INTEGER(I4), ALLOCATABLE, INTENT(OUT) :: MOLLST(:) ! List of molec.indices
    REAL(R4),    ALLOCATABLE, INTENT(OUT) :: OPTLST(:) ! List of opt.depths
!
! LOCAL VARIABLES
    INTEGER(I4) :: IDX    ! Offset index for molecule in arrays OLOG10,WAVENO
    INTEGER(I4) :: ILOG   ! Current value of OLOG10 array
    INTEGER(I4) :: IMAX   ! Maximum value of ILOG for current molecule
    INTEGER(I4) :: IMOL   ! Counter for all HITRAN species (1:MAXMOL)
    INTEGER(I4) :: IREQ   ! Counter for reqd absorbers in decreasing order
    INTEGER(I4) :: ISUM   ! Summation of opt.dep.values for each molecule
    INTEGER(I4) :: IWNO   ! Lower Wno of each molecules opt.dep. data range
    INTEGER(I4) :: IWNOH  ! Integer wavenumber corresponding to WNOH
    INTEGER(I4) :: IWNOL  ! Integer wavenumber corresponding to WNOL 
    INTEGER(I4) :: JWNO   ! Upper Wno of each molecules opt.dep data range
    INTEGER(I4) :: NGAS   ! No.absorbers found (output value of NREQ)
    REAL(R4)    :: FACTOR ! Factor for 'averaging' optical depth
    REAL(R4)    :: OPTMOL(MAXMOL) ! Combined maximum & average optical.depths
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  OPTMOL = 0
  NGAS  = 0                          ! Initialise no.of absorbers found
  IWNOL = INT ( WNOL )
  IWNOH = INT ( WNOH )
  FACTOR = 0.01 / FLOAT ( IWNOH - IWNOL + 1 ) 
  DO IMOL = 1, MAXMOL                 ! MAXMOL in OPTDAT
    IF ( IDXMOL(IMOL) .EQ. 0 ) CYCLE  ! Exclude undefined molecular indices
    IDX = IDXMOL(IMOL)
    JWNO = WAVENO(IDX)             ! Upper limit of first range
    DO WHILE ( JWNO .LT. IWNOL )   ! Step over data ranges until overlap
      IDX = IDX + 1
      JWNO = JWNO + WAVENO(IDX)    ! Upper limit of next data range
    END DO    
    IWNO = JWNO - WAVENO(IDX) + 1  ! Lower wno of current data range   
    ILOG = OLOG10(IDX)             ! Value of Opt.Dep at WNOL
    IMAX = ILOG                    ! Initialise IMAX for this molecule
    ISUM = ( MIN(JWNO,IWNOH) - IWNOL + 1 ) * ILOG 
    DO WHILE ( JWNO .LT. IWNOH )   ! Include any further mol.data ranges
      IDX = IDX + 1
      ILOG = OLOG10(IDX)
      IMAX = MAX ( IMAX, ILOG )           
      IWNO = JWNO + 1              ! Update lower wno of new data range
      JWNO = JWNO + WAVENO(IDX)    ! & Upper wno
      ISUM = ISUM + ( MIN(JWNO,IWNOH) - IWNO + 1 ) * ILOG  ! & Summation
    END DO
    OPTMOL(IMOL) = FLOAT(IMAX) + FACTOR * FLOAT(ISUM) ! Combine Max & Avg opt.d
  END DO
!
  NREQ = COUNT ( OPTMOL .GT. FLOAT(IMIN) )   ! .gt. excludes OPTMOL=0 if IMIN=0
  ALLOCATE ( OPTLST(NREQ) ) 
  ALLOCATE ( MOLLST(NREQ) ) 
  DO IREQ = 1, NREQ 
    IMOL = MAXLOC ( OPTMOL, 1 ) 
    OPTLST(IREQ) = OPTMOL(IMOL)
    MOLLST(IREQ) = IMOL
    OPTMOL(IMOL) = -1.0                      ! Flag so excluded hereafter
  END DO
!
END SUBROUTINE LSTABS
END MODULE LSTABS_SUB
