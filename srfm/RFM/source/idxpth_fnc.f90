MODULE IDXPTH_FNC
CONTAINS
INTEGER(I4) FUNCTION IDXPTH ( ITAN, IATM, IGAS, IDIR )
!
! VERSION
!   06OCT21 AD Store paths to save time. Checked.
!   05MAR19 AD Add TANCOM and use TAN%ITN if no path found for ITAN
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Index in PTHCOM of tan/atm/gas/dir
!   General purpose module
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE PTHCOM_DAT ! Path segment data
    USE TANCOM_DAT ! Tangent path data
!
  IMPLICIT NONE 
  SAVE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: ITAN ! Tangent Ray#
    INTEGER(I4), INTENT(IN) :: IATM ! Atmospheric Layer#
    INTEGER(I4), INTENT(IN) :: IGAS ! Absorber#
    INTEGER(I4), OPTIONAL, &
                 INTENT(IN) :: IDIR ! Direction (ignored unless USEDIR is TRUE)
! LOCAL VARIABLES
    INTEGER(I4) :: IAT        ! Local atmospheric layer index
    INTEGER(I4) :: IATMAX     ! Max.atmos layer# accessed
    INTEGER(I4) :: IATMIN     ! Min.atmos layer# accessed
    INTEGER(I4) :: IGS        ! Local absorber index
    INTEGER(I4) :: IGSMAX     ! Max. absorber# accessed
    INTEGER(I4) :: IPTH       ! Counter
    INTEGER(I4) :: ITN        ! Local tangent path index
    INTEGER(I4) :: ITNMAX = 0 ! Max Tan ray# accessed
    INTEGER(I4) :: JDIR = 1   ! Direction index (1 if not specified)
    INTEGER(I4) :: JTAN       ! Local version of ITAN
    INTEGER(I4) :: LPTH = 0   ! Last value of NPTH
    INTEGER(I4) :: NDIR = 1   ! No. path directions to use 
    INTEGER(I4), ALLOCATABLE :: IDX(:,:,:,:) ! List of saved IDXPTH values
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Assume that if using directions (ie GRA flag) then IDIR=-1 or +1
! which is converted to JDIR = 1 or 2
  IF ( LPTH .NE. NPTH ) THEN ! re-index
    IATMIN = MINVAL ( PTH(:)%IAT ) 
    IATMAX = MAXVAL ( PTH(:)%IAT )
    ITNMAX = MAXVAL ( PTH(:)%ITN ) 
    IGSMAX = MAXVAL ( PTH(:)%IGS )
    IF ( USEDIR ) NDIR = 2
    IF ( ALLOCATED ( IDX ) ) DEALLOCATE ( IDX ) 
    ALLOCATE ( IDX(ITNMAX,IATMIN:IATMAX,IGSMAX,NDIR) )
    IDX = 0
    DO IPTH = 1, NPTH
      ITN = PTH(IPTH)%ITN
      IAT = PTH(IPTH)%IAT
      IGS = PTH(IPTH)%IGS
      IF ( USEDIR ) JDIR = ( PTH(IPTH)%IDR + 3 ) / 2 
      IDX(ITN,IAT,IGS,JDIR) = IPTH
    END DO
    LPTH = NPTH
  END IF
!    
  IF ( USEDIR ) JDIR = ( IDIR + 3 ) / 2    ! IDIR = 0,1 JDIR = 1,2
!
  IF ( ITAN .LE. ITNMAX ) THEN
    IDXPTH = IDX(ITAN,IATM,IGAS,JDIR)
    IF ( IDXPTH .GT. 0 ) RETURN
  END IF
!
  JTAN = TAN(ITAN)%ITN                  ! try alternative (orig) ITAN value
  IDXPTH = IDX(JTAN,IATM,IGAS,JDIR)
  IF ( IDXPTH .GT. 0 ) RETURN
!
! Should always identify a path
  STOP 'F-IDXPTH: Logical error'
!
END FUNCTION IDXPTH
END MODULE IDXPTH_FNC
