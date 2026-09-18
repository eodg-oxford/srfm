MODULE FOVCOM_DAT
!
! VERSION
!   22JAN26 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Field of View data
!   Loaded by FOVFIL.
!   This can either be as relative tangent altitudes or elevation angles.
!   Default NFOV=1 for nadir viewing where no FOV function specified
!   although NFOV>1 is also used to model cloudy fractions of FOV.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
  SAVE
!
  TYPE :: FOVTYP
    REAL(R4) :: ALT ! Rel.alts [km] of FOV fn.
    REAL(R8) :: FNC ! Normalised FOV function
  END TYPE FOVTYP
!
! GLOBAL VARIABLES
    TYPE(FOVTYP), ALLOCATABLE :: FOV(:)
!
    LOGICAL     :: ELEFOV = .FALSE. ! T=stored elevation angles [deg]
    INTEGER(I4) :: NFOV   = 1       ! No. FOV points stored
    INTEGER(I4), ALLOCATABLE :: ITNFOV(:,:) ! [NTAN,NFOV] Tan.ray, FOV
!
END MODULE FOVCOM_DAT
