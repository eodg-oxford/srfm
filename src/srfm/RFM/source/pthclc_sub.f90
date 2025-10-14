MODULE PTHCLC_SUB
CONTAINS
SUBROUTINE PTHCLC 
!
! VERSION
!   20SEP21 AD Replace ADDCLC by ICLPTH. Checked.
!   01MAY17 AD F90 conversion of rfmsca.for. Checked.
!
! DESCRIPTION
!   Determine which path absorption calculations can be scaled.
!   Called by RFMPTH in limb-viewing mode unless CLC option selected.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE PTHCOM_DAT ! Path segment data
    USE TANCOM_DAT ! Tangent path data
    USE ATMCOM_DAT, ONLY: NATM, NVMR ! No. atm.profile levels, species
    USE CLCCOM_DAT, ONLY: NCLC ! No. calculated path segments
!
! SUBROUTINES
    USE ICLPTH_FNC ! Index of corresponding (new) calculated path
    USE IDXPTH_FNC ! Index in PTHCOM of tan/atm/gas/dir
    USE WRTLOG_SUB ! Write text message to log file
!
  IMPLICIT NONE
!
! LOCAL CONSTANTS
    INTEGER(I4), PARAMETER :: NEXTRA = 1 ! No.layers above tangent layer for 
!                                          explicit line-by-line calculation
! LOCAL VARIABLES
    LOGICAL       :: ADDNEW  ! T=force new calc path to be set 
    INTEGER(I4)   :: IATM    ! Atmospheric profile layer
    INTEGER(I4)   :: IDIR    ! Direction pointer
    INTEGER(I4)   :: IPTH    ! Path counter
    INTEGER(I4)   :: ITAN    ! Tangent path counter
    INTEGER(I4)   :: IVMR    ! Gas counter
    INTEGER(I4)   :: NDIR    ! -1=direction irrelevant, 1=both up and down reqd.
    CHARACTER(80) :: MESSGE  ! Info message sent to Log file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  IF ( USEDIR ) THEN 
    NDIR = 1
  ELSE
    NDIR = -1
  END IF
!
! Set all paths in tangent layer + NEXTRA layers for line-by-line calculation
  DO ITAN = 1, MTAN
    IF ( .NOT. TAN(ITAN)%CLC ) CYCLE
    DO IATM = TAN(ITAN)%IAT, NATM-1
      DO IVMR = 1, NVMR
        DO IDIR = -1, NDIR, 2
          IPTH = IDXPTH ( ITAN, IATM, IVMR, IDIR )
          ADDNEW = IATM .LE. TAN(ITAN)%IAT + NEXTRA
          PTH(IPTH)%ICL = ICLPTH ( PTH(IPTH), ADDNEW ) 
        END DO
      END DO
    END DO
  END DO
! 
  IF ( NPTH .GT. 0 ) THEN
    WRITE ( MESSGE, '(A,I8,A,I8,A,I3,A)' ) &
      'I-PTHCLC: Explicit Absorption calc. for', NCLC, &
      ' paths out of', NPTH, '(=', NINT(100.0*NCLC/FLOAT(NPTH)),'%)'
    CALL WRTLOG ( MESSGE ) 
  END IF
!
END SUBROUTINE PTHCLC
END MODULE PTHCLC_SUB

