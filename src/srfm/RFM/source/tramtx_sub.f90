MODULE TRAMTX_SUB
CONTAINS
SUBROUTINE TRAMTX 
!
! VERSION
!   28JUL26 AD Checked.
!   01AUG25 AD Remove LEVCOM_DAT. Calculate KTAN,LTAN locally
!   05MAY24 AD Checked.
!   05MAR19 AD Remove FLXEFN. Checked.
!   01MAY17 AD F90 conversion of part of rfmflx.for. Checked.
!
! DESCRIPTION
!   Calculate transmittance matrix
!   Called by SPCFLX if MTX and (TRA or ABS) flags enabled.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE FULCOM_DAT ! Full grid data
    USE QADCOM_DAT ! Gaussian quadrature data
    USE FINCOM_DAT, ONLY: NFIN ! No. of fine mesh grid points
    USE PHYCON_DAT, ONLY: PI
    USE TANCOM_DAT, ONLY: NTAN ! Tangent path data
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4) :: IQAD   ! Counter for quadrature points
    INTEGER(I4) :: ITAN   ! Counter for output levels
    INTEGER(I4) :: JTAN   ! Secondary counter for output levels
    INTEGER(I4) :: KTAN   ! Index of ray path for matrix element
    INTEGER(I4) :: LTAN   ! Index of ray path for transpose matrix element
    REAL(R8)    :: OPT(NFIN)        ! Cumulative optical path
    REAL(R8)    :: TQAD(NFIN,NQAD)  ! Transmittances for quadrature paths
!       
! EXECUTABLE CODE -------------------------------------------------------------
!
  DO ITAN = 1, NTAN
    DO JTAN = ITAN, NTAN
      KTAN = NTAN + (ITAN-1)*NTAN + JTAN
      OPT = OPTFUL(IFUL1:IFUL2,JTAN) - OPTFUL(IFUL1:IFUL2,ITAN) 
      DO IQAD = 1, NQAD
        TQAD(:,IQAD) = EXP ( - ABS ( OPT ) / XQAD(IQAD) ) / PI
      END DO
      TRAFUL(IFUL1:IFUL2,KTAN) = MATMUL ( TQAD, WQAD ) 
      LTAN = NTAN + (JTAN-1)*NTAN + ITAN
      TRAFUL(IFUL1:IFUL2,LTAN) = TRAFUL(IFUL1:IFUL2,KTAN) 
    END DO
  END DO
!
END SUBROUTINE TRAMTX
END MODULE TRAMTX_SUB
