MODULE CTMCO2_SUB
CONTAINS
SUBROUTINE CTMCO2 ( ILBL ) 
!
! VERSION
!   01SEP24 AD Checked. 
!   16NOV17 AD F90 conversion. Checked.
! 
! DESCRIPTION    
!   CO2 continuum
!   Called by SPCCTM for each path containing CO2.  
!   Continuum is calculated on the widemesh grid.
!     
! REFERENCES
!   CLOUGH S.A.,KNEIZYS F.X.,ROTHMAN L.S.,GALLERY W.O. (1981), 
!     Atmospheric spectral transmittance and radiance: FASCOD1B
!     SPIE 277 Atmospheric Transmission  152-166
!
!   CLOUGH S.A.,KNEIZYS F.X.,ROTHMAN L.S.,ANDERSON G.P.,SHETTLE E.P. (1987)
!     Current issues in infrared atmospheric transparency
!     International meeting on Atmospheric Transparency for Satellite 
!     Applications, 15-19 Sept. 1986 Capri, Italy. Ed. G.V. Silvestrini. CUEN.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE CLCCOM_DAT ! Calculated path segments
    USE CTMCO2_DAT ! CO2 continuum data
    USE WIDCOM_DAT ! Widemesh data
!
! SUBROUTINES
    USE LKPIDX_FNC ! Interpolate by real array index
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: ILBL ! Index of LBL calc path
!
! LOCAL CONSTANTS
    REAL(R4), PARAMETER :: PREREF = 1.0 ! Ref. pressure [atm]
    REAL(R4), PARAMETER :: WGT230 = 1.0/(296.0-230.0)/(260.0-230.0) !  5.05E-4
    REAL(R4), PARAMETER :: WGT260 = 1.0/(296.0-260.0)/(230.0-260.0) ! -9.25E-4
    REAL(R4), PARAMETER :: WGT296 = 1.0/(260.0-296.0)/(230.0-296.0) !  4.21E-4
!
! LOCAL VARIABLES      
    INTEGER(I4) :: ICLC   ! Index of calc path segment
    INTEGER(I4) :: IQAD   ! Counter for quadratic interp.pts (1:3)
    INTEGER(I4) :: IWID   ! Counter for widemesh intervals (1:NWID)
    INTEGER(I4) :: IWD2   ! Counter for half-wide-mesh grid (0:NWID*2)
    REAL(R4)    :: CTW    ! Continuum absorption coeff interp to Temp, Wno.
    REAL(R4)    :: CW230,CW260,CW296 ! Temp. Coeffs interp. to Half-WM points
    REAL(R4)    :: DT230,DT260,DT296 ! Diff.of path temp from each tab. temp
    REAL(R4)    :: XW     ! Position of half-WM pt on CO2 continuum tab. axis
    REAL(R4)    :: CTMPTH(0:NWD2) ! Continuum absorption for current path
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  ICLC = IDXLBL(ILBL)
!
  DT230 = CLC(ICLC)%TEM - 230.0
  DT260 = CLC(ICLC)%TEM - 260.0
  DT296 = CLC(ICLC)%TEM - 296.0
!
  DO IWD2 = 0, NWD2
    XW = 1.0 + ( SNGL ( WNOWD2(IWD2) ) - WNOLOW ) / DELWNO 
!
! Interpolate in wavenumber for each tabulated temperature
    CW296 = LKPIDX ( XW, CO2296 )  ! LKPIDX returns 0 if XW out of range 1:N
    CW260 = LKPIDX ( XW, CO2260 ) 
    CW230 = LKPIDX ( XW, CO2230 )
!
! Lagrangian interpolation in temperature
    CTW   =  WGT230 * DT260 * DT296 * CW230 &
           + WGT260 * DT230 * DT296 * CW260 &
           + WGT296 * DT230 * DT260 * CW296 
    CTMPTH(IWD2) = CTW
  END DO
  CTMPTH = CTMPTH * CLC(ICLC)%AMT * ( CLC(ICLC)%PRE / PREREF )
!
  DO IWID = 1, NWID
    IWD2 = 2 * IWID - 3       ! -1, 1, 3, ...
    DO IQAD = 1, 3
      IWD2 = IWD2 + 1         ! 0,1,2,  2,3,4,  4,5,6
      ABSWID(IQAD,IWID,ILBL) = ABSWID(IQAD,IWID,ILBL) + CTMPTH(IWD2)
      IF ( USECNT ) &
        CNTWID(IQAD,IWID,ILBL) = CNTWID(IQAD,IWID,ILBL) + CTMPTH(IWD2) 
    END DO
  END DO
!
END SUBROUTINE CTMCO2
END MODULE CTMCO2_SUB

