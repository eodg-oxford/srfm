! Copyright (C) 1998-2017 University of Oxford
!
! This source code is licensed under the GNU General Public License (GPL),
! Version 3.  See the file COPYING for more details.


!     General purpose Mie scattering routine for single particles reworked into f90
!
!     Author: R Grainger 1990, 2024


      Subroutine mie_ewp(Dx, SCm, Inp, Dqv, Dqxt, Dqsc, Pf,  Error)

      Implicit None
      INTEGER, PARAMETER :: Imaxx   = 105000  ! say 5 mm radius drop, lambda = .3 um (uv)
      INTEGER, PARAMETER :: Itermax = 106000  ! must be large enough to cope with the
                                              ! largest possible nmx = x * abs(scm) + 15 ((x=10500, scm = [2.5, -10] -> nmx = 105026)
                                              ! or nmx =  Dx + 4.05*Dx**(1./3.) + 2.0  (x=10500 -> nmx = 105193)
!     INPUT
      real(8), intent(in) :: Dx
      complex(8), intent(in) :: SCm
      integer(4), intent(in) :: Inp
      real(8), intent(in) :: Dqv(Inp)

!     OUTPUT
      real(8), intent(out) :: Dqxt, Dqsc
      real(8), intent(out) :: PF(Inp)
      integer(4), intent(out) :: Error

!     LOCAL
      Integer * 4  I
      Integer * 4  NStop
      Integer * 4  NmX
      Integer * 4  N ! N*N > 32767 ie N > 181
  
      Real * 8     Dx2
      Real * 8     Chi,Chi0,Chi1
      Real * 8     APsi,APsi0,APsi1
      Real(8), ALLOCATABLE :: Pi0(:)
      Real(8), ALLOCATABLE :: Pi1(:)
      Real * 8     Taun
      Real * 16    Psi,Psi0,Psi1
      Complex * 8  Ir
      Complex * 16 Cm
      Complex * 16 A,ANM1,APB
      Complex * 16 B,BNM1,AMB
      Complex(8), ALLOCATABLE :: D(:)
      Complex(8), ALLOCATABLE :: Sp(:)
      Complex(8), ALLOCATABLE :: Sm(:)
      complex(8)              :: Xs1, Xs2  
      Complex * 16 Xi,Xi0,Xi1
      Complex * 16 Y

!     ACCELERATOR VARIABLES
      Integer * 4  Tnp1
      Integer * 4  Tnm1
      Real * 16    Dn
      Real * 8     Rnx
      Real * 8     S
      Real * 8     T
      Real * 8     Turbo
      Real * 8     AA
      Real * 8     A2
      Complex * 16 A1

     If ((Dx.Gt.Imaxx)) Then
        Error = 1
        Return
      EndIf
      Cm = SCm
      Ir = 1 / Cm

      Y =  Dx * Cm
      If (Dx.Lt.0.02) Then
         NStop = 2
      Else
         If (Dx.Le.8.0) Then
            NStop = Dx + 4.00*Dx**(1./3.) + 1.0
         Else
            If (Dx.Lt. 4200.0) Then
               NStop = Dx + 4.05*Dx**(1./3.) + 2.0
            Else
               NStop = Dx + 4.00*Dx**(1./3.) + 2.0
            End If
         End If 
      End If
      NmX = Max(Real(NStop),Real(Abs(Y))) + 15.
      If (Nmx .Gt. Itermax) Then
          Error = 2
          Return
      End If
      ALLOCATE(Pi0(InP), Pi1(InP), Sp(InP), Sm(InP), D(Nmx))

      D(NmX) = Dcmplx(0,0)
      Do N = Nmx-1,1,-1
         A1 = (N+1) / Y
         D(N) = A1 - 1/(A1+D(N+1))
      End Do
     
      Sm(:Inp) = dcmplx(0d0,0d0) 
      Sp(:Inp) = dcmplx(0d0,0d0)
 
      Pi0(:Inp) = 0.0_8
      Pi1(:Inp) = 1.0_8

      Psi0 = Cos(Dx)
      Psi1 = Sin(Dx)
      Chi0 =-Sin(Dx)
      Chi1 = Cos(Dx)
      APsi0 = Psi0
      APsi1 = Psi1
      Xi0 = Dcmplx(APsi0,Chi0)
      Xi1 = Dcmplx(APsi1,Chi1)
      Dqsc = 0
      Dqxt = 0
      Tnp1 = 1
      Do N = 1,Nstop
         DN = N
         Tnp1 = Tnp1 + 2
         Tnm1 = Tnp1 - 2
         A2 = Tnp1 / (DN*(DN+1D0))
         Turbo = (DN+1D0) / DN
         Rnx = DN/Dx
         Psi = Dble(Tnm1)*Psi1/Dx - Psi0
         APsi = Psi
         Chi = Tnm1*Chi1/Dx       - Chi0
         Xi = Dcmplx(APsi,Chi)
         A = ((D(N)*Ir+Rnx)*APsi-APsi1) / ((D(N)*Ir+Rnx)*  Xi-  Xi1)
         B = ((D(N)*Cm+Rnx)*APsi-APsi1) / ((D(N)*Cm+Rnx)*  Xi-  Xi1)
         Dqxt   = Tnp1 *      Dble(A + B)          + Dqxt
         Dqsc   = Tnp1 * (A*Conjg(A) + B*Conjg(B)) + Dqsc
         Anm1 = A
         Bnm1 = B
         APB = A2 * (A + B)
         AMB = A2 * (A - B)
         Do concurrent (I = 1:Inp)
            S = Dqv(I) * Pi1(I)
            T = S - Pi0(I)
            Taun = N*T - Pi0(I)
            Sp(I) = APB * (Pi1(I) + Taun) + Sp(I)
            Sm(I) = AMB * (Pi1(I) - Taun) + Sm(I)
            Pi0(I) = Pi1(I)
            Pi1(I) = S + T*Turbo
         End Do
         Psi0 = Psi1
         Psi1 = Psi
         Apsi1 = Psi1
         Chi0 = Chi1
         Chi1 = Chi
         Xi1 = Dcmplx(APsi1,Chi1)
      End Do
      Dx2 = Dx**2
      Dqsc =  2 * Dqsc / Dx2
      Dqxt =  2 * Dqxt / Dx2
      AA = 2 / (Dx2 * Dqsc)
      Do concurrent (I = 1:Inp)
         Xs1 = 0.5 * (Sp(I)+Sm(I)) 
         Xs2 = 0.5 * (Sp(I)-Sm(I)) 
         PF(I) =  AA * Dble( Xs1*Conjg(Xs1) +Xs2*Conjg(Xs2))
      End Do
      DEALLOCATE(Pi0, Pi1, Sp, Sm, D)
      Error = 0
      Return
  end 

