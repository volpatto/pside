c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        Pollution problem
c        ODE of dimension 20
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/pollu.f
c
c     This is revision
c     $Id: pollu.F,v 1.2 2006/10/02 10:29:14 testset Exp $
c
c-----------------------------------------------------------------------
      integer function pidate()
      pidate = 20060828
      return
      end
c-----------------------------------------------------------------------
      subroutine prob(fullnm,problm,type,
     +                neqn,ndisc,t,
     +                numjac,mljac,mujac,
     +                nummas,mlmas,mumas,
     +                ind)
      character*(*) fullnm, problm, type
      integer neqn,ndisc,mljac,mujac,mlmas,mumas,ind(*)
      double precision t(0:*)
      logical numjac, nummas

      fullnm = 'Pollution problem'
      problm = 'pollu'
      type   = 'ODE'
      neqn   = 20
      ndisc  = 0
      t(0)   = 0d0
      t(1)   = 60d0
      numjac = .false.
      mljac  = neqn
      mujac  = neqn

      return
      end
c-----------------------------------------------------------------------
      subroutine init(neqn,t,y,yprime,consis)
      integer neqn
      double precision t,y(neqn),yprime(neqn)
      logical consis

      integer i

      do 10 i=1,neqn
         y(i) = 0d0
   10 continue

      y(2)  = 0.2d0
      y(4)  = 0.04d0
      y(7)  = 0.1d0
      y(8)  = 0.3d0
      y(9)  = 0.01d0
      y(17) = 0.007d0

      return
      end
c-----------------------------------------------------------------------
      subroutine settolerances(neqn,rtol,atol,tolvec)
      integer neqn 
      logical tolvec
      double precision rtol(neqn), atol(neqn)
       
      tolvec  = .false.
      

      return
      end
c-----------------------------------------------------------------------
      subroutine setoutput(neqn,solref,printsolout,
     +                    nindsol,indsol)

      logical solref, printsolout
      integer neqn, nindsol
      integer indsol(neqn)

c the reference solution is available
      solref = .true.  

c output file is required
      printsolout = .true.

c default values if printsolout is .true.
      nindsol = neqn
c only nindsol component of indsol are referenced
      do i=1,nindsol
          indsol(i) = i
      end do  

      return
      end

c-----------------------------------------------------------------------
      subroutine feval(neqn,t,y,yprime,f,ierr,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*)

      double precision k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,
     +                 k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,r(25)
      parameter (k1=.35d0,   k2=.266d2,
     +           k3=.123d5,  k4=.86d-3,
     +           k5=.82d-3,  k6=.15d5,
     +           k7=.13d-3,  k8=.24d5,
     +           k9=.165d5,  k10=.9d4,
     +           k11=.22d-1, k12=.12d5,
     +           k13=.188d1, k14=.163d5,
     +           k15=.48d7,  k16=.35d-3,
     +           k17=.175d-1,k18=.1d9,
     +           k19=.444d12,k20=.124d4,
     +           k21=.21d1,  k22=.578d1,
     +           k23=.474d-1,k24=.178d4,
     +           k25=.312d1)

      r( 1) = k1 *y( 1)
      r( 2) = k2 *y( 2)*y(4)
      r( 3) = k3 *y( 5)*y(2)
      r( 4) = k4 *y( 7)
      r( 5) = k5 *y( 7)
      r( 6) = k6 *y( 7)*y(6)
      r( 7) = k7 *y( 9)
      r( 8) = k8 *y( 9)*y(6)
      r( 9) = k9 *y(11)*y(2)
      r(10) = k10*y(11)*y(1)
      r(11) = k11*y(13)
      r(12) = k12*y(10)*y(2)
      r(13) = k13*y(14)
      r(14) = k14*y( 1)*y(6)
      r(15) = k15*y( 3)
      r(16) = k16*y( 4)
      r(17) = k17*y( 4)
      r(18) = k18*y(16)
      r(19) = k19*y(16)
      r(20) = k20*y(17)*y(6)
      r(21) = k21*y(19)
      r(22) = k22*y(19)
      r(23) = k23*y( 1)*y(4)
      r(24) = k24*y(19)*y(1)
      r(25) = k25*y(20)

      f(1)  = -r(1)-r(10)-r(14)-r(23)-r(24)+
     +        r(2)+r(3)+r(9)+r(11)+r(12)+r(22)+r(25)
      f(2)  = -r(2)-r(3)-r(9)-r(12)+r(1)+r(21)
      f(3)  = -r(15)+r(1)+r(17)+r(19)+r(22)
      f(4)  = -r(2)-r(16)-r(17)-r(23)+r(15)
      f(5)  = -r(3)+r(4)+r(4)+r(6)+r(7)+r(13)+r(20)
      f(6)  = -r(6)-r(8)-r(14)-r(20)+r(3)+r(18)+r(18)
      f(7)  = -r(4)-r(5)-r(6)+r(13)
      f(8)  = r(4)+r(5)+r(6)+r(7)
      f(9)  = -r(7)-r(8)
      f(10) = -r(12)+r(7)+r(9)
      f(11) = -r(9)-r(10)+r(8)+r(11)
      f(12) = r(9)
      f(13) = -r(11)+r(10)
      f(14) = -r(13)+r(12)
      f(15) = r(14)
      f(16) = -r(18)-r(19)+r(16)
      f(17) = -r(20)
      f(18) = r(20)
      f(19) = -r(21)-r(22)-r(24)+r(23)+r(25)
      f(20) = -r(25)+r(24)

      return
      end
c-----------------------------------------------------------------------
      subroutine jeval(ldim,neqn,t,y,yprime,dfdy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfdy(ldim,neqn),rpar(*)

      integer i,j
      double precision k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,
     +                 k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25
      parameter (k1=.35d0,   k2=.266d2,
     +           k3=.123d5,  k4=.86d-3,
     +           k5=.82d-3,  k6=.15d5,
     +           k7=.13d-3,  k8=.24d5,
     +           k9=.165d5,  k10=.9d4,
     +           k11=.22d-1, k12=.12d5,
     +           k13=.188d1, k14=.163d5,
     +           k15=.48d7,  k16=.35d-3,
     +           k17=.175d-1,k18=.1d9,
     +           k19=.444d12,k20=.124d4,
     +           k21=.21d1,  k22=.578d1,
     +           k23=.474d-1,k24=.178d4,
     +           k25=.312d1)

      do 20 j=1,neqn
         do 10 i=1,neqn
            dfdy(i,j) = 0d0
   10    continue
   20 continue

      dfdy(1,1)   = -k1-k10*y(11)-k14*y(6)-k23*y(4)-k24*y(19)
      dfdy(1,11)  = -k10*y(1)+k9*y(2)
      dfdy(1,6)   = -k14*y(1)
      dfdy(1,4)   = -k23*y(1)+k2*y(2)
      dfdy(1,19)  = -k24*y(1)+k22
      dfdy(1,2)   = k2*y(4)+k9*y(11)+k3*y(5)+k12*y(10)
      dfdy(1,13)  = k11
      dfdy(1,20)  = k25
      dfdy(1,5)   = k3*y(2)
      dfdy(1,10)  = k12*y(2)
c
      dfdy(2,4)   = -k2*y(2)
      dfdy(2,5)   = -k3*y(2)
      dfdy(2,11)  = -k9*y(2)
      dfdy(2,10)  = -k12*y(2)
      dfdy(2,19)  = k21
      dfdy(2,1)   = k1
      dfdy(2,2)   = -k2*y(4)-k3*y(5)-k9*y(11)-k12*y(10)
c
      dfdy(3,1)   = k1
      dfdy(3,4)   = k17
      dfdy(3,16)  = k19
      dfdy(3,19)  = k22
      dfdy(3,3)   = -k15
c
      dfdy(4,4)   = -k2*y(2)-k16-k17-k23*y(1)
      dfdy(4,2)   = -k2*y(4)
      dfdy(4,1)   = -k23*y(4)
      dfdy(4,3)   = k15
c
      dfdy(5,5)   = -k3*y(2)
      dfdy(5,2)   = -k3*y(5)
      dfdy(5,7)   = 2d0*k4+k6*y(6)
      dfdy(5,6)   = k6*y(7)+k20*y(17)
      dfdy(5,9)   = k7
      dfdy(5,14)  = k13
      dfdy(5,17)  = k20*y(6)
c
      dfdy(6,6)   = -k6*y(7)-k8*y(9)-k14*y(1)-k20*y(17)
      dfdy(6,7)   = -k6*y(6)
      dfdy(6,9)   = -k8*y(6)
      dfdy(6,1)   = -k14*y(6)
      dfdy(6,17)  = -k20*y(6)
      dfdy(6,2)   = k3*y(5)
      dfdy(6,5)   = k3*y(2)
      dfdy(6,16)  = 2d0*k18
c
      dfdy(7,7)   = -k4-k5-k6*y(6)
      dfdy(7,6)   = -k6*y(7)
      dfdy(7,14)  = k13
c
      dfdy(8,7)   = k4+k5+k6*y(6)
      dfdy(8,6)   = k6*y(7)
      dfdy(8,9)   = k7
c
      dfdy(9,9)   = -k7-k8*y(6)
      dfdy(9,6)   = -k8*y(9)
c
      dfdy(10,10) = -k12*y(2)
      dfdy(10,2)  = -k12*y(10)+k9*y(11)
      dfdy(10,9)  = k7
      dfdy(10,11) = k9*y(2)
c
      dfdy(11,11) = -k9*y(2)-k10*y(1)
      dfdy(11,2)  = -k9*y(11)
      dfdy(11,1)  = -k10*y(11)
      dfdy(11,9)  = k8*y(6)
      dfdy(11,6)  = k8*y(9)
      dfdy(11,13) = k11
c
      dfdy(12,11) = k9*y(2)
      dfdy(12,2)  = k9*y(11)
c
      dfdy(13,13) = -k11
      dfdy(13,11) = k10*y(1)
      dfdy(13,1)  = k10*y(11)
c
      dfdy(14,14) = -k13
      dfdy(14,10) = k12*y(2)
      dfdy(14,2)  = k12*y(10)
c
      dfdy(15,1)  = k14*y(6)
      dfdy(15,6)  = k14*y(1)
c
      dfdy(16,16) = -k18-k19
      dfdy(16,4)  = k16
c
      dfdy(17,17) = -k20*y(6)
      dfdy(17,6)  = -k20*y(17)
c
      dfdy(18,17) = k20*y(6)
      dfdy(18,6)  = k20*y(17)
c
      dfdy(19,19) = -k21-k22-k24*y(1)
      dfdy(19,1)  = -k24*y(19)+k23*y(4)
      dfdy(19,4)  = k23*y(1)
      dfdy(19,20) = k25
c
      dfdy(20,20) = -k25
      dfdy(20,1)  = k24*y(19)
      dfdy(20,19) = k24*y(1)

      return
      end
c-----------------------------------------------------------------------
      subroutine meval(ldim,neqn,t,y,yprime,dfddy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfddy(ldim,neqn),rpar(*)
c
c     dummy subroutine
c
      return
      end
c-----------------------------------------------------------------------
      subroutine solut(neqn,t,y)
      integer neqn
      double precision t,y(neqn)
c
c computed using true double precision RADAU5 on Cray C90
c          uround = work(1) = 1.01d-19
c          rtol = atol = h0 = 1.1d-18
c
      y( 1) = 0.5646255480022769d-01
      y( 2) = 0.1342484130422339d+00
      y( 3) = 0.4139734331099427d-08
      y( 4) = 0.5523140207484359d-02
      y( 5) = 0.2018977262302196d-06
      y( 6) = 0.1464541863493966d-06
      y( 7) = 0.7784249118997964d-01
      y( 8) = 0.3245075353396018d+00
      y( 9) = 0.7494013383880406d-02
      y(10) = 0.1622293157301561d-07
      y(11) = 0.1135863833257075d-07
      y(12) = 0.2230505975721359d-02
      y(13) = 0.2087162882798630d-03
      y(14) = 0.1396921016840158d-04
      y(15) = 0.8964884856898295d-02
      y(16) = 0.4352846369330103d-17
      y(17) = 0.6899219696263405d-02
      y(18) = 0.1007803037365946d-03
      y(19) = 0.1772146513969984d-05
      y(20) = 0.5682943292316392d-04

      return
      end
