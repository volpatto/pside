c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        Problem HIRES
c        ODE of dimension 8
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/hires.f
c
c     This is revision
c     $Id: hires.F,v 1.2 2006/10/02 10:29:14 testset Exp $
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

      fullnm = 'Problem HIRES'
      problm = 'hires'
      type   = 'ODE'
      neqn   = 8
      ndisc  = 0
      t(0)   = 0d0
      t(1)   = 321.8122d0
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

      y(1) = 1d0
      y(2) = 0d0
      y(3) = 0d0
      y(4) = 0d0
      y(5) = 0d0
      y(6) = 0d0
      y(7) = 0d0
      y(8) = 0.0057d0

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

      f(1) = -1.71d0*y(1)+0.43d0*y(2)+8.32d0*y(3)+0.0007d0
      f(2) = 1.71d0*y(1)-8.75d0*y(2)
      f(3) = -10.03d0*y(3)+0.43d0*y(4)+0.035d0*y(5)
      f(4) = 8.32d0*y(2)+1.71d0*y(3)-1.12d0*y(4)
      f(5) = -1.745d0*y(5)+0.43d0*(y(6)+y(7))
      f(6) = -280d0*y(6)*y(8)+0.69d0*y(4)+1.71d0*y(5)-0.43d0*y(6)+
     +        0.69d0*y(7)
      f(7) = 280d0*y(6)*y(8)-1.81d0*y(7)
      f(8) = -f(7)

      return
      end
c-----------------------------------------------------------------------
      subroutine jeval(ldim,neqn,t,y,yprime,dfdy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfdy(ldim,neqn),rpar(*)

      integer i,j

      do 20 j=1,neqn
         do 10 i=1,neqn
            dfdy(i,j)=0d0
   10    continue
   20 continue

      dfdy(1,1) = -1.71d0
      dfdy(1,2) = 0.43d0
      dfdy(1,3) = 8.32d0
      dfdy(2,1) = 1.71d0
      dfdy(2,2) = -8.75d0
      dfdy(3,3) = -10.03d0
      dfdy(3,4) = 0.43d0
      dfdy(3,5) = 0.035d0
      dfdy(4,2) = 8.32d0
      dfdy(4,3) = 1.71d0
      dfdy(4,4) = -1.12d0
      dfdy(5,5) = -1.745d0
      dfdy(5,6) = 0.43d0
      dfdy(5,7) = 0.43d0
      dfdy(6,4) = 0.69d0
      dfdy(6,5) = 1.71d0
      dfdy(6,6) = -280d0*y(8)-0.43d0
      dfdy(6,7) = 0.69d0
      dfdy(6,8) = -280d0*y(6)
      dfdy(7,6) = 280d0*y(8)
      dfdy(7,7) = -1.81d0
      dfdy(7,8) = 280d0*y(6)
      dfdy(8,6) = -280d0*y(8)
      dfdy(8,7) = 1.81d0
      dfdy(8,8) = -280d0*y(6)

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
      y(1) = 0.7371312573325668d-3
      y(2) = 0.1442485726316185d-3
      y(3) = 0.5888729740967575d-4
      y(4) = 0.1175651343283149d-2
      y(5) = 0.2386356198831331d-2
      y(6) = 0.6238968252742796d-2
      y(7) = 0.2849998395185769d-2
      y(8) = 0.2850001604814231d-2

      return
      end
