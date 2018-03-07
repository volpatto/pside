c----------------------------------------------------------------------- 
c 
c     This file is part of the Test Set for IVP solvers 
c     http://www.dm.uniba.it/~testset/ 
c 
c        generic PSIDE driver 
c 
c     DISCLAIMER: see 
c     http://www.dm.uniba.it/~testset/disclaimer.php 
c 
c     The most recent version of this source file can be found at 
c     http://www.dm.uniba.it/~testset/src/drivers/psided.f 
c 
c     This is revision 
c     $Id: psided.f,v 1.6 2006/10/02 10:19:09 testset Exp $ 
c 
c----------------------------------------------------------------------- 
 
      program psided 
 
      integer md 
      parameter (md=400) 
      integer lrwork, liwork 
      parameter (lrwork = 20 + 
     +                    6*(md*4) + 
     +                    (3*md+1)*md*4 + 
     +                    2*(2*md+1)*md + 
     +                    3*md, 
     +           liwork = 20 + 
     +                    md*4) 
 
      integer neqn,ndisc,nlj,nuj,nlm,num,ind(md), 
     +        iwork(20+4*md),ierr,ipar(2),idid 
      double precision y(md),dy(md),t(0:100), 
     +                 h0,rtol(md),atol(md), 
     +                 rwork(lrwork),rpar(1) 
      logical numjac,nummas,consis,tolvec 
 
      integer ldmas,nlmas,numas 
      double precision yref(md) 
      character fullnm*40, problm*8, type*3 
      character driver*8, solver*8 
      parameter (driver = 'psided', solver='PSIDE') 
      external ideg,idej,idem 
      external odeg,odej,odem 
      external daeg,daej,daem 
 
      double precision solmax 
      real gettim, timtim, cputim 
      external gettim 
      double precision scd,mescd 
 
      double precision mas(md*md) 
      common /mascom/ mas,ldmas,nlmas,numas 
 
      character  fileout*140,  filepath*100 
      character formatout*30,  namefile*100 
      logical printsolout, solref, printout 
      integer nindsol, indsol(md) 
 
      integer i 
      do 10 i=1,20 
         iwork(i) = 0 
         rwork(i) = 0d0 
   10 continue 
 
c----------------------------------------------------------------------- 
c     check the problem definition interface date 
c----------------------------------------------------------------------- 
 
      call chkdat(driver,solver,20060828) 
 
c----------------------------------------------------------------------- 
c     get the problem dependent parameters 
c----------------------------------------------------------------------- 
 
      call prob(fullnm,problm,type, 
     +          neqn,ndisc,t, 
     +          numjac,nlj,nuj, 
     +          nummas,nlm,num, 
     +          ind) 
      if (type.eq.'IDE') then 
         do 20 i=1,neqn 
            if (ind(i).eq.0) ind(i)=1 
            if (ind(i).gt.1) iwork(2)=1 
   20    continue 
      elseif (type.eq.'ODE') then 
         nummas = .false. 
         nlm    = 0 
         num    = 0 
      elseif (type.eq.'DAE') then 
         nummas = .false. 
         do 30 i=1,neqn 
            if (ind(i).eq.0) ind(i)=1 
            if (ind(i).gt.1) iwork(2)=1 
   30    continue 
         if (nlm.lt.neqn) then 
            ldmas = nlm + num + 1 
         else 
            ldmas = neqn 
         endif 
         nlmas = nlm 
         numas = num 
         call meval(ldmas,neqn,t(0),y,dy,mas,ierr,rpar,ipar) 
      else 
         print *, 'PSIDED: ERROR: ', 
     +            'unknown Test Set problem type', type 
         stop 
      endif 
 
c----------------------------------------------------------------------- 
c     get the initial values 
c----------------------------------------------------------------------- 
 
      call init(neqn,t(0),y,dy,consis) 
      if (type.eq.'ODE') then 
         ierr =  0 
         call feval(neqn,t(0),y,dy,dy,ierr,rpar,ipar) 
         if (ierr.ne.0) then 
            print *, 'PSIDED: ERROR: ', 
     +               'PSIDED could not evaluate f(t(0))' 
            stop 
         endif 
         consis = .true. 
      endif 
 
      if (.not.consis) then 
         print *, 'PSIDED: ERROR: ', 
     +            'PSIDE can not handle inconsistent initial values' 
         stop 
      endif 
 
c----------------------------------------------------------------------- 
c     read the tolerances 
c----------------------------------------------------------------------- 
 
      call getinp(driver,problm,solver,fullnm, 
     +            tolvec,rtol,atol,h0,solmax) 
      
      call settolerances(neqn,rtol,atol,tolvec) 
       
      if (tolvec) iwork(1)=1 
 
      call  setoutput(neqn,solref,printsolout, 
     +                nindsol,indsol) 
 
 
      if (printsolout) then 
          print *, 'PSIDED: WARNING: ', 
     +            'this option is not yet implemented;', 
     +            'we set printsolout  = .false. ' 
          printsolout = .false. 
      end if 
 
c----------------------------------------------------------------------- 
c     call of the subroutine PSIDE 
c----------------------------------------------------------------------- 
 
      timtim = gettim() 
      timtim = gettim() - timtim 
 
      cputim = gettim() 
 
      do 40 i=0,ndisc 
         rwork(1)=0d0 
         if (type.eq.'IDE') then 
            call pside(neqn,y,dy,ideg, 
     +                 numjac,nlj,nuj,idej, 
     +                 nummas,nlm,num,idem, 
     +                 t(i),t(i+1),rtol,atol,ind, 
     +                 lrwork,rwork,liwork,iwork,rpar,ipar, 
     +                 idid) 
         elseif (type.eq.'ODE') then 
            call pside(neqn,y,dy,odeg, 
     +                 numjac,nlj,nuj,odej, 
     +                 nummas,nlm,num,odem, 
     +                 t(i),t(i+1),rtol,atol,ind, 
     +                 lrwork,rwork,liwork,iwork,rpar,ipar, 
     +                 idid) 
         elseif (type.eq.'DAE') then 
            call pside(neqn,y,dy,daeg, 
     +                 numjac,nlj,nuj,daej, 
     +                 nummas,nlm,num,daem, 
     +                 t(i),t(i+1),rtol,atol,ind, 
     +                 lrwork,rwork,liwork,iwork,rpar,ipar, 
     +                 idid) 
         endif 
         if (idid.eq.-1) then 
            print *, 'PSIDED: ERROR: ', 
     +               'PSIDE has stopped because stepsize has become ', 
     +               'too small' 
            stop 
         elseif (idid.eq.-2) then 
            print *, 'PSIDED: ERROR ' 
            stop 
         elseif (idid.ne.1) then 
            print *, 'PSIDED: ERROR: ', 
     +               'unknown PSIDE return code: IDID = ', idid 
            stop 
         endif 
   40 continue 
 
      cputim = gettim() - cputim - timtim 
 
c----------------------------------------------------------------------- 
c     print numerical solution in endpoint and integration statistics 
c----------------------------------------------------------------------- 
      printout = .true. 
      if (solref) then  
         call solut(neqn,t(ndisc+1),yref) 
         call getscd(mescd,scd,neqn,yref,y,problm,tolvec,atol,rtol, 
     +            printout) 
      else 
        call printsol(neqn,y,problm) 
      end if 
          
      call report( 
     +   driver,problm,solver, 
     +   rtol,atol,h0,solmax, 
     +   iwork,cputim,scd,mescd 
     +) 
      end 
 
C======================================================================= 
C     `Test Set for IVP Solvers' IDE wrappers for PSIDE 
C======================================================================= 
c 
c     since in PSIDE the format of the subroutines for the 
c     function G and its derivatives differ from the format 
c     in the testset, we transform them 
c 
c          G    = f        -> ideg 
c 
c        dG/dy  = df/dy    -> idej 
c        dG/dy' = df/dy'   -> idem 
c 
c----------------------------------------------------------------------- 
      subroutine ideg(neqn,t,y,dy,g,ierr,rpar,ipar) 
      integer neqn,ierr,ipar(*) 
      double precision t,y(neqn),dy(neqn),g(neqn),rpar(*) 
      call feval(neqn,t,y,dy,g,ierr,rpar,ipar) 
      return 
      end 
 
      subroutine idej(ldj,neqn,nlj,nuj,t,y,dy,dgdy, 
     +                rpar,ipar) 
      integer ldj,neqn,nlj,nuj,ierr,ipar(*) 
      double precision t,y(neqn),dy(neqn),dgdy(ldj,neqn),rpar(*) 
      ierr = 0 
      call jeval(ldj,neqn,t,y,dy,dgdy,ierr,rpar,ipar) 
      if (ierr.ne.0) then 
         print *, 'PSIDED: ERROR: ', 
     +            'PSIDE can not handle JEVAL IERR' 
         stop 
      endif 
      return 
      end 
 
      subroutine idem(ldm,neqn,nlm,num,t,y,dy,dgddy, 
     +                rpar,ipar) 
      integer ldm,neqn,nlm,num,ierr,ipar(*) 
      double precision t,y(neqn),dy(neqn),dgddy(ldm,neqn),rpar(*) 
      ierr = 0 
      call meval(ldm,neqn,t,y,dy,dgddy,ierr,rpar,ipar) 
      if (ierr.ne.0) then 
         print *, 'PSIDED: ERROR: ', 
     +            'PSIDE can not handle JEVAL IERR' 
         stop 
      endif 
      return 
      end 
 
C======================================================================= 
C     `Test Set for IVP Solvers' ODE wrappers for PSIDE 
C======================================================================= 
c 
c     since in PSIDE the format of the subroutines for the 
c     function G and its derivatives differ from the format 
c     in the testset, we transform them 
c 
c          G    = f - y'   -> odeg 
c 
c        dG/dy  = df/dy    -> odej 
c        dG/dy' = -I       -> odem 
c 
c----------------------------------------------------------------------- 
      subroutine odeg(neqn,t,y,dy,g,ierr,rpar,ipar) 
      integer neqn,ierr,ipar(*) 
      double precision t,y(neqn),dy(neqn),g(neqn),rpar(*) 
      integer i 
c compute G = f - y' 
      call feval(neqn,t,y,dy,g,ierr,rpar,ipar) 
      do 10 i=1,neqn 
         g(i) = g(i) - dy(i) 
   10 continue 
      return 
      end 
 
      subroutine odej(ldj,neqn,nlj,nuj,t,y,dy,dgdy, 
     +                rpar,ipar) 
      integer ldj,neqn,nlj,nuj,ierr,ipar(*) 
      double precision t,y(neqn),dy(neqn),dgdy(ldj,neqn),rpar(*) 
c compute dG/dy = df/dy 
      ierr = 0 
      call jeval(ldj,neqn,t,y,dy,dgdy,ierr,rpar,ipar) 
      if (ierr.ne.0) then 
         print *, 'PSIDED: ERROR: ', 
     +            'PSIDE can not handle JEVAL IERR' 
         stop 
      endif 
      return 
      end 
 
      subroutine odem(ldm,neqn,nlm,num,t,y,dy,dgddy, 
     +                rpar,ipar) 
      integer ldm,neqn,nlm,num,ipar(*) 
      double precision t,y(neqn),dy(neqn),dgddy(ldm,neqn),rpar(*) 
      integer j 
c initialize dG/dy' = -I 
      do 10 j=1,neqn 
         dgddy(1,j) = -1d0 
   10 continue 
      return 
      end 
C======================================================================= 
C     `Test Set for IVP Solvers' DAE wrappers for PSIDE 
C======================================================================= 
c 
c     since in PSIDE the format of the subroutines for the 
c     function G and its derivatives differ from the format 
c     in the testset, we transform them 
c 
c          G    = f - My'   -> daeg 
c 
c        dG/dy  = df/dy     -> daej 
c        dG/dy' = -M        -> daem 
c 
c----------------------------------------------------------------------- 
      subroutine daeg(neqn,t,y,dy,g,ierr,rpar,ipar) 
      integer neqn,ierr,ipar(*) 
      double precision t,y(neqn),dy(neqn),g(neqn),rpar(*) 
      integer i,j 
      integer md,ldmas,nlmas,numas 
      parameter (md=400) 
      double precision mas(md*md) 
      common /mascom/ mas,ldmas,nlmas,numas 
c compute f in g 
      call feval(neqn,t,y,dy,g,ierr,rpar,ipar) 
c compute g := - M*dy + f 
      if (nlmas.ne.neqn) then 
         do 20 j=1,neqn 
            do 10 i=max(1,j-numas ),min(neqn,j+nlmas) 
               g(i) = g(i) - mas((j-1)*ldmas+numas+1-j+i)*dy(j) 
   10       continue 
   20    continue 
      else 
         do 40 j=1,neqn 
            do 30 i=1,neqn 
               g(i) = g(i) - mas((j-1)*ldmas+i)*dy(j) 
   30       continue 
   40    continue 
      endif 
      return 
      end 
 
      subroutine daej(ldj,neqn,nlj,nuj,t,y,dy,dgdy, 
     +                rpar,ipar) 
      integer ldj,neqn,nlj,nuj,ierr,ipar(*) 
      double precision t,y(neqn),dy(neqn),dgdy(ldj,neqn),rpar(*) 
c compute dG/dy = df/dy 
      ierr = 0 
      call jeval(ldj,neqn,t,y,dy,dgdy,ierr,rpar,ipar) 
      if (ierr.ne.0) then 
         print *, 'PSIDED: ERROR: ', 
     +            'PSIDE can not handle JEVAL IERR' 
         stop 
      endif 
      return 
      end 
 
      subroutine daem(ldm,neqn,nlm,num,t,y,dy,dgddy, 
     +                rpar,ipar) 
      integer ldm,neqn,nlm,num,ipar(*) 
      double precision t,y(neqn),dy(neqn),dgddy(ldm,neqn),rpar(*) 
      integer i,j 
      integer md,ldmas,nlmas,numas 
      parameter (md=400) 
      double precision mas(md*md) 
      common /mascom/ mas,ldmas,nlmas,numas 
c compute dG/dy' = -M 
      if (nlm.ne.neqn) then 
         do 20 j=1,neqn 
            do 10 i=max(1,         num+1 +    1-j), 
     +              min(nlm+num+1, num+1 + neqn-j) 
               dgddy(i,j) = -mas((j-1)*ldmas+i) 
   10       continue 
   20    continue 
      else 
         do 40 j=1,neqn 
            do 30 i=1,neqn 
               dgddy(i,j) = -mas((j-1)*ldmas+i) 
   30       continue 
   40    continue 
      endif 
      return 
      end 
