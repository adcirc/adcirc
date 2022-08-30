C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 41.09                                    *
C                                                                             *
C  mod history                                                                *
C  v43.03     - 05/20/03 - rl - from 43.02 - parallel wind stuff (m.brown)    *
C                                          output buffer flush (m.cobb)       *
C                                          3D fixes (k.dresback)              *
C                                          drop MNPROC in fort.15 (t.campbell)*
C                                          various bug fixes in RBCs          *
c                                          ZSURFBUOY/BCPG calc                *
C  v40.02m002 - 12/22 - jjw/vjp - Vic suggested to avoid compiler conflicts   *
C****************************************************************************** 
C 
      MODULE ITPACKV
      use AMPI_LUN_VIRTUALIZED, only : AMPI_LUN
      use AMPI_LUN_migratable
C
C  vjp  9/19/99
c
c-----------------------------------------------------------------------------
c     version -  itpackv 2d (january 1990)
c
c     code written by - david kincaid, roger grimes, john respess
c                       center for numerical analysis
c                       university of texas
c                       austin, tx  78712
c                       (512) 471-1242
c-----------------------------------------------------------------------------
c

#ifdef CMPI
      USE MESSENGER                !   MPI interface for padcirc
#else
      USE SIZES
#endif


c--------------------data declarations end here-------------------------------c


      CONTAINS

      subroutine jcg (n,ndim,maxnz,jcoef,coef,rhs,u,iwksp,nw,wksp,
     *                iparm,rparm,ier)
      implicit none
      integer n,ndim,maxnz,nw,ier,ib1,ib2,ib3,ib4,ib5,nb,n3,i,nbo,
     *  itmax1,loop,idgts
      integer j, nandex, lnans, gnans, iimax !jgf46.00 add
      integer jcoef(ndim,maxnz),iwksp(3*n),iparm(12)
      real(sz) timi1,timj1,temp,time1,time2,timi2,timj2,digit1,
     *   digit2,tol
      real(sz) coef(ndim,maxnz),rhs(n),u(n),wksp(nw),rparm(12)
c
c *** begin -- itpackv common
c
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
c
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
c
      real(sz) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
c *** end   -- itpackv common
c
c     itpackv 2d main routine  jcg  (jacobi conjugate gradient)
c     each of the main routines --
c           jcg, jsi, sor, ssorcg, ssorsi, rscg, rssi
c     can be used independently of the others
c
c ... function --
c
c          jcg drives the jacobi conjugate gradient algorithm.
c
c ... parameter list --
c
c          n      input integer.  dimension of the matrix.
c          ndim   row dimension of jcoef and coef arrays in calling
c                   routine
c          maxnz  maximum number of nonzeros per row
c          jcoef  integer array for sparse matrix representation.
c          coef   array for sparse matrix representation.
c                 jcoef and coef use the ellpack data structure.
c          rhs    input vector.  contains the right hand side
c                 of the matrix problem.
c          u      input/output vector.  on input, u contains the
c                 initial guess to the solution. on output, it contains
c                 the latest estimate to the solution.
c          iwksp  integer vector workspace of length 3*n
c          nw     input integer.  length of available wksp.  on output,
c                 iparm(8) is amount used.
c          wksp   vector used for working space.  jacobi conjugate
c                 gradient needs this to be in length at least
c                 4*n + 4*itmax.  here itmax = iparm(1) is the
c                 maximum allowable number of iterations.
c          iparm  integer vector of length 12.  allows user to specify
c                 some integer parameters which affect the method.
c          rparm  vector of length 12. allows user to specify some
c                 parameters which affect the method.
c          ier    output integer.  error flag.
c
c ... jcg module references --
c
c         from itpackv    chgcon, determ, dfault, echall,
c                         eigvns, eigvss, eqrt1s, iterm ,
c                         itjcg , parcon, permat, peror,
c                         pervec, pjac  , pmult , prbndx, pstop ,
c                         sbelm , scal  , unscal,
c                         vout  , zbrent
c          system         abs, log10, amax0, amax1, mod, sqrt
c
c ... local itpackv references --
c
c          echall, itjcg , permat,
c          peror, pervec, pjac  , prbndx, sbelm , scal  , unscal
c
c     version -  itpackv 2d (january 1990)
c
c     code written by - david kincaid, roger grimes, john respess
c                       center for numerical analysis
c                       university of texas
c                       austin, tx  78712
c                       (512) 471-1242
c
c     for additional details on the
c          (a) routine    see toms article 1982
c          (b) algorithm  see cna report 150
c
c     based on theory by - david young, david kincaid, lou hageman
c
c     reference the book - applied iterative methods
c                          l. hageman, d. young
c                          academic press, 1981
c
c     **************************************************
c     *               important note                   *
c     *                                                *
c     *      when installing itpackv routines on a     *
c     *  different computer, reset some of the values  *
c     *  in  subroutne dfault.   most important are    *
c     *                                                *
c     *   srelpr      machine relative precision       *
c     *   rparm(1)    stopping criterion               *
c     *                                                *
c     *   also change system-dependent routine         *
c     *                                                *
c     **************************************************
c
c ... variables in common block - itcom1
c
c     in     - iteration number
c     is     - iteration number when parameters last changed
c     isym   - symmetric/nonsymmetric case switch
c     itmax  - maximum number of iterations allowed
c     level  - level of output control switch
c     nout   - output unit number
c
c ... variables in common block - itcom2
c
c     adapt  - fully adaptive procedure switch
c     betadt - switch for adaptive determination of beta
c     caseii - adaptive procedure case switch
c     halt   - stopping test switch
c     partad - partially adaptive procedure switch
c
c ... variables in common block - itcom3
c
c     bdelnm - two norm of b times delta-super-n
c     betab  - estimate for the spectral radius of lu matrix
c     cme    - estimate of largest eigenvalue
c     delnnm - inner product of pseudo-residual at iteration n
c     delsnm - inner product of pseudo-residual at iteration s
c     ff     - adaptive procedure damping factor
c     gamma  - acceleration parameter
c     omega  - overrelaxation parameter for sor and ssor
c     qa     - pseudo-residual ratio
c     qt     - virtual spectral radius
c     rho    - acceleration parameter
c     rrr    - adaptive parameter
c     sige   - parameter sigma-sub-e
c     sme    - estimate of smallest eigenvalue
c     specr  - spectral radius estimate for ssor
c     srelpr - machine relative precision
c     stptst - stopping parameter
c     udnm   - two norm of u
c     zeta   - stopping criterion
c
c ... initialize common blocks
c
      ier = 0
      level = iparm(2)
      nout = iparm(4)
      if (level.ge.1) write (AMPI_LUN(nout),10)
   10 format (///1x,'i t p a c k      j c g      ')
      if (iparm(1).le.0) go to 370
c     if (iparm(11).eq.0) timj1 = timer(0.0)
      call echall (n,ndim,maxnz,jcoef,coef,rhs,iparm,rparm,1)
C      temp = 500.0*srelpr  !jgf46.00 commented out
      temp = 512.0D0*srelpr !jgf46.00 added
      if (zeta.ge.temp) go to 30
      if (level.ge.1) write (AMPI_LUN(nout),20) zeta,srelpr,temp
   20 format (/1x,'*** w a r n i n g ************'//1x,
     *      '    in itpackv routine jcg'/1x,
     *      '    rparm(1) =',e10.3,' (zeta)'/1x,
     *      '    a value this small may hinder convergence '/1x,
     *      '    since machine precision srelpr =',e10.3/1x,
     *      '    zeta reset to ',e10.3)
      zeta = temp
   30 continue
      time1 = rparm(9)
      time2 = rparm(10)
      digit1 = rparm(11)
      digit2 = rparm(12)
c
c ... verify n
c
      if (n.gt.0) go to 50
      ier = 11
      if (level.ge.0) write (AMPI_LUN(nout),40) n
   40 format (/1x,'*** f a t a l     e r r o r ************'//1x,
     *        '    called from itpackv routine jcg '/1x,
     *        '    invalid matrix dimension, n =',i8)
      go to 370
c
c ... scale linear system, u, and rhs by the square root of the
c ... diagonal elements.
c
   50 continue
      call scal (n,ndim,maxnz,jcoef,coef,rhs,u,wksp,ier)
      if (ier.eq.0) go to 70
      if (level.ge.0) write (AMPI_LUN(nout),60) ier
   60 format (/1x,'*** f a t a l     e r r o r ************'//1x,
     *   '    called from itpackv routine jcg '/1x,
     *   '    error detected in routine  scal  '/1x,
     *   '    which scales the system   '/1x,
     *   '    ier = ',i5)
      go to 370
c
c ... remove rows and columns if requested
c
   70 continue
      if (iparm(10).eq.0) go to 80
      tol = rparm(8)
      call sbelm (n,ndim,maxnz,jcoef,coef,rhs,wksp,tol)
c
c ... initialize wksp base addresses.
c
   80 ib1 = 1
      ib2 = ib1+n
      ib3 = ib2+n
      ib4 = ib3+n
      ib5 = ib4+n
c
c ... permute to  red-black system if requested
c
      nb = iparm(9)
      if (nb.ge.0) go to 110
      if (nb.le.-2) go to 170
      n3 = n*3
      do 90 i = 1,n3
         iwksp(i) = 0
   90 continue
      call prbndx (n,ndim,maxnz,jcoef,iwksp,iwksp(ib2),nb,level,nout,
     *   ier)
      if (ier.eq.0) go to 110
      if (level.ge.0) write (AMPI_LUN(nout),100) ier,nb
  100 format (/1x,'*** f a t a l     e r r o r ************'//1x,
     *   '    called from itpackv routine jcg  '/1x,
     *   '    error detected in routine  prbndx'/1x,
     *   '    which computes the red-black indexing'/1x,
     *   '    ier = ',i5,' iparm(9) = ',i5,' (nb)')
      go to 350
  110 if (nb.ge.0.and.nb.le.n) go to 130
      ier = 14
      if (level.ge.0) write (AMPI_LUN(nout),120) ier,nb
  120 format (/1x,'*** f a t a l     e r r o r ************'//1x,
     *   '    called from itpackv routine jcg      '/1x,
     *   '    error detected in red-black subsystem index'/1x,
     *   '    ier = ',i5,' iparm(9) =',i5,' (nb)')
      go to 350
  130 if (nb.ne.0.and.nb.ne.n) go to 150
      nbo = nb
      nb = n/2
      if (level.ge.2) write (AMPI_LUN(nout),140) nbo,nb
  140 format (/10x,' nb = ',i5,' implies matrix is diagonal'/10x,
     *   ' nb reset to ',i5)
c
c ... permute matrix and rhs
c
  150 if (level.ge.2) write (AMPI_LUN(nout),160) nb
  160 format (/10x,'order of black subsystem = ',i5,' (nb)')
      if (iparm(9).ge.0) go to 170
      call permat (n,ndim,maxnz,jcoef,coef,iwksp,wksp,iwksp(ib3))
      call pervec (n,iwksp,rhs,wksp)
      call pervec (n,iwksp,u,wksp)
c
c ... check for sufficient workspace.
c
  170 iparm(8) = 4*n+4*itmax
      if (nw.ge.iparm(8)) go to 190
      ier = 12
      if (level.ge.0) write (AMPI_LUN(nout),180) nw,iparm(8)
  180 format (/1x,'*** f a t a l     e r r o r ************'//1x,
     *   '    called from itpackv routine jcg '/1x,
     *   '    not enough workspace at ',i10/1x,
     *   '    set iparm(8) =',i10,' (nw)')
      go to 330
c
  190 continue
      if (level.le.2) go to 220
      write (AMPI_LUN(nout),200)
  200 format (///1x,'in the following, rho and gamma are',
     *   ' acceleration parameters')
      if (adapt) write (AMPI_LUN(nout),210)
  210 format (1x,'cme is the estimate of the largest eigenvalue of',
     *   ' the jacobi matrix')
  220 continue
c     if (iparm(11).eq.0) timi1 = timer(0.0)
c
c ... compute initial pseudo-residual
c
      do 230 i = 1,nw
         wksp(i) = 0.0d0 !jgf46.00 added d0
  230 continue
      call scopy (n,rhs,1,wksp(ib2),1)
      call pjac (n,ndim,maxnz,jcoef,coef,u,wksp(ib2))
      do 240 i = 1,n
         wksp(n+i) = wksp(n+i)-u(i)
  240 continue
c
c ... iteration sequence
c
      itmax1 = itmax+1
      do 260 loop = 1,itmax1
         in = loop-1
         if (mod(in,2).eq.1) go to 250
c
c ... code for the even iterations.
c
c     u           = u(in)             wksp(ib2) = del(in)
c     wksp(ib1)   = u(in-1)           wksp(ib3) = del(in-1)
c
         call itjcg (n,ndim,maxnz,jcoef,coef,u,wksp(ib1),wksp(ib2),
     *      wksp(ib3),wksp(ib4),wksp(ib5))
c
         if (halt) go to 290
         go to 260
c
c ... code for the odd iterations.
c
c     u           = u(in-1)           wksp(ib2) = del(in-1)
c     wksp(ib1)   = u(in)             wksp(ib3) = del(in)
c
  250    call itjcg (n,ndim,maxnz,jcoef,coef,wksp(ib1),u,wksp(ib3),
     *      wksp(ib2),wksp(ib4),wksp(ib5))
c
         if (halt) go to 290
  260 continue
c
c..... itmax has been reached
c
      if (iparm(11).ne.0) go to 270
c     timi2 = timer(0.0)
c     time1 = timi2-timi1
  270 ier = 13
      if (level.ge.1) write (AMPI_LUN(nout),280) itmax
  280 format (/1x,'*** w a r n i n g ************'//1x,
     *   '    in itpackv routine jcg'/1x,
     *   '    failure to converge in',i5,' iterations')
      if (iparm(3).eq.0) rparm(1) = stptst
      go to 320
c
c ... method has converged
c
  290 if (iparm(11).ne.0) go to 300
c     timi2 = timer(0.0)
c     time1 = timi2-timi1
  300 if (level.ge.1) write (AMPI_LUN(nout),310) in
  310 format (/1x,'jcg  has converged in ',i5,' iterations')
c
c ... put solution into u if not already there.
c
  320 continue
      if (mod(in,2).eq.1) call scopy (n,wksp,1,u,1)
c
c ... un-permute matrix,rhs, and solution
c
  330 if (iparm(9).ne.-1) go to 340
      call permat (n,ndim,maxnz,jcoef,coef,iwksp(ib2),wksp(ib4),
     *   iwksp(ib3))
      call pervec (n,iwksp(ib2),rhs,wksp(ib4))
      call pervec (n,iwksp(ib2),u,wksp(ib4))
      if (ier.eq.12) go to 350
c
c ... optional error analysis
c
  340 idgts = iparm(12)
      if (idgts.lt.0) go to 350
      if (iparm(2).le.0) idgts = 0
      call peror (n,ndim,maxnz,jcoef,coef,rhs,u,wksp,digit1,digit2,
     *   idgts)
c
c ... unscale the matrix, solution, and rhs vectors.
c
  350 continue
      call unscal (n,ndim,maxnz,jcoef,coef,rhs,u,wksp)
c
c ... set return parameters in iparm and rparm
c
      iparm(8) = iparm(8)-4*(itmax-in)
C      if (iparm(11).ne.0) go to 360 !jgf46.00 commented out
c     timj2 = timer(0.0)
c     time2 = timj2-timj1
C      time2 = 0.0                   !jgf46.00 commented out
  360 if (iparm(3).ne.0) go to 370
      iparm(1) = in
      iparm(9) = nb
      rparm(2) = cme
      rparm(3) = sme
C      rparm(9) = time1              !jgf46.00 commented out
C      rparm(10) = time2             !jgf46.00 commented out
      rparm(11) = digit1
      rparm(12) = digit2
c
  370 continue
      if (level.ge.3) call echall (n,ndim,maxnz,jcoef,coef,rhs,iparm,
     *   rparm,2)
      if (ier.eq.0.and.level.ge.1) write (AMPI_LUN(nout),380)
  380 format (/1x,'execution successful')
c
      return
      end subroutine


      subroutine itjcg (n,ndim,maxnz,jcoef,coef,u,u1,d,d1,dtwd,tri)
      use global,only:nperseg
      implicit none
C
      integer n,ndim,maxnz,i,j
      integer jcoef(ndim,maxnz)
      real(sz) coef(ndim,maxnz)
      real(sz) u(n),u1(n),d(n),d1(n),dtwd(n),tri(*)
      real(sz) gamold,rhoold,rhotmp,dnrm,con,dtnrm,c1,c2,c3,c4
      real(sz) dumy1(1),dumy2(1),del3nrms(3),unorm
      real(sz) vic1, vic2 !jgf46.00 added
      logical q1
c
c *** begin -- itpackv common
c
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
c
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
c
      real(sz) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
c *** end   -- itpackv common
c
c ... itjcg performs one iteration of the jacobi conjugate gradient
c     algorithm.  it is called by jcg.
c
c ... parameter list --
c
c          n      input integer.  dimension of the matrix.
c          ndim   row dimension of jcoef and coef arrays in calling
c                   routine
c          maxnz  maximum number of nonzeros per row
c          jcoef  integer sparse matrix representation
c          coef   sparse matrix representation
c          u      input vector.  contains the value of the
c                 solution vector at the end of in iterations.
c          u1     input/output vector.  on input, it contains
c                 the value of the solution at the end of the in-1
c                 iteration.  on output, it will contain the newest
c                 estimate for the solution vector.
c          d      input vector.  contains the pseudo-residual
c                 vector after in iterations.
c          d1     input/output vector.  on input, d1 contains
c                 the pseudo-residual vector after in-1 iterations.  on
c                 output, it will contain the newest pseudo-residual
c                 vector.
c          dtwd   array.  used in the computations of the
c                 acceleration parameter gamma and the new pseudo-
c                 residual.
c          tri    array.  stores the tridiagonal matrix associated
c                 with the eigenvalues of the conjugate gradient
c                 polynomial.
c
c ... local itpackv references --
c          chgcon, iterm , parcon, pjac  , pstop 
c
c     description of variables in common blocks in routine jcg
c
c ... compute new estimate for cme if adapt = .true.
c
      save

      if (adapt) call chgcon (itmax,tri,gamold,rhoold,1)
c
C   ...UPDATE PSEUDO-RESIDUAL VECTOR "d" OF SIZE "n"
C   ...BEFORE PERFORMING ONE JACOBI ITERATION
C

#ifdef CMPI
      dumy1(1) = 0.0d0 !jgf46.00 added
      dumy2(1) = 0.0d0 !jgf46.00 added
      if(nperseg.gt.0)then
          CALL UPDATER_W_PERBC(d,dumy1,dumy2,1) ;
      else
          CALL UPDATER(d,dumy1,dumy2,1)
      endif
#endif
 
      do 10 i = 1,n
         dtwd(i) = 0.0d0
   10 continue
      call pjac (n,ndim,maxnz,jcoef,coef,d,dtwd)

C       
c ... test for stopping
c
       if (q1 .OR.  ((in .gt.5) .AND. 
     1        (mod(in,5).ne. 0))) then

#ifdef CMPI
        call ps2dots(n,d,dtwd,del3nrms)
        delnnm = del3nrms(1)
        dtnrm  = del3nrms(2)
        unorm  = del3nrms(3)
#else
      delnnm = sdot(n,d,1,d,1)
      dtnrm = sdot(n,d,1,dtwd,1)
      unorm = 1.0d0 !jgf46.00 added
#endif

        else

#ifdef CMPI
        call ps3dots(n,d,dtwd,u,del3nrms)
        delnnm = del3nrms(1)
        dtnrm  = del3nrms(2)
        unorm  = del3nrms(3)
#else
      delnnm = sdot(n,d,1,d,1)
      dtnrm = sdot(n,d,1,dtwd,1)
      unorm = sdot(n,u,1,u,1)
#endif
       end if

      dnrm = delnnm
      con = cme
      call pstop_nrms (n,unorm,dnrm,con,1,q1)
      if (halt) go to 50
 


       ! ... compute rho and gamma - acceleration parameters
 

   20 call parcon (dtnrm,c1,c2,c3,c4,gamold,rhoold,1)
c
c ... compute u(in+1) and d(in+1)
c
   30 do 40 i = 1,n
         u1(i) = c1*d(i)+c2*u(i)+c3*u1(i)
         d1(i) = c1*dtwd(i)+c4*d(i)+c3*d1(i)
   40 continue

c
c ... output intermediate information
c
   50 call iterm (n,coef,u,dtwd,1)
c
      return
      end subroutine


      subroutine peror (n,ndim,maxnz,jcoef,coef,rhs,u,work,
     *                   digit1,digit2,idgts)
      implicit none
      integer i,n,ndim,maxnz,idgts
      integer jcoef(ndim,maxnz)
      real(sz)  digit1,digit2,bnrm,rnrm,temp
      real(sz) coef(ndim,maxnz),rhs(n),u(n),work(n)
c
c *** begin -- itpackv common
c
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
c
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
c
      real(sz) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
c *** end   -- itpackv common
c
c     peror computes the residual, r = rhs - a*u.  the user
c     also has the option of printing the residual and/or the
c     unknown vector depending on idgts.
c
c ... parameter list --
c
c          n      dimension of matrix
c          ndim   row dimension of jcoef and coef in calling routine
c          maxnz  maximum number of nonzeros per row
c          jcoef  integer array of sparse matrix representation
c          coef   array of sparse matrix representation
c          rhs    right hand side of matrix problem
c          u      latest estimate of solution
c          work   workspace vector of length 2*n
c          digit1 output - measure of accuracy of stopping test
c          digit2 output - measure of accuracy of solution
c          idgts   parameter controlling level of output
c                    if idgts < 1 or idgts > 4, then no output.
c                            = 1, then number of digits is printed, pro-
c                                 vided level .ge. 1
c                            = 2, then solution vector is printed, pro-
c                                 vided level .ge. 1
c                            = 3, then residual vector is printed, pro-
c                                 vided level .ge. 1
c                            = 4, then both vectors are printed, pro-
c                                 vided level .ge. 1
c
c ... local itpackv references --
c
c          pmult , vout
c
c ... specifications for arguments
c
c
c     description of variables in common block in main routine
c
      digit1 = 0.0
      digit2 = 0.0
      if (n.le.0) go to 70
c
      digit1 = -log10(abs(srelpr))
      if (stptst.gt.0.0) digit1 = -log10(abs(stptst))
      do 10 i = 1,n
         work(i) = rhs(i)/coef(i,1)
   10 continue

#ifdef CMPI
      bnrm = psdot(n,work,work)
#else
      bnrm = sdot(n,work,1,work,1)
#endif

      if (bnrm.eq.0.0) go to 30
      call pmult (n,ndim,maxnz,jcoef,coef,u,work)
      do 20 i = 1,n
         work(i) = (rhs(i)-work(i))/coef(i,1)
   20 continue

#ifdef CMPI
      rnrm = psdot(n,work,work)
#else
      rnrm = sdot(n,work,1,work,1)
#endif


      temp = rnrm/bnrm
      if (temp.eq.0.0) go to 30
      digit2 = -log10(abs(temp))/2.0d0
      go to 40
c
   30 digit2 = -log10(abs(srelpr))
c
   40 if ((idgts.lt.1).or.(level.le.0)) go to 70
      write (AMPI_LUN(nout),50) digit1,digit2
   50 format (/10x,'approx. no. of digits in stopping test =',
     *                f5.1,'  (digit1)'
     *        /10x,'approx. no. of digits in ratio test    =',
     *                f5.1,'  (digit2)')
c
      if (idgts.le.1.or.idgts.gt.4) go to 70
      if (idgts.ge.3) call vout (n,work,1,nout)
      do 60 i = 1,n
         work(i) = u(i)*coef(i,1)
   60 continue
      if (idgts.ne.3) call vout (n,work,2,nout)
c
   70 continue
      return
      end subroutine

      subroutine pstop (n,u,dnrm,ccon,iflag,q1)
      implicit none
      integer n,iflag
      logical q1
      real(sz) u(n),dnrm,con,ccon,uold,tr,tl
c
c *** begin -- itpackv common
c
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
c
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
c
      real(sz) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
c *** end   -- itpackv common
c
c     pstop performs a test to see if the iterative
c     method has converged to a solution inside the error
c     tolerance, zeta.
c
c ... parameter list --
c
c          n      order of system
c          u      present solution estimate
c          dnrm   inner product of pseudo-residuals at preceding
c                    iteration
c          con    stopping test parameter (= ccon)
c          iflag  stopping test integer flag
c                    iflag = 0,  sor iteration zero
c                    iflag = 1,  non-rs method
c                    iflag = 2,  rs method
c          q1     stopping test logical flag
c
c
c     description of variables in common block in main routine
c
      con = ccon
      halt = .false.
c
c     special procedure for zeroth iteration
c
      if (in.ge.1) go to 10
      q1 = .false.
      udnm = 1.0d0
      stptst = 1000.0d0
      if (iflag.le.0) return
c
c ... test if udnm needs to be recomputed
c
   10 continue
      if (q1) go to 20
      if ((in.gt.5).and.(mod(in,5).ne.0)) go to 20
      uold = udnm

#ifdef CMPI
      udnm = psdot(n,u,u)
#else
      udnm = sdot(n,u,1,u,1)
#endif

      if (udnm.eq.0.0) udnm = 1.0d0
      if ((in.gt.5).and.(abs(udnm-uold).le.udnm*zeta)) q1 = .true.
c
c ... compute stopping test
c
   20 tr = sqrt(udnm)
      tl = 1.0d0
      if (con.eq.1.0d0) go to 40
      if (iflag.eq.2) go to 30
      tl = sqrt(dnrm)
      tr = tr*(1.0d0-con)
      go to 40
   30 tl = sqrt(2.0d0*dnrm)
      tr = tr*(1.0d0-con*con)
   40 stptst = tl/tr
      if (tl.ge.tr*zeta) return
      halt = .true.
c
      return
      end subroutine


      subroutine pstop_nrms (n,unrm,dnrm,ccon,iflag,q1)
      implicit none
      integer n,iflag
      logical q1
      real(sz) unrm,dnrm,con,ccon,uold,tr,tl
c
c *** begin -- itpackv common
c
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
c
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
c
      real(sz) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
c *** end   -- itpackv common
c
c     pstop performs a test to see if the iterative
c     method has converged to a solution inside the error
c     tolerance, zeta.
c
c ... parameter list --
c
c          n      order of system
c          u      present solution estimate
c          dnrm   inner product of pseudo-residuals at preceding
c                    iteration
c          con    stopping test parameter (= ccon)
c          iflag  stopping test integer flag
c                    iflag = 0,  sor iteration zero
c                    iflag = 1,  non-rs method
c                    iflag = 2,  rs method
c          q1     stopping test logical flag
c
c
c     description of variables in common block in main routine
c
      con = ccon
      halt = .false.
c
c     special procedure for zeroth iteration
c
      if (in.ge.1) go to 10
      q1 = .false.
      udnm = 1.0d0
      stptst = 1000.0d0
      if (iflag.le.0) return
c
c ... test if udnm needs to be recomputed
c
   10 continue
      if (q1) go to 20
      if ((in.gt.5).and.(mod(in,5).ne.0)) go to 20
      uold = udnm

      udnm = unrm

      if (udnm.eq.0.0d0) udnm = 1.0d0 !jgf46.00 added d0 to 0.0
      if ((in.gt.5).and.(abs(udnm-uold).le.udnm*zeta)) q1 = .true.
c
c ... compute stopping test
c
   20 tr = sqrt(udnm)
      tl = 1.0d0
      if (con.eq.1.0d0) go to 40
      if (iflag.eq.2) go to 30
      tl = sqrt(dnrm)
      tr = tr*(1.0d0-con)
      go to 40
   30 tl = sqrt(2.0d0*dnrm)
      tr = tr*(1.0d0-con*con)
   40 stptst = tl/tr
      if (tl.ge.tr*zeta) return
      halt = .true.
c
      return
      end subroutine


      subroutine chgcon (ldt,tri,gamold,rhoold,ibmth)
      implicit none
      integer ldt,ibmth,ip,ier
      real(sz) tri(ldt,4)
      real(sz) gamold,rhoold,cmold,start,end1
c
c *** begin -- itpackv common
c
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
c
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
c
      real(sz) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
c *** end   -- itpackv common
c
c ... chgcon computes the new estimate for the largest eigenvalue
c     for conjugate gradient acceleration.
c
c ... parameter list --
c
c          ldt    leading dimension of tri
c          tri    tridiagonal matrix associated with the eigenvalues
c                    of the conjugate gradient polynomial
c          gamold
c            and
c          rhoold previous values of acceleration parameters
c          ibmth  indicator of basic method being accelerated by cg
c                      ibmth = 1,  jacobi
c                            = 2,  reduced system
c                            = 3,  ssor
c
c ... local itpackv references --
c
c          eigvns, eigvss
      save
      go to (10,20,30), ibmth
c
c ... jacobi conjugate gradient
c
   10 start = cme
      ip = in
      go to 40
c
c ... reduced system cg
c
   20 start = cme**2
      ip = in
      go to 40
c
c ... ssor cg
c
   30 if (adapt) start = spr
      if (.not.adapt) start = specr
      ip = in-is
c
c ... define the matrix
c
   40 if (ip.ge.2) go to 60
      if (ip.eq.1) go to 50
c
c ... ip = 0
c
      end1 = 0.0d0   !jgf46.00 added d0
      cmold = 0.0d0  !jgf46.00 added d0
      go to 110
c
c ... ip = 1
c
   50 end1 = 1.0D0-1.0D0/gamma
      tri(1,1) = end1
      tri(1,2) = 0.0D0
      go to 110
c
c ... ip > 1
c
   60 if (abs(start-cmold).le.zeta*start) go to 120
      cmold = start
c
c ... compute the largest eigenvalue
c
      tri(ip,1) = 1.0d0-1.0d0/gamma
      tri(ip,2) = (1.0d0-rho)/(rho*rhoold*gamma*gamold)
      if (isym.ne.0) go to 80
      end1 = eigvss(ip,tri,start,zeta,itmax,ier)
      if (ier.eq.0) go to 100
      if (level.ge.2) write (AMPI_LUN(nout),70) ier
   70 format (/10x,'difficulty in computation of maximum eigenvalue'/
     *         15x,'of iteration matrix'/
     *         10x,'routine zbrent returned ier =',i5)
      go to 100
   80 continue
      end1 = eigvns(ldt,ip,tri,tri(1,3),tri(1,4),ier)
      if (ier.eq.0) go to 100
      if (level.ge.2) write (AMPI_LUN(nout),90) ier
   90 format (/10x,'difficulty in computation of maximum eigenvalue'/
     *         15x,'of iteration matrix'/
     *         10x,'routine eqrt1s returned ier =',i5)
  100 continue
      if (ier.ne.0) go to 130
c
c ... set spectral radius for the various methods
c
  110 if (ibmth.eq.1) cme = end1
      if (ibmth.eq.2) cme = sqrt(abs(end1))
      if (ibmth.eq.3.and.adapt) spr = end1
      if (ibmth.eq.3.and..not.adapt) specr = end1
      return
c
c ... relative change in cme is less than zeta.  therefore stop
c     changing.
c
  120 adapt = .false.
      partad = .false.
      return
c
c ... estimate for cme > one.  therefore need to stop adaptive
c     procedure and keep old value of cme.
c
  130 adapt = .false.
      partad = .false.
      if (level.ge.2) write (AMPI_LUN(nout),140) in,start
  140 format (/10x,'estimate of maximum eigenvalue of jacobi   '/15x,
     *   'matrix (cme) not accurate'/10x,
     *   'adaptive procedure turned off at iteration ',i5/10x,
     *   'final estimate of maximum eigenvalue =',e15.7/)
c
      return
      end  subroutine

      real(sz) function determ (ldt,n,tri,wk1,wk2,xlmda)
      implicit none
      integer ldt,n,i,l
      real(sz) xlmda
      real(sz) tri(ldt,2),wk1(n),wk2(n)
c
c     determ computes the determinant of a symmetric
c     tridiagonal matrix given by tri. det(tri - xlmda*i) = 0
c
c ... parameter list --
c
c          ldt    leading dimension of array tri
c          n      order of tridiagonal system
c          tri    symmetric tridiagonal matrix of order n
c          wk1,   workspace vectors of length n
c           wk2
c          xlmda  argument for characteristic equation
c
c
      do 10 i = 1,n
         wk1(i) = tri(i,1)-xlmda
   10 continue
      wk2(n) = wk1(n)
      wk2(n-1) = wk1(n-1)*wk2(n)+tri(n,2)
      if (n.eq.2) go to 30
c
c ... beginning of loop
c
      do 20 l = n-2,1,-1
         wk2(l) = wk1(l)*wk2(l+1)+tri(l+1,2)*wk2(l+2)
   20 continue
c
c     wk2(1) = solrn (n,wk1(-1),-1,tri(0,2),-1,wk2,-1)
c
c ... determinant computed
c
   30 determ = wk2(1)
c
      return
      end  function


      subroutine dfault (iparm,rparm)
      implicit none
      integer iparm(12)
      real(sz) rparm(12), temp
c
c *** begin -- itpackv common
c
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
c
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
c
      real(sz) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
c *** end   -- itpackv common
c
c ... dfault sets the default values of iparm and rparm.
c
c ... parameter list --
c
c          iparm
c           and
c          rparm  arrays specifying options and tolerances
c
c ... specifications for arguments
c
c
c     description of variables in common blocks in main routine
c
c     srelpr  - computer precision (approx.)
c     if installer of package does not know srelpr value,
c     an approximate value can be determined from a simple
c     fortran program such as
c
c     srelpr = 1.0d0
c   2 srelpr = 0.5d0*srelpr
c     temp = srelpr + 1.0d0
c     if (temp .gt. 1.0d0)  go to 2
c     srelpr = 2.0d0*srelpr
c     write (6,3) srelpr
c   3 format (5x,e15.8)
c     stop
c     end
c
c     some values are-
c
c     srelpr = 7.1e-15   for cray x-mp, y-mp  (approx.) 2**-47
c          = 1.49e-8   for dec 10  (approx.) 2**-26
c          = 1.192e-7  for vax 11/780 (approx) 2**-23
c          = 1.192e-7  for Sun Spark Station 2 (approx) 2**-23
c          = 1.192e-7  for IBM RISC 6000 2 (approx) 2**-23
c          = 4.768e-7  for ibm 370/158
C          = 0.5960E-7 for Lahey fortran on ALR (486 PC) rl 12/92
c          = 0.5960E-7 for Vax 6000
c             *** should be changed for other machines ***
c
c     to facilitate convergence, rparm(1) should be set to
c          500.*srelpr or larger
c
cvjp--Determine macheps 
c
      srelpr = 1.0d0
    2 srelpr = 0.5d0*srelpr
      temp = srelpr + 1.0D0
      if (temp .gt. 1.0d0)  go to 2
      srelpr = 2.0d0*srelpr
c
      iparm(1) = 100
      iparm(2) = -1 !jgf46.00 changed to -1, was 0
      iparm(3) = 0
      iparm(4) = 6
      iparm(5) = 0
      iparm(6) = 1
      iparm(7) = 1
      iparm(8) = 0
      iparm(9) = -2
      iparm(10) = 0
      iparm(11) = 0
      iparm(12) = 3 !jgf46.00 changed to -3, was 0
c
      rparm(1) = 5.0d-6 !jgf46.00 was 512.d0*srelpr
      rparm(2) = 0.0d0
      rparm(3) = 0.0d0
      rparm(4) = .75d0
      rparm(5) = 1.0d0
      rparm(6) = 0.0d0
      rparm(7) = .25d0
      rparm(8) = 128.0d0*srelpr !jgf46.00 was 100.0d0*srelpr
      rparm(9) = 0.0d0
      rparm(10) = 0.0d0
      rparm(11) = 0.0d0
      rparm(12) = 0.0d0
c
      return
      end  subroutine

      subroutine echall (n,ndim,maxnz,jcoef,coef,rhs,iparm,rparm,icall)
      implicit none
      integer n,i,j,ndim,maxnz,icall
      integer jcoef(ndim,maxnz),iparm(12)
      real(sz) coef(ndim,maxnz),rhs(n),rparm(12)
c
c *** begin -- itpackv common
c
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
c
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
c
      real(sz) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
c *** end   -- itpackv common
c
c ... echall initializes the itpackv common blocks from the
c ... information contained in iparm and rparm. echall also prints the
c ... values of all the parameters in iparm and rparm.
c
c ... parameter list --
c
c          iparm
c           and
c          rparm  arrays of parameters specifying options and
c                    tolerances
c          icall  indicator of which parameters are being printed
c                    icall = 1,  initial parameters
c                    icall = 2,  final parameters
c
c ... specifications for arguments
c
c     description of variables in common blocks in main routine
c
      if (icall.ne.1) go to 120
c
c ... initialize itpackv common
c
      zeta = rparm(1)
      cme = rparm(2)
      sme = rparm(3)
      ff = rparm(4)
      omega = rparm(5)
      specr = rparm(6)
      betab = rparm(7)
      itmax = iparm(1)
      level = iparm(2)
      isym = iparm(5)
c
      adapt = .false.
      partad = .false.
      betadt = .false.
      if (iparm(6).eq.1.or.iparm(6).eq.3) adapt = .true.
      if (iparm(6).eq.1) betadt = .true.
      if (iparm(6).eq.2) partad = .true.
c
      caseii = .false.
      if (iparm(7).eq.2) caseii = .true.
      if (caseii) sme = -cme
      if (.not.caseii.and.sme.eq.0.0d0) sme = -1.0d0 !jgf46.00 added d0
      spr = sme
c
c ... set rest of common variables to zero
c
      in = 0
      is = 0
      halt = .false.
      bdelnm = 0.0d0 !jgf46.00 added d0 to all these
      delnnm = 0.0d0
      delsnm = 0.0d0
      gamma = 0.0d0
      qa = 0.0d0
      qt = 0.0d0
      rho = 0.0d0
      rrr = 0.0d0
      sige = 0.0d0
      stptst = 0.0d0
      udnm = 0.0d0
c
      if (level.le.4) go to 100
c
c     this section of echall causes printing of the linear system and
c     the iterative parameters
c
      write (AMPI_LUN(nout),10)
   10 format (///5x,'the linear system is as follows')
      write (AMPI_LUN(nout),20)
   20 format (/2x,'jcoef array')
      do 30 i = 1,n
         write (AMPI_LUN(nout),40) (jcoef(i,j),j=1,maxnz)
   30 continue
   40 format (1x,8(1x,i8))
      write (AMPI_LUN(nout),50)
   50 format (/2x,'coef array')
      do 60 i = 1,n
         write (AMPI_LUN(nout),70) (coef(i,j),j=1,maxnz)
   60 continue
   70 format (1x,5(2x,g14.6))
      write (AMPI_LUN(nout),80)
   80 format (/2x,'rhs array')
      write (AMPI_LUN(nout),90) (rhs(i),i=1,n)
   90 format (1x,5g16.6)
  100 if (level.le.2) return
      write (AMPI_LUN(nout),110)
  110 format (///5x,'initial iterative parameters'/)
      go to 140
  120 write (AMPI_LUN(nout),130)
  130 format (///5x,'final iterative parameters'/)
  140 write (AMPI_LUN(nout),150) iparm(1),level,iparm(3),nout,isym,iparm(6)
  150 format (10x,'iparm(1)  =',i15,4x,'(itmax)'/
     *        10x,'iparm(2)  =',i15,4x,'(level)'/
     *        10x,'iparm(3)  =',i15,4x,'(ireset)'/
     *        10x,'iparm(4)  =',i15,4x,'(nout)'/
     *        10x,'iparm(5)  =',i15,4x,'(isym)'/
     *        10x,'iparm(6)  =',i15,4x,'(iadapt)')
      write (AMPI_LUN(nout),160) iparm(7),iparm(8),iparm(9),iparm(10),iparm(11),
     *  iparm(12)
  160 format (10x,'iparm(7)  =',i15,4x,'(icase)'/
     *        10x,'iparm(8)  =',i15,4x,'(nwksp)'/
     *        10x,'iparm(9)  =',i15,4x,'(nb)'/
     *        10x,'iparm(10) =',i15,4x,'(iremove)'/
     *        10x,'iparm(11) =',i15,4x,'(itime)'/
     *        10x,'iparm(12) =',i15,4x,'(idgts)')
      write (AMPI_LUN(nout),170) zeta,cme,sme,ff,omega,specr
  170 format (10x,'rparm(1)  =',e15.8,4x,'(zeta)'/
     *        10x,'rparm(2)  =',e15.8,4x,'(cme)'/
     *        10x,'rparm(3)  =',e15.8,4x,'(sme)'/
     *        10x,'rparm(4)  =',e15.8,4x,'(ff)'/
     *        10x,'rparm(5)  =',e15.8,4x,'(omega)'/
     *        10x,'rparm(6)  =',e15.8,4x,'(specr)')
      write (AMPI_LUN(nout),180) betab,rparm(8),rparm(9),rparm(10),rparm(11),
     *  rparm(12)
  180 format (10x,'rparm(7)  =',e15.8,4x,'(betab)'/
     *        10x,'rparm(8)  =',e15.8,4x,'(tol)'/
     *        10x,'rparm(9)  =',e15.8,4x,'(time1)'/
     *        10x,'rparm(10) =',e15.8,4x,'(time2)'/
     *        10x,'rparm(11) =',e15.8,4x,'(digit1)'/
     *        10x,'rparm(12) =',e15.8,4x,'(digit2)')
c
      return
      end  subroutine


      real(sz) function eigvns (ldt,n,tri,d,e2,ier)
      implicit none
      integer ldt,n,i,ier
      real(sz) tri(ldt,*),d(n),e2(n)
c
c ... eigvns computes the largest eigenvalue of a symmetric
c     tridiagonal matrix for conjugate gradient acceleration.
c
c ... parameter list --
c
c          ldt    leading dimension of tri
c          n      order of tridiagonal system
c          tri    symmetric tridiagonal matrix of order n
c          d      array for eqrt1s (negative diagonal elements)
c          e2     array for eqrt1s (super diagonal elements)
c          ier    error flag -- on return, ier=0 indicates that
c                    the largest eigenvalue of tri was found.
c
c ... local itpackv references --
c
c          eqrt1s
c
c ... specifications for arguments
c
c
      eigvns = 0.0
c
      d(1) = -tri(1,1)
      do 10 i = 2,n
         d(i) = -tri(i,1)
         e2(i) = abs(tri(i,2))
   10 continue
c
      call eqrt1s (d,e2,n,1,0,ier)
      eigvns = -d(1)
c
      return
      end function 


      real(sz) function eigvss (n,tri,start,zeta,itmax,ier)
      implicit none
      integer n,itmax,itmp,ier,maxfn,nsig
      real(sz) tri(*),start,eps,a,b,zeta
c
c ... eigvss computes the largest eigenvalue of a symmetric
c     tridiagonal matrix for conjugate gradient acceleration.
c     modified imsl routine zbrent used.
c
c ... parameter list --
c
c          n      order of tridiagonal system
c          tri    symmetric tridiagonal matrix of order n
c          start  initial lower bound of interval containing root
c          zeta   stopping criteria
c          ier    error flag -- on return, ier = 0 indicates that
c                    the largest eigenvalue of tri was found.
c
c ... local itpackv references --
c
c          zbrent
c
c ... specifications for arguments
c
c
      eigvss = 0.0d0
C      itmp = int(-log10(abs(zeta))) !jgf46.00 commented out
cvjp 2/5/06 added kind parameter for ibm platform
#ifdef IBM
      itmp = int(-log10(abs(zeta)),KIND(0.0d0)) !jgf46.00 added
#else
      itmp = int(-log10(abs(zeta))) !jgf46.16 added
#endif
      nsig = max0(itmp,4)
      maxfn = max0(itmax,50)
      eps = 0.0d0
      a = start
      b = 1.0d0
      call zbrent (n,tri,eps,nsig,a,b,maxfn,ier)
      eigvss = b
c
      return
      end  function


      subroutine eqrt1s (d,e2,n,m,isw,ier)
      implicit none
      integer n,m,i,j,k,isw,ier,ii,k1
      real(sz) d(n),e2(n),err,s,tot,p,q,qp,r,delta,f,ep,dlam
c
c *** begin -- itpackv common
c
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
c
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
c
      real(sz) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
c *** end   -- itpackv common
c
c   modified imsl routine name   - eqrt1s
c
c-----------------------------------------------------------------------
c
c   computer            - cdc/single
c
c   latest revision     - june 1, 1980
c
c   purpose             - smallest or largest m eigenvalues of a
c                           symmetric tridiagonal matrix
c
c   usage               - call eqrt1s (d,e2,n,m,isw,ier)
c
c   arguments    d      - input vector of length n containing
c                           the diagonal elements of the matrix.  the
c                           computed eigenvalues replace the first m
c                           components of the vector d in non-
c                           decreasing sequence, while the remaining
c                           components are lost.
c                e2     - input vector of length n containing
c                           the squares of the off-diagonal elements
c                           of the matrix.  input e2 is destroyed.
c                n      - input scalar containing the order of the
c                           matrix.
c                m      - input scalar containing the number of
c                           smallest eigenvalues desired (m is
c                           less than or equal to n).
c                isw    - input scalar meaning as follows -
c                           isw=1 means that the matrix is known to be
c                             positive definite.
c                           isw=0 means that the matrix is not known
c                             to be positive definite.
c                ier    - error parameter. (output)
c                           warning error
c                             ier = 601 indicates that successive
c                               iterates to the k-th eigenvalue were not
c                               monotone increasing. the value k is
c                               stored in e2(1).
c                           terminal error
c                             ier = 602 indicates that isw=1 but matrix
c                               is not positive definite
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c   remarks      as written, the routine computes the m smallest
c                eigenvalues. to compute the m largest eigenvalues,
c                reverse the sign of each element of d before and
c                after calling the routine. in this case, isw must
c                equal zero.
c
c   copyright           - 1980 by imsl, inc. all rights reserved.
c
c   warranty            - imsl warrants only that imsl testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c-----------------------------------------------------------------------
c
c
c                                  srelpr = machine precision
c                                  first executable statement
c
      ier = 0
      dlam = 0.0
      err = 0.0
      s = 0.0
c
c                                  look for small sub-diagonal entries
c                                  define initial shift from lower
c                                  gerschgorin bound.
c
      tot = d(1)
      q = 0.0d0
      j = 0
      do 30 i = 1,n
         p = q
         if (i.eq.1) go to 10
         if (p.gt.srelpr*(abs(d(i))+abs(d(i-1)))) go to 20
   10    e2(i) = 0.0d0
c
c                                  count if e2(i) has underflowed
c
   20    if (e2(i).eq.0.0d0) j = j+1
         q = 0.0d0
         if (i.ne.n) q = sqrt(abs(e2(i+1)))
         tot = min(d(i)-p-q,tot)
   30 continue
      if (isw.eq.1.and.tot.lt.0.0) go to 50
      do 40 i = 1,n
         d(i) = d(i)-tot
   40 continue
      go to 60
   50 tot = 0.0d0
   60 do 190 k = 1,m
c
c                                  next qr transformation
c
   70    tot = tot+s
         delta = d(n)-s
         i = n
         f = abs(srelpr*tot)
         if (dlam.lt.f) dlam = f
         if (delta.gt.dlam) go to 90
         if (delta.ge.(-dlam)) go to 160
         ier = 602
         if (level.ge.1) write (AMPI_LUN(nout),80)
   80    format (/1x,'*** w a r n i n g ************'/1x,
     *      '    in itpackv routine eqrt1s'/1x,
     *      '    parameter isw = 1 but matrix'/1x,
     *      '    not positive definite')
         go to 200
c
c                                  replace small sub-diagonal squares
c                                  by zero to reduce the incidence of
c                                  underflows
c
   90    if (k.eq.n) go to 110
         k1 = k+1
         do 100 j = k1,n
            if (e2(j).le.(srelpr*(d(j)+d(j-1)))**2) e2(j) = 0.0
  100    continue
  110    f = e2(n)/delta
         qp = delta+f
         p = 1.0
         if (k.eq.n) go to 140
         k1 = n-k
         do 130 ii = 1,k1
            i = n-ii
            q = d(i)-s-f
            r = q/qp
            p = p*r+1.0d0
            ep = f*r
            d(i+1) = qp+ep
            delta = q-ep
            if (delta.gt.dlam) go to 120
            if (delta.ge.(-dlam)) go to 160
            ier = 602
            if (level.ge.1) write (AMPI_LUN(nout),80)
            go to 200
  120       f = e2(i)/q
            qp = delta+f
            e2(i+1) = qp*ep
  130    continue
  140    d(k) = qp
         s = qp/p
         if (tot+s.gt.tot) go to 70
         ier = 601
         e2(1) = k
         if (level.ge.1) write (AMPI_LUN(nout),150) k
  150    format (/1x,'*** w a r n i n g ************'//1x,
     *      '    in itpackv routine eqrt1s'/1x,
     *      '    successive iterates to the',i10/1x,
     *      '    eigenvalue were not monotone increasing')
c
c                                  set error -- irregular end
c                                  deflate minimum diagonal element
c
         s = 0.0
         i = ismin(n-k+1,d(k),1)
         delta = min(qp,d(i))
c
c                                  convergence
c
  160    if (i.lt.n) e2(i+1) = e2(i)*f/qp
         if (i.eq.k) go to 180
         do 170 j = i-1,k,-1
            d(j+1) = d(j)-s
            e2(j+1) = e2(j)
  170    continue
  180    d(k) = tot
         err = err+abs(delta)
         e2(k) = err
  190 continue
      if (ier.eq.0) go to 210
  200 continue
  210 return
      end  subroutine

      subroutine iterm (n,coef,u,wk,imthd)
      implicit none
      integer n,i,ip,imthd
      real(sz) coef(n,*),u(n),wk(n),qtff
c
c *** begin -- itpackv common
c
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
c
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
c
      real(sz) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
c *** end   -- itpackv common
c
c     iterm produces the iteration summary line at the end
c     of each iteration. if level .ge. 4, the latest approximation
c     to the solution will be printed.
c
c ... parameter list --
c
c          n      order of system or, for reduced system
c                    routines, order of black subsystem
c          coef   iteration matrix
c          u      solution estimate
c          wk     work array of length n
c          imthd  indicator of method
c                    imthd = 1,  jcg
c                    imthd = 2,  jsi
c                    imthd = 3,  sor
c                    imthd = 4,  ssorcg
c                    imthd = 5,  ssorsi
c                    imthd = 6,  rscg
c                    imthd = 7,  rssi
c
c ... specifications for arguments
c
c
c ... print various parameters after each iteration
c
      if (level.lt.2) return
      go to (10,100,140,170,50,10,100), imthd
   10 if (in.gt.0) go to 30
c
c ... print header for jcg and rscg
c
      write (AMPI_LUN(nout),20)
   20 format (////5x,'intermediate output after each iteration'//
     *   ' number of',3x,'convergence',5x,'cme ',10x,'rho',8x,'gamma'/
     *   ' iterations',5x,'test '//)
c
c ... print summary line
c
   30 write (AMPI_LUN(nout),40) in,stptst,cme,rho,gamma
   40 format (3x,i5,3x,5e13.5)
      if (level.ge.4) go to 200
c
      return
c
   50 if (in.gt.0) go to 70
c
c ... print header for ssor-si
c
      write (AMPI_LUN(nout),60)
   60 format (////5x,'intermediate output after each iteration'//
     *   ' number of',3x,'convergence',3x,'parameter change test',8x,
     *   'rho',8x,'gamma'/' iterations',5x,'test ',6x,'lhs(qa)',4x,
     *   'rhs(qt**ff)'//)
c
c ... print summary line
c
   70 ip = in-is
      if (imthd.eq.7) ip = 2*ip
      if (ip.lt.3) go to 80
      qtff = qt**ff
      write (AMPI_LUN(nout),40) in,stptst,qa,qtff,rho,gamma
      if (level.ge.4) go to 200
      return
c
   80 write (AMPI_LUN(nout),90) in,stptst,rho,gamma
   90 format (3x,i5,3x,e13.5,26x,2e13.5)
      if (level.ge.4) go to 200
      return
c
  100 if (in.gt.0) go to 120
c
c ... print header for j-si and rs-si
c
      write (AMPI_LUN(nout),110)
  110 format (////5x,'intermediate output after each iteration'//
     *   ' number of',3x,'convergence',3x,'parameter change test',8x,
     *   'rho'/' iterations',5x,'test ',6x,'lhs(qa)',4x,'rhs(qt**ff)'//)
c
c ... print summary line
c
  120 ip = in-is
      if (imthd.eq.7) ip = 2*ip
      if (ip.lt.3) go to 130
      qtff = qt**ff
      write (AMPI_LUN(nout),40) in,stptst,qa,qtff,rho
      if (level.ge.4) go to 200
      return
c
  130 write (AMPI_LUN(nout),90) in,stptst,rho
      if (level.ge.4) go to 200
      return
c
c ... print various parameters after each iteration for sor.
c
  140 if (in.gt.0) go to 160
c
c ... print header for sor
c
      write (AMPI_LUN(nout),150)
  150 format (////5x,'intermediate output after each iteration'//
     *   ' number of',3x,'convergence',5x,'cme ',8x,'omega',7x,
     *   'spectral'/' iterations',5x,'test',34x,'radius'//)
c
c ... print summary line for sor
c
  160 continue
      write (AMPI_LUN(nout),40) in,stptst,cme,omega,specr
      if (level.ge.4) go to 200
c
      return
c
c ... print various parameters after each iteration for ssor-cg.
c
  170 if (in.gt.0) go to 190
c
c ... print header for ssor-cg
c
      write (AMPI_LUN(nout),180)
  180 format (////5x,'intermediate output after each iteration'//
     *   ' number of',3x,'convergence',2x,' spectral',5x,'s-prime',9x,
     *   'rho',8x,'gamma'/' iterations',5x,'test ',7x,'radius'//)
c
c ... print summary line for ssor-cg
c
  190 continue
      write (AMPI_LUN(nout),40) in,stptst,specr,spr,rho,gamma
      if (level.ge.4) go to 200
      return
c
  200 if (imthd.gt.5) go to 220
      write (AMPI_LUN(nout),210) in
  210 format (/1x,2x,'estimate of solution at iteration ',i5)
      go to 240
  220 write (AMPI_LUN(nout),230) in
  230 format (/1x,2x,'estimate of solution at black points ',
     *   'at iteration ',i5)
  240 do 250 i = 1,n
         wk(i) = u(i)*coef(i,1)
  250 continue
      write (AMPI_LUN(nout),260) (wk(i),i=1,n)
  260 format (1x,5g16.7)
      write (AMPI_LUN(nout),270)
  270 format (//)
c
      return
      end  subroutine

      subroutine parcon (dtnrm,c1,c2,c3,c4,gamold,rhotmp,ibmth)
      implicit none
      integer ip,ibmth
      real(sz) dtnrm,c1,c2,c3,c4,gamold,rhotmp,rhoold
c
c *** begin -- itpackv common
c
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
c
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
c
      real(sz) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
c *** end   -- itpackv common
c
c ... parcon computes acceleration parameters for conjugate gradient
c     acceleration methods.
c
c ... parameter list --
c
c          dtnrm  inner product of residuals
c          c1     output -- rho*gamma
c          c2     output -- rho
c          c3     output -- 1-rho
c          c4     output -- rho*(1-gamma)
c          gamold output -- value of gamma at preceding iteration
c          rhotmp last estimate for value of rho
c          ibmth  indicator of basic method being accelerated by cg
c                      ibmth = 1,   jacobi
c                            = 2,   reduced system
c                            = 3,   ssor
c
c     description of variables in common blocks in main routine
c
      ip = in-is
c
c ... set rhoold and gamold
c
      rhoold = rho
      gamold = gamma
c
c ... compute gamma (in+1) 
c
c ... for jacobi or reduced system cg
c
C      if (ibmth.le.2) gamma = 1.0d0/(1.0d0-dtnrm/delnnm) !jgf46.00 comm. out
      gamma = 1.0d0/(1.0d0-dtnrm/delnnm) !jgf46.00 added
c
c ... for ssor cg
c
C      if (ibmth.eq.3) gamma = delnnm/dtnrm !jgf46.00 commented out
c
c ... compute rho (in+1)
c
      rho = 1.0d0
      if (ip.eq.0) go to 20
C      if (isym.eq.0) go to 10                 !jgf46.00 commented out
C      rho = 1.0d0/(1.0d0-gamma*rhotmp/delsnm) !jgf46.00 commented out
C      go to 20                                !jgf46.00 commented out
   10 rho = 1.0d0/(1.0d0-gamma*delnnm/(gamold*delsnm*rhoold))
c
c ... compute constants c1, c2, c3, and c4
c
   20 delsnm = delnnm
      rhotmp = rhoold
      c1 = rho*gamma
      c2 = rho
      c3 = 1.0d0-rho
      c4 = rho*(1.0d0-gamma)
c
      return
      end  subroutine



      subroutine permat (n,ndim,maxnz,jcoef,coef,p,work,iwork)
      implicit none
c
c ... permat takes the sparse matrix representation
c     of the matrix stored in the arrays jcoef and coef and
c     permutes both rows and columns, overwriting the previous
c     structure.
c
c ... parameter list --
c
c          n         order of system
c          ndim      row dimension of arrays jcoef and coef in
c                       the calling routine
c          maxnz     maximum number of nonzero entries per row
c          jcoef     integer array for data
c          coef      array for data structure coefficients
c          p         permutation vector
c          work      workspace of length n
c          iwork     integer workspace of length n
c
c ... it is assumed that the i-th entry of the permutation vector
c     p indicates the row the i-th row gets mapped into.  (i.e.
c     if ( p(i) = j ) row i gets mapped into row j)
c
c     *** note ***  this routine is to be called after routine scal.
c
c ... specifications for arguments
c
      integer ndim, maxnz, n, i, j
      integer jcoef(ndim,maxnz),p(*),iwork(n)
      real(sz) coef(ndim,maxnz),work(n)
c
      if (n.le.0) return
      do 50 j = 1,maxnz
         call scopy (n,coef(1,j),1,work,1)
         do 10 i = 1,n
            iwork(i) = jcoef(i,j)
   10    continue
         do 20 i = 1,n
            coef(p(i),j) = work(i)
   20    continue
         do 30 i = 1,n
            jcoef(p(i),j) = iwork(i)
   30    continue
         do 40 i = 1,n
            jcoef(i,j) = p(jcoef(i,j))
   40    continue
   50 continue
      return
      end  subroutine


      subroutine pervec (n,p,v,work)
      implicit none
      integer i,n
      integer p(n)
      real(sz) v(n),work(n)
c
c ... pervec permutes a vector as dictated by the permutation
c ... vector p.  if p(i) = j, then v(j) gets v(i).
c
c ... parameter list --
c
c          n       length of vectors p, v, and work
c          p       integer permutation vector
c          v       vector to be permuted
c          work    workspace vector of length n
c
c
      call scopy (n,v,1,work,1)
      do 10 i = 1,n
         v(p(i)) = work(i)
   10 continue
      return
      end  subroutine


      subroutine pjac (n,ndim,maxnz,jcoef,coef,u,rhs)
      implicit none
C      integer n,ndim,maxnz,maxm1 !jgf46.00 commented out
      integer i,j,n,ndim,maxnz    !jgf46.00 added
      integer jcoef(ndim,maxnz)
      real(sz) coef(ndim,maxnz),u(n),rhs(n)
c
c ... pjac performs one jacobi iteration.
c
c ... parameter list --
c
c         n       dimension of matrix
c         ndim    row dimension of jcoef and coef arrays in calling
c                   routine
c         maxnz   maximum number of nonzeros per row
c         jcoef   integer data structure for coefficient columns
c         coef    data structure for array coefficients
c         u       estimate of solution of a matrix problem
c         rhs     on input -- contains the right hand side of the
c                             matrix problem
c                 on output -- contains b*u + rhs  where b = i - a
c                              and a has been scaled to have a unit
c                              diagonal
c
c
C      maxm1 = maxnz-1  !jgf46.00 commented out
C      call ymasx2 (ndim,n,maxm1,coef(1,2),jcoef(1,2),rhs,u)!jgf46.00 comm.out
C     jgf46.00 Begin add.
      do j = 2,maxnz ! jgf50.03: st3 from maxnz-1, Bug fix, 06.02.2010 
      do i = 1,n
         rhs(i) = rhs(i) - coef(i,j)*u(jcoef(i,j))
      enddo
      enddo
C     jgf46.00 End add.
      return
      end  subroutine


      subroutine pmult (n,ndim,maxnz,jcoef,coef,b,c)
      implicit none
C      integer n,ndim,maxnz,maxm1 !jgf46.00 commented out
      integer i,j,n,ndim,maxnz    !jgf46.00 added
      integer jcoef(ndim,maxnz)
      real(sz) coef(ndim,maxnz),b(n),c(n)
c
c ... pmult computes c = a*b, a matrix-vector product.  matrix
c     a is assumed to be stored in the coef, jcoef ellpack
c     data structure and all entries in the column array jcoef
c     are assumed to be between 1 and n, inclusive.
c     a is assumed to have a unit diagonal.
c
c ... parameter list --
c
c          n        dimension of matrix
c          ndim     row dimension of coef and jcoef in calling routine
c          maxnz    maximum number of nonzeros per row
c          jcoef    integer array for coefficient columns
c          coef     array for coefficients
c          b        multiplying vector of length n
c          c        product vector of length n
c
c
      call scopy (n,b,1,c,1)
C      maxm1 = maxnz-1 !jgf46.00 commented out
C      call ypasx2 (ndim,n,maxm1,coef(1,2),jcoef(1,2),c,b) !jgf46.00 comm.out
C     jgf46.00 Begin add.
      do j = 2,maxnz !jgf50.03: st3 from maxnz-1, Bug fix, 06.02.2010
      do i = 1,n
         c(i) = c(i) - coef(i,j)*b(jcoef(i,j))
      enddo
      enddo
C     jgf46.00 End add.
      return
      end  subroutine

      subroutine prbndx (n,ndim,maxnz,jcoef,p,ip,nblack,level,nout,ier)
      implicit none
      integer i,n,ndim,maxnz,nblack,level,nout,ier,ibgn,next,last,
     *  j,k,nxttyp,jcol,l,nred
      integer jcoef(ndim,maxnz),p(n),ip(n)
      integer first,old,young,curtyp,type
c
c**************************************************************
c
c     prbndx computes the red-black permutation
c     vectors p ( and its inverse ip ) if possible.
c
c     the algorithm is to mark the first node as red (arbitrary).
c     all of its adjacent nodes are marked black and placed in
c     a stack.  the remainder of the code pulls the first node
c     off the top of the stack and tries to type its adjacent nodes.
c     the typing of the adjacent point is a five way case statement
c     which is well commented below (see do loop 100).
c
c     the array p is used both to keep track of the color of a node
c     (red node is positive, black is negative) but also the father
c     node that caused the color marking of that point.  since
c     complete information on the adjacency structure is hard to come
c     by this forms a link to enable the color change of a partial
c     tree when a recoverable color conflict occurs.
c
c     the array ip is used as a stack to point to the set of nodes
c     left to be typed that are known to be adjacent to the current
c     father node.
c
c     *** note ***  this routine is to be called after routine scal.
c
c*********************************************************************
c
c     input parameters --
c
c        n      number of nodes.  (integer, scalar)
c
c        ndim   row dimension of jcoef in calling routine.
c
c        maxnz  maximum number of nonzeros per row
c
c        jcoef  array of column indices.  it is assumed
c               that for every row where only one element is
c               stored that element corresponds to the diagonal
c               entry.  the diagonal must be the first entry stored.
c                 (integer, arrays)
c
c        level  switch for printing
c
c        nout   output tape number
c
c     output parameters --
c
c        nblack number of black nodes.  number of red nodes is
c               n - nblack.  (integer, scalar)
c
c        p, ip  permutation and inverse permutation vectors.
c               (integer, arrays each of length n)
c
c        ier    error flag. (integer, scalar)
c
c               ier = 0, normal return.  indexing performed
c                        successfully
c               ier = 201, red-black indexing not possible.
c
c********************************************************************
c
c
      ier = 0
      if (n.le.0) return
      do 10 i = 1,n
         p(i) = 0
         ip(i) = 0
   10 continue
c
c ... handle the first set of points until some adjacent points
c ... are found
c
      first = 1
c
   20 p(first) = first
      if (maxnz.gt.1) go to 40
c
c ... search for next entry that has not been marked
c
      if (first.eq.n) go to 120
      ibgn = first+1
      do 30 i = ibgn,n
         if (p(i).ne.0) go to 30
         first = i
         go to 20
   30 continue
      go to 120
c
c ... first set of adjacent points found
c
   40 next = 1
      last = 1
      ip(1) = first
c
c ... loop over labeled points indicated in the stack stored in
c ... the array ip
c
   50 k = ip(next)
      curtyp = p(k)
      nxttyp = -curtyp
      do 100 j = 2,maxnz
         jcol = jcoef(k,j)
         if (jcol.eq.k) go to 100
         type = p(jcol)
c
c==================================================================
c
c     the following is a five way case statement dealing with the
c     labeling of the adjacent node.
c
c ... case i.  if the adjacent node has already been labeled with
c              label equal to nxttyp, then skip to the next adjacent
c              node.
c
         if (type.eq.nxttyp) go to 100
c
c ... case ii.  if the adjacent node has not been labeled yet label
c               it with nxttyp and enter it in the stack
c
         if (type.ne.0) go to 60
         last = last+1
         ip(last) = jcol
         p(jcol) = nxttyp
         go to 100
c
c ... case iii.  if the adjacent node has already been labeled with
c                opposite color and the same father seed, then there
c                is an irrecoverable color conflict.
c
   60    if (type.eq.curtyp) go to 140
c
c ... case iv.  if the adjacent node has the right color and a different
c               father node, then change all nodes of the youngest fathe
c               node to point to the oldest father seed and retain the
c               same colors.
c
         if (type*nxttyp.lt.1) go to 80
         old = min0(iabs(type),iabs(nxttyp))
         young = max0(iabs(type),iabs(nxttyp))
         do 70 l = young,n
            if (iabs(p(l)).eq.young) p(l) = isign(old,p(l))
   70    continue
         curtyp = p(k)
         nxttyp = -curtyp
         go to 100
c
c ... case v.  if the adjacent node has the wrong color and a different
c              father node, then change all nodes of the youngest father
c              node to point to the oldest father node along with
c              changing their colors.  since until this time the
c              youngest father node tree has been independent no other
c              color conflicts will arise from this change.
c
   80    old = min0(iabs(type),iabs(nxttyp))
         young = max0(iabs(type),iabs(nxttyp))
         do 90 l = young,n
            if (iabs(p(l)).eq.young) p(l) = isign(old,-p(l))
   90    continue
         curtyp = p(k)
         nxttyp = -curtyp
c
c ... end of case statement
c
c==================================================================
c
  100 continue
c
c ... advance to next node in the stack
c
      next = next+1
      if (next.le.last) go to 50
c
c ... all nodes in the stack have been removed
c
c ... check for nodes not labeled.  if any are found
c ... start the labeling process again at the first
c ... node found that is not labeled.
c
      ibgn = first+1
      do 110 i = ibgn,n
         if (p(i).ne.0) go to 110
         first = i
         go to 20
  110 continue
c
c===================================================================
c
c ... all nodes are now typed either red or black
c
c ... generate permutation vectors
c
  120 call whenige (n,p,1,0,ip,nred)
      call whenilt (n,p,1,0,ip(nred+1),nblack)
      do 130 i = 1,n
         p(ip(i)) = i
  130 continue
c
c ... successful red-black ordering completed
c
      return
c
c ...... type conflict
c
  140 ier = 201
      if (level.ge.0) write (AMPI_LUN(nout),150)
  150 format (//1x,'*** f a t a l     e r r o r ************'//1x,
     *   '    in itpackv routine prbndx  '/1x,
     *   '    red-black indexing not possible')
      return
      end  subroutine

      subroutine sbelm (n,ndim,maxnz,jcoef,coef,rhs,work,tol)
      implicit none
      integer i,j,jcol,n,ndim,maxnz
      integer jcoef(ndim,maxnz)
      real(sz) coef(ndim,maxnz),rhs(n),work(n),tol
c
c ... sbelm is designed to remove rows of the matrix for which
c ... all off-diagonal elements are very small (less than tol).
c ... this is to take care of matrices arising from finite
c ... element discretizations of partial differential equations
c ... with dirichlet boundary conditions.  any such rows and
c ... corresponding columns are then eliminated (set to the
c ... identity after correcting the rhs).
c ... *** note ***  this routine is to be called after routine scal.
c
c ... parameter list --
c
c         n       dimension of matrix
c         ndim    row dimension of arrays jcoef and coef in the
c                    calling program
c         maxnz   maximum number of nonzero entries per row
c         jcoef   integer array of matrix representation
c         coef    array of sparse matrix representation
c         rhs     right hand side of matrix problem
c         work    work array of length n
c         tol     tolerance factor
c
c
      if (n.le.0.or.maxnz.lt.2) return
c
c ... find maximum off-diagonal elements in absolute value.
c
      do 10 i = 1,n
         work(i) = 0.0
   10 continue
      do 30 j = 2,maxnz
         do 20 i = 1,n
            work(i) = max(work(i),abs(coef(i,j)))
   20    continue
   30 continue
c
c ... eliminate desired rows and columns.
c
      do 60 j = 2,maxnz
         do 50 i = 1,n
            if (work(i).lt.tol) go to 40
            jcol = jcoef(i,j)
            if (work(jcol).ge.tol) go to 50
            rhs(i) = rhs(i)-coef(i,j)*rhs(jcol)
   40       coef(i,j) = 0.0
            jcoef(i,j) = i
   50    continue
   60 continue
      return
      end  subroutine

      subroutine scal (n,ndim,maxnz,jcoef,coef,rhs,u,work,ier)
      implicit none
      integer i,j,n,ndim,maxnz,nsgncg,ier
      integer jcoef(ndim,maxnz)
      real(sz) coef(ndim,maxnz),rhs(n),u(n),work(n),save
c
c *** begin -- itpackv common
c
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
c
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
c
      real(sz) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
c *** end   -- itpackv common
c
c ... scal scales original matrix to a unit diagonal matrix.  rhs
c ... and u vectors are scaled accordingly.  the data
c ... structure is adjusted to have diagonal entries in
c ... column 1.  zero entries in jcoef array are changed to
c ... positive integers between 1 and n.
c
c ... parameter list --
c
c         n       dimension of matrix
c         ndim    row dimension of arrays jcoef and coef in the
c                    calling program
c         maxnz   maximum number of nonzero entries per row
c         jcoef   integer array of matrix representation
c         coef    array of sparse matrix representation
c         rhs     right hand side of matrix problem
c         u       latest estimate of solution
c         work    work array of length n
c         ier     error flag -- on return, nonzero values mean
c                    401 -- zero diagonal element
c                    402 -- nonexistent diagonal element
c
c
c     description of variables in common block in main routine
c
c ... check for positive diagonal entries for each row.
c ... put diagonal entries in column 1.  replace zeros in
c ... row i of jcoef with i.
c
      ier = 0
      nsgncg = 0
      if (n.le.0) return
      do 110 i = 1,n
         if (jcoef(i,1).eq.i) go to 50
         if (maxnz.lt.2) go to 20
         do 10 j = 2,maxnz
            if (jcoef(i,j).eq.i) go to 40
   10    continue
c
c ... fatal error -- no diagonal entry for row i.
c
   20    ier = 402
         if (level.ge.0) write (AMPI_LUN(nout),30) i
   30    format (//1x,'*** f a t a l     e r r o r ************'//1x,
     *      '    in itpackv routine scal    '/1x,
     *      '    no diagonal entry in row',i10)
         return
c
c ... shift row i so that diagonal element is in column 1.
c
   40    save = coef(i,j)
         coef(i,j) = coef(i,1)
         jcoef(i,j) = jcoef(i,1)
         coef(i,1) = save
         jcoef(i,1) = i
c
c ... check sign of diagonal entry.  if negative, change signs of
c ... all row coefficients and corresponding rhs element.
c
   50    if (coef(i,1)) 60 , 90 , 110
   60    do 70 j = 1,maxnz
            coef(i,j) = -coef(i,j)
   70    continue
         rhs(i) = -rhs(i)
         nsgncg = nsgncg+1
         if (level.ge.5) write (AMPI_LUN(nout),80) i
   80    format (//1x,'*** n o t e ***'//1x,
     *      '    in itpackv routine scal'/1x,
     *      '    equation ',i10,' has been negated')
         go to 110
c
c ... fatal error -- zero diagonal element for row i.
c
   90    ier = 401
         if (level.ge.0) write (AMPI_LUN(nout),100) i
  100    format (//1x,'*** f a t a l     e r r o r ************'//1x,
     *      '    in itpackv routine scal    '/1x,
     *      '    diagonal entry in row ',i10,' is zero')
         return
  110 continue
c
c ... change zero elements of jcoef array.
c
      if (maxnz.lt.2) go to 140
      do 130 j = 2,maxnz
         do 120 i = 1,n
            if (jcoef(i,j).le.0) jcoef(i,j) = i
  120    continue
  130 continue
c
c ... scale rhs and u arrays.  store reciprocal square roots
c ... of diagonal entries in column 1 of coef.
c
  140 do 150 i = 1,n
         work(i) = sqrt(coef(i,1))
  150 continue
      do 160 i = 1,n
         u(i) = u(i)*work(i)
  160 continue
      do 170 i = 1,n
         work(i) = 1.0/work(i)
  170 continue
      call scopy (n,work,1,coef,1)
      do 180 i = 1,n
         rhs(i) = rhs(i)*work(i)
  180 continue
c
c ... scale matrix.
c
      if (maxnz.lt.2) return
      do 200 j = 2,maxnz
         do 190 i = 1,n
            coef(i,j) = coef(i,j)*work(i)*work(jcoef(i,j))
  190    continue
  200 continue
c
c ... adjust isym if the  0 .lt. nsgncg .lt. n
c
      if (nsgncg.gt.0.and.nsgncg.lt.n) isym = 1
c
      return
      end  subroutine


      subroutine unscal (n,ndim,maxnz,jcoef,coef,rhs,u,work)
      implicit none
      integer i,j,n,ndim,maxnz
      integer jcoef(ndim,maxnz)
      real(sz) coef(ndim,maxnz),rhs(n),u(n),work(n)
c
c ... unscal reverses the scaling done in routine scal.
c
c ... parameter list --
c
c         n       dimension of matrix
c         ndim    row dimension of arrays jcoef and coef in the
c                    calling program
c         maxnz   maximum number of nonzero entries per row
c         jcoef   integer array of matrix representation
c         coef    array of sparse matrix representation
c         rhs     right hand side of matrix problem
c         u       latest estimate of solution
c         work    work array of length n
c
c
c ... unscale u and rhs arrays.
c
      call scopy (n,coef,1,work,1)
      do 10 i = 1,n
         u(i) = u(i)*work(i)
   10 continue
      do 20 i = 1,n
         work(i) = 1.0/work(i)
   20 continue
      do 30 i = 1,n
         rhs(i) = rhs(i)*work(i)
   30 continue
c
c ... unscale matrix.
c
      if (maxnz.lt.2) go to 80
      do 50 j = 2,maxnz
         do 40 i = 1,n
            coef(i,j) = coef(i,j)*work(i)*work(jcoef(i,j))
   40    continue
   50 continue
c
c ... put original zeros back in icoef array and restore unscaled
c ... diagonal entries to column one.
c
      do 70 j = 2,maxnz
         do 60 i = 1,n
            if (jcoef(i,j).eq.i) jcoef(i,j) = 0
   60    continue
   70 continue
   80 do 90 i = 1,n
         coef(i,1) = work(i)**2
   90 continue
      return
      end  subroutine


      subroutine vout (n,v,iswt,nout)
      implicit none
c
c     vout effects printing of residual and solution
c     vectors - called from peror
c
c ... parameter list --
c
c          v      vector of length n
c          iswt   labelling information
c          nout   output device number
c
c ... specifications for arguments
c
      integer n, iswt, nout, k, kupper, i, j, jm1
      real(sz)  v(n)
c
c        if (n .le. 0) return
c
      kupper = min0(n,4)
      if (iswt.eq.1) write (AMPI_LUN(nout),10)
   10 format (//5x,'residual vector')
      if (iswt.eq.2) write (AMPI_LUN(nout),20)
   20 format (//5x,'solution vector')
      write (AMPI_LUN(nout),30) (i,i=1,kupper)
   30 format (10x,4i15)
      write (AMPI_LUN(nout),40)
   40 format (10x,65('-')/)
c
      do 60 j = 1,n,4
         kupper = min0(j+3,n)
         jm1 = j-1
         write (AMPI_LUN(nout),50) jm1,(v(k),k=j,kupper)
   50    format (4x,i5,'+  ',4e15.5)
   60 continue
c
      return
      end  subroutine

      subroutine zbrent (n,tri,eps,nsig,a,b,maxfn,ier)
      USE SIZES
      implicit none
      integer n,ier,nsig,maxfn,ib3,ib4,ic
      real(sz) a,b,c,d,e,eps,p,q,r,s,t,fa,fb,fc,tol,rm,rone,temp
      real(sz) tri(*)
c
c *** begin -- itpackv common
c
      integer in,is,isym,itmax,level,nout,numwav
      common /itcom1/ in,is,isym,itmax,level,nout,numwav
c
      logical adapt,betadt,caseii,halt,partad
      common /itcom2/ adapt,betadt,caseii,halt,partad
c
      real(sz) bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
      common /itcom3/ bdelnm,betab,cme,delnnm,delsnm,ff,gamma,omega,qa,
     *   qt,rho,rrr,sige,sme,specr,spr,srelpr,stptst,udnm,zeta
c
c *** end   -- itpackv common
c
c   modified imsl routine name   - zbrent
c
c-----------------------------------------------------------------------
c
c   computer            - cdc/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - zero of a function which changes sign in a
c                           given interval (brent algorithm)
c
c   usage               - call zbrent (f,eps,nsig,a,b,maxfn,ier)
c
c   arguments    tri    - a tridiagonal matrix of order n
c                eps    - first convergence criterion (input).  a root,
c                           b, is accepted if abs(f(b)) is less than or
c                           equal to eps.  eps may be set to zero.
c                nsig   - second convergence criterion (input).  a root,
c                           b, is accepted if the current approximation
c                           agrees with the true solution to nsig
c                           significant digits.
c                a,b    - on input, the user must supply two points, a
c                           and b, such that f(a) and f(b) are opposite
c                           in sign.
c                           on output, both a and b are altered.  b
c                           will contain the best approximation to the
c                           root of f. see remark 1.
c                maxfn  - on input, maxfn should contain an upper bound
c                           on the number of function evaluations
c                           required for convergence.  on output, maxfn
c                           will contain the actual number of function
c                           evaluations used.
c                ier    - error parameter. (output)
c                         terminal error
c                           ier = 501 indicates the algorithm failed to
c                             converge in maxfn evaluations.
c                           ier = 502 indicates f(a) and f(b) have the
c                             same sign.
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c   remarks  1.  let f(x) be the characteristic function of the matrix
c                tri evaluated at x. function determ evaluates f(x).
c                on exit from zbrent, when ier=0, a and b satisfy the
c                following,
c                f(a)*f(b) .le. 0,
c                abs(f(b)) .le. abs(f(a)), and
c                either abs(f(b)) .le. eps or
c                abs(a-b) .le. max(abs(b),0.1)*10.0**(-nsig).
c                the presence of 0.1 in this error criterion causes
c                leading zeroes to the right of the decimal point to be
c                counted as significant digits. scaling may be required
c                in order to accurately determine a zero of small
c                magnitude.
c            2.  zbrent is guaranteed to reach convergence within
c                k = (log((b-a)/d)+1.0)**2 function evaluations where
c                  d=min(over x in (a,b) of
c                    max(abs(x),0.1)*10.0**(-nsig)).
c                this is an upper bound on the number of evaluations.
c                rarely does the actual number of evaluations used by
c                zbrent exceed sqrt(k). d can be computed as follows,
c                  p = amin1(abs(a),abs(b))
c                  p = amax1 (0.1,p)
c                  if ((a-0.1)*(b-0.1).lt.0.0) p = 0.1
c                  d = p*10.0**(-nsig)
c
c   copyright           - 1977 by imsl, inc. all rights reserved.
c
c   warranty            - imsl warrants only that imsl testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c-----------------------------------------------------------------------
c
c
c     description of variables in common block in main routine
c
c ... local itpackv references --
c
c          determ
c
c                                  first executable statement
c
      ier = 0
      ib3 = 2*itmax+1
      ib4 = 3*itmax+1
      t = 10.0d0**(-nsig) !jgf46.00 added d0
      ic = 2
      fa = determ(itmax,n,tri,tri(ib3),tri(ib4),a)
      fb = determ(itmax,n,tri,tri(ib3),tri(ib4),b)
      s = b
c
c                                  test for same sign
c
C      if (fa*fb.gt.0.0) go to 110 jgf46.00 commented out
C     jgf46.00 Begin add.
c cvjp 2/10/06  replaced because produced NaNs
c     if (fa*fb.gt.0.0d0) go to 110
            
      if ( (fa > 0.0d0 .and. fb > 0.0d0) .or.
     &     (fa < 0.0d0 .and. fb < 0.0d0) ) then
        go to 110
      endif
C     jgf46.00 End add.

   10 c = a
      fc = fa
      d = b-c
      e = d
   20 if (abs(fc).ge.abs(fb)) go to 30
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
   30 continue

#ifdef REAL8
      tol = t*max(abs(b),0.1d0)
#else
      tol = t*max(abs(b),0.1e0)
#endif

      rm = (c-b)*0.5d0
c
c                                  test for first convergence criteria
c
      if (abs(fb).le.eps) go to 80
c
c                                  test for second convergence criteria
c
      if (abs(c-b).le.tol) go to 80
c
c                                  check evaluation counter
c
      if (ic.ge.maxfn) go to 90
c
c                                  is bisection forced
c
      if (abs(e).lt.tol) go to 60
      if (abs(fa).le.abs(fb)) go to 60
      s = fb/fa
      if (a.ne.c) go to 40
c
c                                  linear interpolation
c
      p = (c-b)*s
      q = 1.0d0-s
      go to 50
c
c                                  inverse quadratic interpolation
c
   40 q = fa/fc
      r = fb/fc
      rone = r-1.0d0
      p = s*((c-b)*q*(q-r)-(b-a)*rone)
      q = (q-1.0d0)*rone*(s-1.0d0)
   50 if (p.gt.0.0d0) q = -q
      if (p.lt.0.0d0) p = -p
      s = e
      e = d
c
c                                  if abs(p/q).ge.75*abs(c-b) then
c                                     force bisection
c
      if (p+p.ge.3.0d0*rm*q) go to 60
c
c                                  if abs(p/q).ge..5*abs(s) then force
c                                     bisection. s = the value of p/q
c                                     on the step before the last one
c
      if (p+p.ge.abs(s*q)) go to 60
      d = p/q
      go to 70
c
c                                  bisection
c
   60 e = rm
      d = e
c
c                                  increment b
c
   70 a = b
      fa = fb
      temp = d
cjjw/vjpm002 - modified/added the following 5 lines
#ifdef REAL8
      if (abs(temp).le.0.5d0*tol) temp = sign(0.5d0*tol,rm)
#else
      if (abs(temp).le.0.5e0*tol) temp = sign(0.5e0*tol,rm)
#endif
      b = b+temp
      s = b
      fb = determ(itmax,n,tri,tri(ib3),tri(ib4),s)
      ic = ic+1
C      if (fb*fc.le.0.0) go to 20 !jgf46.00 commented out
C     jgf46.00 Begin add. 
cvjp 2/10/06 replace because produced NaNs
c     if (fb*fc.le.0.0d0) go to 20

      if (fb == 0.0d0 .or. fc == 0.0d0) goto 20

      if ( (fb > 0.0d0 .and. fc < 0.0d0) .or.
     &     (fb < 0.0d0 .and. fc > 0.0d0) ) then
        go to 20
      endif
C     jgf46.00 End add.

      go to 10
c
c                                  convergence of b
c
   80 a = c
      maxfn = ic
      go to 130
c
c                                  maxfn evaluations
c
   90 ier = 501
      a = c
      maxfn = ic
      if (level.ge.1) write (AMPI_LUN(nout),100) maxfn
  100 format (/1x,'*** w a r n i n g ************'//1x,
     *   '    in itpackv routine zbrent'/1x,
     *   '    algorithm failed to converge'/1x,
     *   '    in',i6,' iterations ')
      go to 130
c
c                                  terminal error - f(a) and f(b) have
c                                  the same sign
c
  110 ier = 502
      maxfn = ic
      if (level.ge.1) write (AMPI_LUN(nout),120)
  120 format (/1x,'*** w a r n i n g ************'//1x,
     *   '    in itpackv routine zbrent  '/1x,
     *   '    f(a) and f(b) have same sign   ')
  130 continue
      return
      end  subroutine

C     jgf46.00 deleted the following subroutine
C      subroutine ypasx2 (ndim,n,m,a,ja,y,x)

C     jgf46.00 deleted the following subroutine
C      subroutine ymasx2 (ndim,n,m,a,ja,y,x)


      subroutine whenige (n,p,inc,itarg,ip,npt)
      implicit none
      integer n,inc,itarg,npt,i
      integer p(n), ip(n)
c
      npt = 0
      do 10 i = 1,n
         if (p(i) .lt. itarg) go to 10
         npt = npt + 1
         ip(npt) = i
 10   continue
      return
      end subroutine


      subroutine whenilt (n,p,inc,itarg,ip,npt)
      implicit none
      integer n,inc,itarg,npt,i
      integer p(n), ip(n)
c
      npt = 0
      do 10 i = 1,n
         if (p(i) .ge. itarg) go to 10
         npt = npt + 1
         ip(npt) = i
 10   continue
      return
      end subroutine


C
C*****************************************************************
C
C  Modified BLAS for ITPACK2D
C
C*****************************************************************
C

      real(sz) function sdot(n,sx,incx,sy,incy)
      implicit none
      integer i,n,incx,incy
      real(sz) sx(*),sy(*)
      REAL*8 ddot
c
      ddot = 0.0D0
      sdot = ddot
c
      if (n.gt.0) then  
        do i = 1,n
           ddot = ddot + sx(i)*sy(i)
        enddo 
        sdot = ddot
      endif
c
      return
      end function

      integer function ismin (n,sx,incx)
      implicit none
      integer n,incx,ns,i,ii
      real(sz) sx(*),smin,xval
c
c     find smallest index of minimum value of single precision sx.
c
      ismin = 0
      if (n.le.0) return
      ismin = 1
      if (n.le.1) return
      if (incx.eq.1) go to 30
c
c        code for increments not equal to 1.
c
      smin = sx(1)
      ns = n*incx
      ii = 1
      do 20 i = 1,ns,incx
         xval = sx(i)
         if (xval.ge.smin) go to 10
         ismin = ii
         smin = xval
   10    ii = ii+1
   20 continue
      return
c
c        code for increments equal to 1.
c
   30 smin = sx(1)
      do 40 i = 2,n
         xval = sx(i)
         if (xval.ge.smin) go to 40
         ismin = i
         smin = xval
   40 continue
      return
      end function



      subroutine  scopy(n,sx,incx,sy,incy)
      implicit none
      integer n,incx,incy,ix,iy,i,m,mp1,ns
      real(sz) sx(*),sy(*)
c
c     copy single precision sx to single precision sy.
c
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
c
c        code for unequal or nonpositive increments.
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop so remaining vector length is a multiple of 7.
c
   20 m = n - (n/7)*7
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
c
c        code for equal, positive, nonunit increments.
c
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          sy(i) = sx(i)
   70     continue
      return
      end subroutine


      END MODULE ITPACKV
