C
C  compute M(c) in Jiang Cheng (2012)
C
      subroutine Mofcall(c,kern,nk,ng,n,e,l,fvmofc,res)
      implicit logical (a-z)
      integer nk,ng,n
      real*8 c(nk,n),kern(nk,nk,ng),e(ng,n),l(nk),fvmofc(n),res(ng,n)
      integer i
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(c,kern,nk,ng,n,e,l,fvmofc,res)
C$OMP& PRIVATE(i)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         call Mofcres(c(1,i),kern,nk,ng,e(1,i),l,fvmofc(i),res(1,i))
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(res,fvmofc)
      RETURN
      END
C
C  compute M(c) in Jiang Cheng (2012)
C
      subroutine Mofc(c,kern,nk,ng,e,l,fvmofc)
      implicit logical (a-z)
      integer nk,ng
      real*8 c(nk),kern(nk,nk,ng),e(ng),l(nk),fvmofc
C diagonal matrix in l
      integer i,j,k
      real*8 z,ci,ci2,zz
C  first the penalty term
      z=0.d0
      DO i=1,nk
         ci=c(i)
         z=z+ci*ci*l(i)
      END DO
      DO k=1,ng
C now c^T K(g(k)) C  for gradient k
         zz=c(1)*c(1)*kern(1,1,k)
         DO i=2,nk
            ci=c(i)
            zz=zz+ci*ci*kern(i,i,k)
            ci2=2.d0*ci
            Do j=1,i-1
               zz=zz+ci2*c(j)*kern(j,i,k)
            END DO
         END DO
C  get residual in zz and penalized sum of squares in z
         zz=zz-e(k)
         z=z+zz*zz
      END DO
      fvmofc=z/2.d0
      RETURN
      END
      subroutine Mofcres(c,kern,nk,ng,e,l,fvmofc,res)
      implicit logical (a-z)
      integer nk,ng
      real*8 c(nk),kern(nk,nk,ng),e(ng),l(nk),fvmofc,res(ng)
C diagonal matrix in l
      integer i,j,k
      real*8 z,ci,ci2,zz
C  first the penalty term
      z=0.d0
      DO i=1,nk
         ci=c(i)
         z=z+ci*ci*l(i)
      END DO
      DO k=1,ng
         zz=c(1)*c(1)*kern(1,1,k)
         DO i=2,nk
            ci=c(i)
            zz=zz+ci*ci*kern(i,i,k)
            ci2=2.d0*ci
            Do j=1,i-1
               zz=zz+ci2*c(j)*kern(j,i,k)
            END DO
         END DO
         zz=zz-e(k)
         res(k)=zz
         z=z+zz*zz
      END DO
      fvmofc=z/2.d0
      RETURN
      END
C
C  compute dM(c)/dc in Jiang Cheng (2012)
C
      subroutine Dmofcdc(c,kern,nk,ng,e,l,w,mofcdc)
      implicit logical (a-z)
      integer nk,ng
      real*8 c(nk),kern(nk,nk,ng),e(ng),l(nk),w(nk),mofcdc(nk)
C diagonal matrix in l
      integer k
      real*8 z,ddot
      external ddot
C  first the penalty term  Lambda%*%c
      DO k=1,nk
         mofcdc(k) = c(k)* l(k)
      END DO      
C  now the sum
      DO k=1,ng
C  compute kern(,,k)%*%c  in w
         call DGEMV("N",nk,nk,1.d0,kern(1,1,k),nk,c,1,0.d0,w,1)
C  compute   2*(c^T%*%kern(,,k)%*%c-e(k)) in z
         z= 2.d0*(DDOT(nk,c,1,w,1)-e(k))
C  add 2*(c^T%*%kern(,,k)%*%c-e(k))*kern(,,k)%*%c to mofcdc
         call DAXPY(nk,z,w,1,mofcdc,1)
      END DO
      RETURN
      END
C
C  compute nabla M(c) in Jiang Cheng (2012)
C
      subroutine nablmofc(c,kern,nk,ng,e,l,w,nablam)
      implicit logical (a-z)
      integer nk,ng
      real*8 c(nk),kern(nk,nk,ng),e(ng),l(nk),w(nk),nablam(nk)
      real*8 z,ddot
      external ddot
C  first get  dM(c)/dc in   nablam
      call Dmofcdc(c,kern,nk,ng,e,l,w,nablam)
C  -c^T%*%dM(c)/dc in z
      z = -DDOT(nk,c,1,nablam,1)
      call DAXPY(nk,z,c,1,nablam,1)
      RETURN
      END
C
C  compute Exp_c(dt*v)  in Jiang Cheng (2012)
C
      subroutine expcv(c,v,nk,dt,r) 
      implicit logical (a-z)
      integer nk
      real*8 v(nk),c(nk),r(nk),dt
      integer i
      real*8 sv,cv,vn,dnrm2
      external dnrm2
      vn = dnrm2(nk,v,1)
      cv = dcos(dt*vn)
      sv = dsin(dt*vn)/vn
      DO i=1,nk
         r(i) = c(i)*cv-v(i)*sv
      END DO
      RETURN
      END
C
C  compute c^(k+!) from c^k in Jiang Cheng (2012)
C
      subroutine getnewck(c,kern,nk,ng,m0,e,l,w,nablam,ck1,mdt0,ierr)
      implicit logical (a-z)
      integer nk,ng
      real*8 c(nk),kern(nk,nk,ng),e(ng),l(nk),w(nk),nablam(nk),
     1       ck1(nk)
      real*8 dt0,dnrm2,m0,mdt0,normck
      integer i,ierr
      external dnrm2
      ierr=0
C first get nabla M(c) in nablam
      call nablmofc(c,kern,nk,ng,e,l,w,nablam)
C get norm of it 
C      dt0 = dnrm2(nk,nablam,1)
      dt0=.25d0
C perform line search 
C first M(c) in m0
C now M(dt*v/||v||) in mdt
      call expcv(c,nablam,nk,dt0,ck1) 
      call Mofc(ck1,kern,nk,ng,e,l,mdt0)
      winit = .TRUE.
C try to find valid start
      DO while (winit) 
         IF(mdt0.ge.m0) THEN
            dt0 = dt0*.6
            call expcv(c,nablam,nk,dt0,ck1) 
            call Mofc(ck1,kern,nk,ng,e,l,mdt0)
         ELSE
            winit = .FALSE.
         END IF
         if(dt0.lt.1d-6) THEN
C no descent found, keep c(nk)
            winit = .FALSE.
            ierr=1
         END IF
      END DO
C haben jetzt mdt0 <= m0 
      normck = dnrm2(nk,ck1,1)   
      if(abs(normck-1.d0).gt.1d-6) THEN
C          call dblepr("normck",6,normck,1)
          Do i = 1,nk
             ck1(i)=ck1(i)/normck
          End do
      END IF
      RETURN
      END
C
C   estimate coefficients
C      
      subroutine sqrteap(ei,kern,nk,ng,n,l,w,nablam,ck,ck1,cres)
      implicit logical (a-z)
      integer nk,ng,n
      real*8 ck(nk,1),kern(nk,nk,ng),ei(ng,n),l(nk),w(nk,1),
     1       nablam(nk,1),ck1(nk,1),cres(nk,n)
      logical ndone
      integer i,k,ierr,thrnr
      real*8 mck,mck1
      integer omp_get_thread_num
      external omp_get_thread_num
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ei,kern,nk,ng,n,l,w,nablam,ck,ck1,cres)
C$OMP& PRIVATE(i,k,ierr,mck,mck1,ndone,thrnr)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         thrnr = omp_get_thread_num()+1
         ck(1,thrnr)=1.d0
         DO k=2,nk
            ck(k,thrnr)=0.d0
         END DO
         call Mofc(ck(1,thrnr),kern,nk,ng,ei(1,i),l,mck)
         ndone=.TRUE.
         DO WHILE (ndone)
            call getnewck(ck(1,thrnr),kern,nk,ng,mck,ei(1,i),l,
     1              w(1,thrnr),nablam(1,thrnr),ck1(1,thrnr),mck1,ierr)
C            call Mofc(ck1,kern,nk,ng,ei(1,i),l,mck1)
            if(ierr.eq.0) THEN
               call DCOPY(nk,ck1(1,thrnr),1,ck(1,thrnr),1)
               ndone = (mck-mck1)/mck.ge.1d-7
               mck = mck1
            ELSE
               ndone = .FALSE.
C  kein Abswtiegspunkt gefunden
            END IF
         END DO
         call DCOPY(nk,ck(1,thrnr),1,cres(1,i),1)
C         call rchkusr()
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(cres)
      RETURN
      END
      
      
