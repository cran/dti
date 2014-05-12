      subroutine nlrdtirg(s,nb,n,mask,b,sdcoef,th0,D,niter,eps,
     1                    res,rss,varinv)
C
C  this is based on a regularized tensor similar to Koay et.al. (2006)
C
      implicit logical (a-z)
      integer nb,n,niter
      real*8 s(nb,n)
      logical mask(n)
      real*8 D(6,n),b(6,nb),res(nb,n),th0(n),eps,rss(n),sdcoef(4),
     1       varinv(nb)
      integer i,j
      DO i=1,n
         if(mask(i)) THEN
            call islvdti(s(1,i),nb,b,sdcoef,varinv,th0(i),D(1,i),
     1                       res(1,i),niter,eps,rss(i))
         ELSE
            DO j=1,6
               D(j,i)=0.d0
            END DO
            rss(i)=0.d0
         END IF
      call rchkusr()         
      END DO
      RETURN
      END
      subroutine nlrdtirp(s,nb,n,b,sdcoef,th0,niter,eps,res,npar,
     1                    varinv)
C
C  this is based on a regularized tensor similar to Koay et.al. (2006)
C  modified to use parallelization  with function pmatrix
C  npar = 8+nb
C  res(1,) used for th0
C  res(2:7,) used for D
C  res(8,) used for rss
C  res(9:npar,) used for residuals
      implicit logical (a-z)
      integer nb,n,npar,niter
      real*8 s(nb,n)
      real*8 b(6,nb),res(npar,n),th0(n),eps,sdcoef(4),varinv(nb)
      integer i
      DO i=1,n
            res(1,i) = th0(i)
            call islvdti(s(1,i),nb,b,sdcoef,varinv,res(1,i),res(2,i),
     1                       res(9,i),niter,eps,res(8,i))
      call rchkusr()         
      END DO
      RETURN
      END
C
C
C
      subroutine islvdtir(s,nb,b,varinv,th0,D,F,niter,eps,rss)
C
C  Implements the regularized Gauss-Newton Algortithm (10.2.8)
C  from Schwetlick (1979)
C
      implicit logical (a-z)
      integer nb,niter
      real*8 s(nb)
      real*8 D(6),rho(6),b(6,nb),th0,F(nb),eps
      integer i,j,k,info,iter,icount
      real*8 z,gamma,alpha,delta,varinv(nb),
     1       dg(7),pk(7),ak(7,7),ck(7,7),rss,nrss,crss,maxabsdg,
     2       oldrss,relrss,rhon(6),Dn(6),res,X(7),th0n,XX(6),xzvarinv
      logical negdefin
      delta=0.25D0
      gamma=1.d0
      alpha=0.7d0
      oldrss=1.d50
      rss=0.d0
C      call regD(D,D0)
C   for test purposes only
      call D2rho(D,rho)
C includes regularization of D
      DO i=1,nb
         z=b(1,i)*D(1)
         DO j=2,6
            z=z+b(j,i)*D(j)
         END DO
         z=exp(-z)
         res=s(i)-th0*z
         rss=rss+res*res*varinv(i)
         F(i)=res
      END DO
      th0n = th0
      nrss = rss
      DO iter=1,niter
         DO j=1,7
            dg(j)=0.d0
            DO k=j,7
               ak(j,k)=0.d0
            END DO
         END DO
         DO i=1,nb
            z=b(1,i)*D(1)
            DO j=2,6
               z=z+b(j,i)*D(j)
            END DO
            z=exp(-z)
            X(7)= -z
            z=z*th0
            DO j=1,6
               XX(j)=b(j,i)*z
            END DO
            X(1)=2.d0*rho(1)*XX(1)+rho(2)*XX(2)+rho(3)*XX(3)
            X(2)=rho(1)*XX(2)+2.d0*rho(2)*XX(4)+rho(3)*XX(5)
            X(3)=rho(1)*XX(3)+rho(2)*XX(5)+2.d0*rho(3)*XX(6)
            X(4)=2.d0*rho(4)*XX(4)+rho(5)*XX(5)
            X(5)=rho(4)*XX(5)+2.d0*rho(5)*XX(6)
            X(6)=2*rho(6)*XX(6)
            DO j=1,7
               xzvarinv=X(j)*varinv(i)
               dg(j)=dg(j)+xzvarinv*F(i)
               DO k=j,7
                  ak(j,k)=ak(j,k)+xzvarinv*X(k)
               END DO
            END DO 
         END DO
         maxabsdg=0.d0
         DO j=1,7
            maxabsdg=max(maxabsdg,abs(dg(j)))
         END DO
         relrss = (oldrss-rss)/rss
         IF(maxabsdg.lt.eps.or.relrss.lt.1d-6) THEN
C  prepare things for return if gradient is close to 0
            call regularD(D,negdefin)
            RETURN
         END IF
         IF(iter.gt.1.and.abs(rho(1)*rho(4)*rho(6)).lt.1.d-10) THEN
C  prepare things for return if gradient is close to 0
            call regularD(D,negdefin)
            RETURN
         END IF
         gamma=min(gamma/alpha,1.d0)
C  End of step 3
         notacc=.TRUE.
         icount = 10
         DO WHILE (notacc.and.icount.gt.0) 
            icount = icount-1
            call rchkusr()
            IF(gamma.lt.1.d0) THEN
               DO j=1,7
                  DO k=j,7
                     ck(j,k)=gamma*ak(j,k)
                  END DO
                  ck(j,j)=ck(j,j)+1.d0-gamma
               END DO
            ELSE
C   we may still need ak and dg so copy them to pk and ck
               DO j=1,7
                  DO k=j,7
                     ck(j,k)=ak(j,k)
                  END DO
               END DO
            END IF
            DO j=1,7
               pk(j)=dg(j)
            END DO
C   Now solve  ak%*%dtheta= dg
            call dposv("U",7,1,ck,7,pk,7,info)
C  Step 4 we have pk 
            IF(info.ne.0) THEN
               gamma=alpha*gamma
C  thats step 6
            ELSE
C  comute things needed for decision in step 5 
C  if successful F, nrss, and theta will be reused in the  
C  next iteration
               
               DO j=1,6
                  rhon(j)=rho(j)-gamma*pk(j)
               END DO
               th0n=th0-gamma*pk(7)
               nrss=0.d0
               call rho2D0(rhon,Dn)
               DO i=1,nb
                  z=b(1,i)*Dn(1)
                  DO j=2,6
                     z=z+b(j,i)*Dn(j)
                  END DO
                  res=(s(i)-th0n*exp(-z))
                  nrss=nrss+res*res*varinv(i)
                  F(i)=res
               END DO
               crss=dg(1)*pk(1)
               DO j=2,7
                  crss=crss+dg(j)*pk(j)
               END DO
               crss=rss-delta*gamma*crss
               IF(nrss.le.crss) THEN
                  notacc=.FALSE.
C  accept new estimate, prepare for next iteration
               ELSE
                  gamma=alpha*gamma
C  decrease gamma and try new regularization
               END IF
            END IF
         END DO
         th0=th0n
         DO j=1,6
            rho(j)=rhon(j)
            D(j)=Dn(j)
         END DO
         oldrss=rss
         rss=nrss
         call rchkusr()
      END DO
      call regularD(D,negdefin)
      RETURN
      END
      subroutine dslvdtir(s,nb,b,varinv,th0,D,F,niter,eps,rss)
C
C  Implements the regularized Gauss-Newton Algortithm (10.2.8)
C  from Schwetlick (1979)
C
      implicit logical (a-z)
      integer nb,niter
      real*8 s(nb),D(6),rho(6),b(6,nb),varinv(nb),th0,F(nb),eps
      integer i,j,k,info,iter,icount
      real*8 z,gamma,alpha,delta,xzvarinv,
     1       dg(7),pk(7),ak(7,7),ck(7,7),rss,nrss,crss,maxabsdg,
     2       oldrss,relrss,rhon(6),Dn(6),res,X(7),th0n,XX(6)
      logical negdefin
      delta=0.25D0
      gamma=1.d0
      alpha=0.7d0
      oldrss=1.d50
      rss=0.d0
      call D2rho(D,rho)
      DO i=1,nb
         z=b(1,i)*D(1)
         DO j=2,6
            z=z+b(j,i)*D(j)
         END DO
         z=exp(-z)
         res=s(i)-th0*z
         rss=rss+res*res*varinv(i)
         F(i)=res
      END DO
      th0n = th0
      nrss = rss
      DO iter=1,niter
         DO j=1,7
            dg(j)=0.d0
            DO k=j,7
               ak(j,k)=0.d0
            END DO
         END DO            
         DO i=1,nb
            z=b(1,i)*D(1)
            DO j=2,6
               z=z+b(j,i)*D(j)
            END DO
            z=exp(-z)
            X(7)= -z
            z=z*th0
            DO j=1,6
               XX(j)=b(j,i)*z
            END DO
            X(1)=2.d0*rho(1)*XX(1)+rho(2)*XX(2)+rho(3)*XX(3)
            X(2)=rho(1)*XX(2)+2.d0*rho(2)*XX(4)+rho(3)*XX(5)
            X(3)=rho(1)*XX(3)+rho(2)*XX(5)+2.d0*rho(3)*XX(6)
            X(4)=2.d0*rho(4)*XX(4)+rho(5)*XX(5)
            X(5)=rho(4)*XX(5)+2.d0*rho(5)*XX(6)
            X(6)=2*rho(6)*XX(6)
            DO j=1,7
               xzvarinv=X(j)*varinv(i)
               dg(j)=dg(j)+xzvarinv*F(i)
               DO k=j,7
                  ak(j,k)=ak(j,k)+xzvarinv*X(k)
               END DO
            END DO 
         END DO
         maxabsdg=0.d0
         DO j=1,7
            maxabsdg=max(maxabsdg,abs(dg(j)))
         END DO
         relrss = (oldrss-rss)/rss
         IF(maxabsdg.lt.eps.or.relrss.lt.1d-6) THEN
C  prepare things for return if gradient is close to 0
            call regularD(D,negdefin)
            RETURN
         END IF
         IF(iter.gt.1.and.abs(rho(1)*rho(4)*rho(6)).lt.1.d-10) THEN
C  prepare things for return if gradient is close to 0
            call regularD(D,negdefin)
            RETURN
         END IF
         gamma=min(gamma/alpha,1.d0)
C  End of step 3
         notacc=.TRUE.
         icount = 10
         DO WHILE (notacc.and.icount.gt.0)
            icount=icount-1 
            IF(gamma.lt.1.d0) THEN
               DO j=1,7
                  DO k=j,7
                     ck(j,k)=gamma*ak(j,k)
                  END DO
                  ck(j,j)=ck(j,j)+1.d0-gamma
               END DO
            ELSE
C   we may still need ak and dg so copy them to pk and ck
               DO j=1,7
                  DO k=j,7
                     ck(j,k)=ak(j,k)
                  END DO
               END DO
            END IF
            DO j=1,7
               pk(j)=dg(j)
            END DO
C   Now solve  ak%*%dtheta= dg
            call dposv("U",7,1,ck,7,pk,7,info)
C  Step 4 we have pk 
            IF(info.ne.0) THEN
               gamma=alpha*gamma
C  thats step 6
            ELSE
C  comute things needed for decision in step 5 
C  if successful F, nrss, and theta will be reused in the  
C  next iteration
               
               DO j=1,6
                  rhon(j)=rho(j)-gamma*pk(j)
               END DO
               th0n=th0-gamma*pk(7)
               nrss=0.d0
               call rho2D0(rhon,Dn)
               DO i=1,nb
                  z=b(1,i)*Dn(1)
                  DO j=2,6
                     z=z+b(j,i)*Dn(j)
                  END DO
                  res=s(i)-th0n*exp(-z)
                  nrss=nrss+res*res*varinv(i)
                  F(i)=res
               END DO
               crss=dg(1)*pk(1)
               DO j=2,7
                  crss=crss+dg(j)*pk(j)
               END DO
               crss=rss-delta*gamma*crss
               IF(nrss.le.crss) THEN
                  notacc=.FALSE.
C  accept new estimate, prepare for next iteration
               ELSE
                  gamma=alpha*gamma
C  decrease gamma and try new regularization
               END IF
            END IF
         END DO
         th0=th0n
         DO j=1,6
            rho(j)=rhon(j)
            D(j)=Dn(j)
         END DO
         oldrss=rss
         rss=nrss
         call rchkusr()
      END DO
      call regularD(D,negdefin)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     get D from rho without regularization
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rho2D0(rho,D)
      implicit logical(a-z)
      real*8 D(6),rho(6),r1,r2,r3,r4,r5,r6
      r1=rho(1)
      r2=rho(2)
      r3=rho(3)
      r4=rho(4)
      r5=rho(5)
      r6=rho(6)
      D(1)=r1*r1
      D(2)=r1*r2
      D(3)=r1*r3
      D(4)=r2*r2+r4*r4
      D(5)=r2*r3+r4*r5
      D(6)=r3*r3+r5*r5+r6*r6
      RETURN
      END
C
C
C   compute Least Squares Criterion based on positive definit tensor parametrization D = R^T R
C
C     
!       subroutine opttensR(par,si,s0,grad,ngrad,sdcoef,erg)
!       implicit logical (a-z)
!       integer ngrad
!       real*8 par(6),si(ngrad),s0,grad(3,ngrad),sdcoef(4),erg
!       integer i
!       real*8 low,up,g1,g2,g3,s,D(6),z,w
!       low=sdcoef(1)+sdcoef(3)*sdcoef(2)
!       up=sdcoef(1)+sdcoef(4)*sdcoef(2)
!       call rho2D0(par,D)
!       s=0.d0
!       DO i = 1,ngrad
!          g1 = grad(1,i)
!          g2 = grad(2,i)
!          g3 = grad(3,i)
!          z =  g1*g1*D(1)+g2*g2*D(4)+g3*g3*D(6)+
!      1        2.d0*(g1*g2*D(2)+g1*g3*D(3)+g2*g3*D(5))
!          w = max(low,(min(up,sdcoef(1)+si(i)*sdcoef(2))))
!          z = (si(i)-s0*exp(-z))/w
!          s = s+z*z
!       END DO
!       erg = s
!       RETURN
!       END      
!       subroutine tensRres(par,si,s0,grad,ngrad,res)
!       implicit logical (a-z)
!       integer ngrad
!       real*8 par(6),si(ngrad),s0,grad(3,ngrad)
!       integer i
!       real*8 g1,g2,g3,D(6),res(ngrad),z
!       call rho2D0(par,D)
!       DO i = 1,ngrad
!          g1 = grad(1,i)
!          g2 = grad(2,i)
!          g3 = grad(3,i)
!          z =  g1*g1*D(1)+g2*g2*D(4)+g3*g3*D(6)+
!      1        2.d0*(g1*g2*D(2)+g1*g3*D(3)+g2*g3*D(5))
!          res(i) = (si(i)-s0*exp(-z))
!       END DO
!       RETURN
!       END
