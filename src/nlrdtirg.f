      subroutine nlrdtirg(s,nb,n1,n2,n3,mask,b,th0,D,niter,eps,
     1                    res,rss)
C
C  this is based on a regularized tensor similar to Koay et.al. (2006)
C
      implicit logical (a-z)
      integer nb,n1,n2,n3,s(nb,n1,n2,n3),niter
      logical mask(n1,n2,n3)
      real*8 D(6,n1,n2,n3),b(6,nb),res(nb,n1,n2,n3),
     1    th0(n1,n2,n3),eps,rss(n1,n2,n3)
      integer i1,i2,i3,j
      call intpr("niter",5,niter,1)
      DO i3=1,n3
         call intpr("Nonlinear regression for slice No:",34,i3,1)
         DO i2=1,n2
            DO i1=1,n1
C               call intpr("i1",2,i1,1)
               if(mask(i1,i2,i3)) THEN
                  call islvdti(s(1,i1,i2,i3),nb,b,th0(i1,i2,i3),
     1                       D(1,i1,i2,i3),
     2                       res(1,i1,i2,i3),niter,eps,rss(i1,i2,i3))
               ELSE
                  DO j=1,6
                     D(j,i1,i2,i3)=0.d0
                  END DO
                  rss(i1,i2,i3)=0.d0
               END IF
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine islvdtir(s,nb,b,th0,D,F,niter,eps,rss)
C
C  Implements the regularized Gauss-Newton Algortithm (10.2.8)
C  from Schwetlick (1979)
C
      implicit logical (a-z)
      integer nb,s(nb),niter
      real*8 D(6),rho(6),b(6,nb),th0,F(nb),eps
      integer i,j,k,info,iter,indvar,icount
      real*8 z,gamma,alpha,delta,
     1       dg(7),pk(7),ak(7,7),ck(7,7),rss,nrss,crss,maxabsdg,
     2       oldrss,relrss,rhon(6),Dn(6),res,X(7),th0n,XX(6)
      external indvar
C      call intpr("slvdtirg",8,1,1)
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
         z=0.d0
         DO j=1,6
            z=z+b(j,i)*D(j)
         END DO
         z=exp(-z)
         res=s(i)-th0*z
         rss=rss+res*res
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
            z=0.d0
            DO j=1,6
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
               dg(j)=dg(j)+X(j)*F(i)
               DO k=j,7
                  ak(j,k)=ak(j,k)+X(j)*X(k)
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
C      call dblepr("D0",2,D0,6)
C      call dblepr("Da",2,D,6)
C            call testreg(D,1)
C            call dblepr("rhon1",5,rhon,6)
C            call intpr("iter",4,iter,1)
            RETURN
         END IF
         IF(iter.gt.1.and.abs(rhon(1)*rhon(4)*rhon(6)).lt.1.d-10) THEN
C  prepare things for return if gradient is close to 0
            call regularD(D,negdefin)
C      call dblepr("D0",2,D0,6)
C      call dblepr("Db",2,D,6)
C            call testreg(D,1)
C            call dblepr("rhon3",5,rhon,6)
C            call intpr("iter",4,iter,1)
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
               call rho2D(rhon,Dn)
               DO i=1,nb
                  z=0.d0
                  DO j=1,6
                     z=z+b(j,i)*Dn(j)
                  END DO
                  res=s(i)-th0n*exp(-z)
                  nrss=nrss+res*res
                  F(i)=res
               END DO
               crss=0.d0
               DO j=1,7
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
C      call dblepr("D0",2,D0,6)
C      call dblepr("Dc",2,D,6)
C      call testreg(D,2)
C      call dblepr("rhon2",5,rhon,6)
C      call intpr("niter",5,niter,1)
      RETURN
      END
      subroutine dslvdtir(s,nb,b,th0,D,F,niter,eps,rss)
C
C  Implements the regularized Gauss-Newton Algortithm (10.2.8)
C  from Schwetlick (1979)
C
      implicit logical (a-z)
      integer nb,niter
      real*8 s(nb),D(6),rho(6),b(6,nb),th0,F(nb),eps
      integer i,j,k,info,iter,indvar,icount
      real*8 z,gamma,alpha,delta,
     1       dg(7),pk(7),ak(7,7),ck(7,7),rss,nrss,crss,maxabsdg,
     2       oldrss,relrss,rhon(6),Dn(6),res,X(7),th0n,XX(6)
      external indvar
      delta=0.25D0
      gamma=1.d0
      alpha=0.7d0
      oldrss=1.d50
      rss=0.d0
C      call regularD(D)
      call D2rho(D,rho)
      DO i=1,nb
         z=0.d0
         DO j=1,6
            z=z+b(j,i)*D(j)
         END DO
         z=exp(-z)
         res=s(i)-th0*z
         rss=rss+res*res
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
            z=0.d0
            DO j=1,6
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
               dg(j)=dg(j)+X(j)*F(i)
               DO k=j,7
                  ak(j,k)=ak(j,k)+X(j)*X(k)
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
C            call testreg(D,1)
            RETURN
         END IF
         IF(iter.gt.1.and.abs(rhon(1)*rhon(4)*rhon(6)).lt.1.d-10) THEN
C  prepare things for return if gradient is close to 0
            call regularD(D,negdefin)
C            call testreg(D,1)
C            call dblepr("rhon3",5,rhon,6)
C            call intpr("iter",4,iter,1)
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
               call rho2D(rhon,Dn)
               DO i=1,nb
                  z=0.d0
                  DO j=1,6
                     z=z+b(j,i)*Dn(j)
                  END DO
                  res=s(i)-th0n*exp(-z)
                  nrss=nrss+res*res
                  F(i)=res
               END DO
               crss=0.d0
               DO j=1,7
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
C      call testreg(D,2)
      RETURN
      END
