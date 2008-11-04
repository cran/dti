CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     regularize D 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine regD0(D,negdefin)
      implicit logical(a-z)
      logical negdefin
      real*8 D(6),ew(3),ev(3,3),ewmax
      integer ierr
      call eigen3(D,ew,ev,ierr)
      if(ew(1).le.1.d-4) THEN
C  first regularize
         negdefin=.TRUE.
         ewmax=max(1.d-4,(ew(1)+ew(2)+ew(3))/3.d0)
         D(1)=ewmax
         D(2)=0.d0
         D(3)=0.d0
         D(4)=ewmax
         D(5)=0.d0
         D(6)=ewmax
      ELSE
         negdefin=.FALSE.
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     regularize D 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine regularD(D,negdefin)
      implicit logical(a-z)
      logical negdefin
      real*8 D(6),ew(3),ev(3,3),ewmax
      integer ierr
      call eigen3(D,ew,ev,ierr)
      if(ew(1).le.1.d-4) THEN
C  first regularize
         negdefin=.TRUE.
         ewmax=max(1.d-4,ew(3))
         ew(1)=max(.1*ewmax,ew(1))
         ew(2)=max(.1*ewmax,ew(2))
         ew(3)=ewmax
         D(1)=ew(1)*ev(1,1)*ev(1,1)+ew(2)*ev(1,2)*ev(1,2)+
     1        ew(3)*ev(1,3)*ev(1,3)
         D(2)=ew(1)*ev(1,1)*ev(2,1)+ew(2)*ev(1,2)*ev(2,2)+
     1        ew(3)*ev(1,3)*ev(2,3)
         D(3)=ew(1)*ev(1,1)*ev(3,1)+ew(2)*ev(1,2)*ev(3,2)+
     1        ew(3)*ev(1,3)*ev(3,3)
         D(4)=ew(1)*ev(2,1)*ev(2,1)+ew(2)*ev(2,2)*ev(2,2)+
     1        ew(3)*ev(2,3)*ev(2,3)
         D(5)=ew(1)*ev(2,1)*ev(3,1)+ew(2)*ev(2,2)*ev(3,2)+
     1        ew(3)*ev(2,3)*ev(3,3)
         D(6)=ew(1)*ev(3,1)*ev(3,1)+ew(2)*ev(3,2)*ev(3,2)+
     1        ew(3)*ev(3,3)*ev(3,3)
      ELSE
         negdefin=.FALSE.
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     get rho from D 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine D2rho(D,rho)
      implicit logical(a-z)
      real*8 D(6),rho(6),eps
      eps=1.d-6
      rho(1)=sqrt(D(1)+eps)
      rho(2)=D(2)/rho(1)
      rho(3)=D(3)/rho(1)
      rho(4)=sqrt(D(4)-rho(2)*rho(2)+eps)
      rho(5)=(D(5)-rho(2)*rho(3))/rho(4)
      rho(6)=sqrt(D(6)-rho(3)*rho(3)-rho(5)*rho(5)+eps)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     get D from rho 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rho2D(rho,D)
      implicit logical(a-z)
      real*8 D(6),rho(6),eps
      eps=0.d-6
      D(1)=rho(1)*rho(1)+eps
      D(2)=rho(1)*rho(2)
      D(3)=rho(1)*rho(3)
      D(4)=rho(2)*rho(2)+rho(4)*rho(4)+eps
      D(5)=rho(2)*rho(3)+rho(4)*rho(5)
      D(6)=rho(3)*rho(3)+rho(5)*rho(5)+rho(6)*rho(6)+eps
C      call testreg2(D,rho,111)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     calculate fitted values of si
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sihat(th0i,Di,btb,F,nb)
      implicit logical (a-z)
      integer nb
      real*8 th0i,Di(6),btb(6,nb),F(nb)
      integer j,k
      real*8 z
      DO k=1,nb
         z=0.d0
         DO j=1,6
            z=z+btb(j,k)*Di(j)
         END DO
         F(k)=th0i*exp(-z)
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   3D anisotropic smoothing of diffusion tensor data
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsrgdti(si,siest,sipred,nb,n1,n2,n3,mask,sbind,btb,
     1                    sdcoef,th0,th0n,D,Dn,bi,ani,andir,det,
     2                    sigma2r,vol,niter,vext,rho0,lambda,swsi,
     3                    swsi2,F,varinv,eps,rician,nw,nriter,sisel,
     4                    isel,wselect,s2)
C
C   si       -  observed diffusion weighted images
C   nb       -  number of gradients (including zero gradients)
C   n1,n2,n3 -  dimensions of data cube
C   mask     -  logical: 1 for voxel inside region of interest
C   btb      -  b^T%*%b obtained from matrix of gradients b
C   sigma2   -  voxelwise varinviance estimates (data)
C   th0      -  current estimate of mean s_0 values
C   th0n     -  new  estimate of mean s_0 values (output)
C   D        -  current estimate of tensor
C   Dn       -  new estimate of tensor (output)
C   bi       -  sum of weights
C   ani      -  anisotropy index 
C   dir      -  direction of main anisotropy 
C   det      -  det(D)
C   h        -  actual bandwidth
C   niter    -  number of iterations for imbedded nl-regression
C   zext     -  
C   rho      -  regularization parameter for anisotropic neighborhoods
C               (X,y,z) ( A(theta)+ rho/bi I ) (X,y,z)^T  = h^2  defines the elloispid 
C   lambda   -  scale factor in the statistical penalty
C   swsi     -  auxiliary for sum_j w_j si(j)
C   F        -  auxiliary for function values in LSE
C   eps      -  something small and positive
      implicit logical (a-z)
      integer n1,n2,n3,nb,si(nb,n1,n2,n3),niter,siest(nb,n1,n2,n3),
     1       nw,nriter(nb),sisel(nb,nw),isel(3,nw)
      real*8 btb(6,nb),swsi2(nb),sdcoef(4),th0(n1,n2,n3),
     1       th0n(n1,n2,n3),sigma2r(n1,n2,n3),D(6,n1,n2,n3),
     2       Dn(6,n1,n2,n3),sipred(nb,n1,n2,n3),bi(n1,n2,n3),
     3       ani(n1,n2,n3),andir(3,n1,n2,n3),s2(nb),det(n1,n2,n3),h,
     4       rho0,vext(3),lambda,swsi(nb),F(nb),varinv(nb),eps,rss,
     5       wselect(nw),vol
      logical mask(n1,n2,n3),rician,sbind(nb)
      integer i1,j1,j1a,j1e,jj1,i2,j2,j2a,j2e,jj2,i3,j3,j3a,j3e,jj3,
     1        ierr,k,center,l,ns0
      real*8 wij,adist,sw,h2,thi(7),bii,sqrbii,ew(3),ev(3,3),rssj,
     1       mew,z1,z2,z3,sij,deti,z,sew,sw2,Di(6),dtidisrg,
     2       th0i,s2hat,ssigma2,rssi,h0,h1,squot,shat,low,up,zsd
      integer nselect
      real*8 x(10000),fw(10000)
      external adist,dtidisrg
      logical aws
      aws=lambda.lt.1e20
      if(rician) THEN
         call besselq(x,10000,fw)
      END IF
      h1=exp(log(vol*vext(1)*vext(2)*vext(3))/3.d0)
      h0=exp(log(vol*vext(1)*vext(2)*vext(3))/3.d0)/1.4
C  first fill predicted 
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               call sihat(th0(i1,i2,i3),D(1,i1,i2,i3),btb,
     1                    sipred(1,i1,i2,i3),nb)
            END DO
         END DO
      END DO
C  now anisotropic smoothing 
      ns0=0
      DO k=1,nb
         if(sbind(k)) THEN
            ns0=ns0+1
         END IF
      END DO
      low=sdcoef(1)+sdcoef(3)*sdcoef(2)
      up=sdcoef(1)+sdcoef(4)*sdcoef(2)
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(.not.mask(i1,i2,i3)) CYCLE
               DO k=1,nb
                  z=si(k,i1,i2,i3)
                  zsd=sdcoef(1)+z*sdcoef(2)
                  if(z.lt.sdcoef(3)) zsd=low
                  if(z.gt.sdcoef(4)) zsd=up
                  varinv(k)=1.d0/zsd/zsd
               END DO
               rssi = dtidisrg(siest(1,i1,i2,i3),
     1                         sipred(1,i1,i2,i3),varinv,nb)
               sw=0.d0
               sw2=0.d0
               DO k=1,nb
                  swsi(k)=0.d0
                  swsi2(k)=0.d0
               END DO
               deti=exp(log(det(i1,i2,i3))/3)
               bii=bi(i1,i2,i3)
               sqrbii=1.d0/sqrt(bii)
               th0i=th0(i1,i2,i3)
               th0n(i1,i2,i3)=th0i
C    used as initial values
               DO k=1,6
                  Di(k)=D(k,i1,i2,i3)
                  Dn(k,i1,i2,i3)=Di(k)
C    used as initial values
                  thi(k)=Di(k)/deti
               END DO
               thi(1)=thi(1)+rho0*sqrbii
               thi(4)=thi(4)+rho0*sqrbii
               thi(6)=thi(6)+rho0*sqrbii
C  this is scale invarinviant sice sqrbii scales with sqrt(sigma2) (standard deviation)
               call eigen3(thi,ew,ev,ierr)
               if(ierr.ne.0) THEN
                  call intpr("ierr",4,ierr,1)
                  thi(1)=1
                  thi(2)=0
                  thi(3)=0
                  thi(4)=1
                  thi(5)=0
                  thi(6)=1
               ELSE
                  if(ew(1).lt.1.d-5) THEN
                     ew(1)=max(1.d-6,ew(1))
                     ew(2)=max(1.d-6,ew(2))
                     ew(3)=max(1.d-6,ew(3))
                  END IF
                  sew=ew(1)*ew(2)*ew(3)
                  sew=exp(log(sew)/3.d0)
                  ew(1)=ew(1)/sew
                  ew(2)=ew(2)/sew
                  ew(3)=ew(3)/sew
                  thi(1)=ev(1,1)*ev(1,1)/ew(1)+ev(1,2)*ev(1,2)/ew(2)+
     1                   ev(1,3)*ev(1,3)/ew(3)
                  thi(2)=ev(1,1)*ev(2,1)/ew(1)+ev(1,2)*ev(2,2)/ew(2)+
     1                   ev(1,3)*ev(2,3)/ew(3)
                  thi(3)=ev(1,1)*ev(3,1)/ew(1)+ev(1,2)*ev(3,2)/ew(2)+
     1                   ev(1,3)*ev(3,3)/ew(3)
                  thi(4)=ev(2,1)*ev(2,1)/ew(1)+ev(2,2)*ev(2,2)/ew(2)+
     1                   ev(2,3)*ev(2,3)/ew(3)
                  thi(5)=ev(2,1)*ev(3,1)/ew(1)+ev(2,2)*ev(3,2)/ew(2)+
     1                   ev(2,3)*ev(3,3)/ew(3)
                  thi(6)=ev(3,1)*ev(3,1)/ew(1)+ev(3,2)*ev(3,2)/ew(2)+
     1                   ev(3,3)*ev(3,3)/ew(3)
               END IF
C   now get bandwidth such that the ellopsoid has specified volume
               call gethani(h0,h1,vol,thi,vext,1.d-2,h)
               h2=h*h
C   create needed estimates of s_i
               nselect=0
               bii=bii/lambda
               call rangex(thi,h,j1a,j1e,vext)
               DO j1=j1a,j1e
                  jj1=i1+j1
                  if(jj1.le.0.or.jj1.gt.n1) CYCLE
                  call rangey(thi,j1,h,j2a,j2e,vext)
                  DO j2=j2a,j2e
                     jj2=i2+j2
                     if(jj2.le.0.or.jj2.gt.n2) CYCLE
                     call rangez(thi,j1,j2,h,j3a,j3e,vext)
                     DO j3=j3a,j3e
                        jj3=i3+j3
                        if(jj3.le.0.or.jj3.gt.n3) CYCLE
                        if(.not.mask(jj1,jj2,jj3)) CYCLE
                        wij=adist(thi,j1,j2,j3,vext)
C     triangular location kernel
                        if(wij.ge.h2) CYCLE
                        center=abs(i1-jj1)+abs(i2-jj2)+abs(i3-jj3)
                        wij = (1.d0 - wij/h2)
                        IF(aws) THEN
                           rssj = dtidisrg(siest(1,i1,i2,i3),
     1                               sipred(1,jj1,jj2,jj3),varinv,nb)
                           sij = bii*max(0.d0,rssj-rssi)
                           if(center.eq.0) sij=0.d0
                           if(sij.gt.1.d0) CYCLE
                           wij=wij*(1.d0-sij)
                        END IF
                        sw=sw+wij
                        sw2=sw2+wij*wij
                        if(rician) THEN
                           DO k=1,nb
                              z=si(k,jj1,jj2,jj3)
                              swsi(k)=swsi(k)+wij*z
                              z=z*z
                              swsi2(k)=swsi2(k)+wij*z
                           END DO
                          if(nselect.gt.nw) THEN
                              call intpr("nselect>nw",10,nselect,1)
                              CYCLE
                           ENDIF 
                           nselect=nselect+1
                           isel(1,nselect)=jj1
                           isel(2,nselect)=jj2
                           isel(3,nselect)=jj3
                           wselect(nselect)=wij
                        ELSE
                           DO k=1,nb
                              swsi(k)=swsi(k)+wij*si(k,jj1,jj2,jj3)
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
               bi(i1,i2,i3)=sw
               Do k=1,nb
                  swsi(k)=swsi(k)/sw
               END DO                  
               if(rician.and.sw.gt.1.d0) THEN
                  ssigma2=0.d0
                  DO k=1,nb
                     if(sbind(k)) THEN
                        s2hat = swsi2(k)/sw-swsi(k)*swsi(k)
                        s2(k) = s2hat
                        ssigma2=ssigma2+s2hat
                        ns0=ns0+1
                     END IF
                  END DO
                  s2hat = sw*sw/(sw*sw-sw2)*ssigma2/ns0
                  shat = sqrt(s2hat)
C   this also adjusts for eliminating \theta by combining the second and 4th moment
                 DO k=1,nb
                     squot=swsi(k)/shat
                     if(squot.lt.1.d1) THEN
                        nriter(k)=6
                        if(squot.gt.2.25d0) nriter(k)=2
                        if(squot.gt.3.25d0) nriter(k)=1
                     ELSE
                        nriter(k)=0
                     END IF
                     DO l=1,nselect
                        sisel(k,l)=si(k,isel(1,l),isel(2,l),isel(3,l))
                     END DO
                  END DO
                  call ricecorr(sisel,wselect,nselect,nb,sbind,ns0,
     1                          nriter,sw,swsi,s2,s2hat,fw)
                  sigma2r(i1,i2,i3)=s2hat
               END IF
               IF(sw.gt.0.d0.and.sw.lt.1.d20) THEN
                  call dslvdti(swsi,nb,btb,sdcoef,varinv,th0n(i1,i2,i3),
     1                          Dn(1,i1,i2,i3),F,niter,eps,rss)
                  DO k=1,nb
                     siest(k,i1,i2,i3)=swsi(k)
                  END DO
               END IF
               call eigen3(Dn(1,i1,i2,i3),ew,ev,ierr)
               IF(ew(1).lt.0.d0) call dblepr("C0",2,ew,3)
               mew=(ew(1)+ew(2)+ew(3))/3.d0
               z1=ew(1)-mew
               z2=ew(2)-mew
               z3=ew(3)-mew
               z=3.d0*(z1*z1+z2*z2+z3*z3)
               z1=ew(1)
               z2=ew(2)
               z3=ew(3)
               mew=2.d0*(z1*z1+z2*z2+z3*z3)
               if(mew.le.eps*eps) mew=1.d0
               ani(i1,i2,i3)=sqrt(z/mew)
               det(i1,i2,i3)=ew(1)*ew(2)*ew(3)
               DO k=1,3
                  andir(k,i1,i2,i3)=ev(k,3)
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
      real*8 function dtidisrg(si,sj,varinv,nb)
      implicit logical (a-z)
      integer nb,si(nb)
      real*8 sj(nb),varinv(nb)
      real*8 z,z1
      integer i
      z=0.d0
      DO i=1,nb
         z1=si(i)-sj(i)
         z=z+z1*z1*varinv(i)
      END DO
      dtidisrg=z
      RETURN
      END
      subroutine islvdti(s,nb,b,sdcoef,varinv,th0,D,F,niter,eps,rss)
C
C  Implements the regularized Gauss-Newton Algortithm (10.2.8)
C  from Schwetlick (1979)
C  same as solvedti except that s is double
C
      implicit logical (a-z)
      integer nb,s(nb),niter
      real*8 D(6),b(6,nb),varinv(nb),th0,F(nb),eps,sdcoef(4)
      integer i,j,k,info,iter,icount
      logical negdefin
      real*8 z,gamma,alpha,delta,xzvarinv,
     1       dg(7),pk(7),ak(7,7),ck(7,7),rss,nrss,crss,maxabsdg,
     2       oldrss,relrss,Dn(6),res,X(7),th0n,zsd,low,up
C  first check if D defines a positive definite densor
      call regularD(D,negdefin)
      delta=0.25D0
      gamma=1.d0
      alpha=0.7d0
      oldrss=1.d50
      rss=0.d0
      low=sdcoef(1)+sdcoef(3)*sdcoef(2)
      up=sdcoef(1)+sdcoef(4)*sdcoef(2)
      DO i=1,nb
         zsd=sdcoef(1)+s(i)*sdcoef(2)
         if(s(i).lt.sdcoef(3)) zsd=low
         if(s(i).gt.sdcoef(4)) zsd=up
         varinv(i)=1.d0/zsd/zsd
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
               X(j)=b(j,i)*z
            END DO
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
            if(negdefin) THEN
C  estimate using reparametrization
               call islvdtir(s,nb,b,varinv,th0,D,F,niter,eps,rss)
            ENDIF
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
                  Dn(j)=D(j)-gamma*pk(j)
               END DO
               th0n=th0-gamma*pk(7)
               nrss=0.d0
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
C  check if tensor is positive definite
         DO j=1,6
            D(j)=Dn(j)
         END DO
         oldrss=rss
         rss=nrss
         call rchkusr()
      END DO
      call regularD(D,negdefin)
      if(negdefin) THEN
C  estimate using reparametrization
         call islvdtir(s,nb,b,varinv,th0,D,F,niter,eps,rss)
      ENDIF
      RETURN
      END
      subroutine dslvdti(s,nb,b,sdcoef,varinv,th0,D,F,niter,eps,rss)
C
C  Implements the regularized Gauss-Newton Algortithm (10.2.8)
C  from Schwetlick (1979)
C  same as solvedti except that s is double
C
      implicit logical (a-z)
      integer nb,niter
      real*8 s(nb),D(6),b(6,nb),varinv(nb),th0,F(nb),eps,sdcoef(4)
      integer i,j,k,info,iter,icount
      logical negdefin
      real*8 z,gamma,alpha,delta,xzvarinv,
     1       dg(7),pk(7),ak(7,7),ck(7,7),rss,nrss,crss,maxabsdg,
     2       oldrss,relrss,Dn(6),res,X(7),th0n,zsd,low,up
C  first check if D defines a positive definite densor
      call regularD(D,negdefin)
      delta=0.25D0
      gamma=1.d0
      alpha=0.7d0
      oldrss=1.d50
      rss=0.d0
      low=sdcoef(1)+sdcoef(3)*sdcoef(2)
      up=sdcoef(1)+sdcoef(4)*sdcoef(2)
      DO i=1,nb
         zsd=sdcoef(1)+s(i)*sdcoef(2)
         if(s(i).lt.sdcoef(3)) zsd=low
         if(s(i).gt.sdcoef(4)) zsd=up
         varinv(i)=1.d0/zsd/zsd
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
               X(j)=b(j,i)*z
            END DO
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
            if(negdefin) THEN
C  estimate using reparametrization
               call dslvdtir(s,nb,b,varinv,th0,D,F,niter,eps,rss)
            ENDIF
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
                  Dn(j)=D(j)-gamma*pk(j)
               END DO
               th0n=th0-gamma*pk(7)
               nrss=0.d0
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
C  check if tensor is positive definite
         DO j=1,6
            D(j)=Dn(j)
         END DO
         oldrss=rss
         rss=nrss
         call rchkusr()
      END DO
      call regularD(D,negdefin)
      if(negdefin) THEN
C  estimate using reparametrization
         call dslvdtir(s,nb,b,varinv,th0,D,F,niter,eps,rss)
C      else
C         call testreg(D,4)
      ENDIF
      RETURN
      END
