CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   3D anisotropic smoothing of diffusion tensor data
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awssidti(s0,si,mask,th,bi,ani,andir,det,bcov,
     1                    solvebtb,sigma2,sigma2h,n1,n2,n3,ngrad,h,
     2                    vext,rho,lambda,thnew,sigma2n,swsi,eps)
C
C   y        -  observed diffusion tensor data
C   th       -  smoothed diffusion tensor data
C   bi       -  voxelwise sum of weights 
C   ani      -  anisotropy index 
C   dir      -  direction of main anisotropy 
C   det      -  det(A)
C   n1,n2,n3 -  spatial dimensions
C   rho      -  regularization parameter for anisotropic neighborhoods
C               (X,y,z) ( A(theta)+ rho/bi I ) (X,y,z)^T  = h^2  defines the elloispid 
C   lambda   -  scale factor in the statistical penalty
C   thnew    -  new smoothed diffusion tensor data
      implicit logical (a-z)
      integer n1,n2,n3,ngrad,s0(n1,n2,n3),si(n1,n2,n3,ngrad)
      real*8 th(6,n1,n2,n3),thnew(6,n1,n2,n3),h,rho,
     1       lambda,bi(n1,n2,n3),ani(n1,n2,n3),andir(3,n1,n2,n3),
     2       det(n1,n2,n3),bcov(6,6),sigma2(n1,n2,n3),vext(3),
     3       solvebtb(6,ngrad),swsi(ngrad),
     4       sigma2h(n1,n2,n3),sigma2n(n1,n2,n3)
      logical mask(n1,n2,n3)
      integer i1,j1,j1a,j1e,jj1,i2,j2,j2a,j2e,jj2,i3,j3,j3a,j3e,jj3,
     1        ierr,k,l
      real*8 wij,adist,sw,sws0,h2,thi(6),bii,sqrbii,ew(3),ev(3,3),
     1       mew,z1,z2,z3,dtidist2,sij,deti,z,sew,eps,eps3,ss2,sw0
      external adist,dtidist2
      logical aws
      aws=lambda.lt.1e20
      h2=h*h
      eps3=eps*eps*eps
C  now anisotropic smoothing 
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(.not.mask(i1,i2,i3)) CYCLE
               sw=0.d0
               sw0=0.d0
               sws0=0.d0
               ss2=0.d0
               DO k=1,ngrad
                  swsi(k)=0.d0
               END DO
               deti=exp(log(det(i1,i2,i3))/3)
               bii=bi(i1,i2,i3)
               sqrbii=sigma2h(i1,i2,i3)/sqrt(bii)
               DO k=1,6
                  thi(k)=th(k,i1,i2,i3)/deti
               END DO
               thi(1)=thi(1)+rho*sqrbii
               thi(4)=thi(4)+rho*sqrbii
               thi(6)=thi(6)+rho*sqrbii
C  this is scale invariant sice sqrbii scales with sqrt(sigma2) (standard deviation)
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
                        if(wij.gt.h2) CYCLE
                        wij = (1.d0 - wij/h2)
                        IF(aws) THEN
                           sij=dtidist2(th(1,i1,i2,i3),
     1                          th(1,jj1,jj2,jj3),bcov)*bii/lambda
                           if(sij.gt.1.d0) CYCLE
C  use Plateau kernel
                           if(sij.gt.0.5d0) THEN
                              wij=wij*2.d0*(1.d0-sij)
                           END IF
                        END IF
                        ss2=ss2+wij*sigma2(jj1,jj2,jj3)
                        sw0=sw0+wij
                        wij=wij/sigma2h(jj1,jj2,jj3)
                        sw=sw+wij
                        sws0=sws0+wij*s0(jj1,jj2,jj3)
                        DO k=1,ngrad
                           swsi(k)=swsi(k)+wij*si(jj1,jj2,jj3,k)
                        END DO
                     END DO
                  END DO
               END DO
               bi(i1,i2,i3)=sw
               IF(sws0.gt.0.d0) THEN
                  IF(sw.gt.0.d0.and.sw.lt.1.d20) THEN
                     sigma2n(i1,i2,i3)=ss2/sw0
                   ELSE
                     sigma2n(i1,i2,i3)=sigma2h(i1,i2,i3)
                  END IF
                  sws0=dlog(sws0)
                  DO k=1,ngrad
                     IF(swsi(k).gt.0.d0) THEN
                        swsi(k)=sws0-log(swsi(k))
                     ELSE
                        swsi(k)=0.d0
                     END IF
                  END DO
               ELSE
                  DO k=1,ngrad
                     swsi(k)=0.d0
                  END DO 
                  sigma2n(i1,i2,i3)=sigma2h(i1,i2,i3)
               END IF
               DO k=1,6
                  z=0.d0
                  DO l=1,ngrad
                     z=z+solvebtb(k,l)*swsi(l)
                  END DO
                  thnew(k,i1,i2,i3)=z
               END DO
               call eigen3(thnew(1,i1,i2,i3),ew,ev,ierr)
               IF(ierr.ne.0) THEN
C
C       error in eigenvalue decomposition of tensor, set voxel inactive
C
                  thnew(1,i1,i2,i3)=eps
                  thnew(2,i1,i2,i3)=0.d0
                  thnew(3,i1,i2,i3)=0.d0
                  thnew(4,i1,i2,i3)=eps
                  thnew(5,i1,i2,i3)=0.d0
                  thnew(6,i1,i2,i3)=eps
               ELSE IF(dmin1(ew(1),ew(2),ew(3)).lt.eps) THEN
C
C       negative eigenvalue of tensor, project to space of positive definite tensors
C
                  ew(1)=max(ew(1),eps)
                  ew(2)=max(ew(2),eps)
                  ew(3)=max(ew(3),eps)
                  thnew(1,i1,i2,i3)=ew(1)*ev(1,1)*ev(1,1)+
     1                     ew(2)*ev(1,2)*ev(1,2)+ew(3)*ev(1,3)*ev(1,3)
                  thnew(2,i1,i2,i3)=ew(1)*ev(1,1)*ev(2,1)+
     1                     ew(2)*ev(1,2)*ev(2,2)+ew(3)*ev(1,3)*ev(2,3)
                  thnew(3,i1,i2,i3)=ew(1)*ev(1,1)*ev(3,1)+
     1                     ew(2)*ev(1,2)*ev(3,2)+ew(3)*ev(1,3)*ev(3,3)
                  thnew(4,i1,i2,i3)=ew(1)*ev(2,1)*ev(2,1)+
     1                     ew(2)*ev(2,2)*ev(2,2)+ew(3)*ev(2,3)*ev(2,3)
                  thnew(5,i1,i2,i3)=ew(1)*ev(2,1)*ev(3,1)+
     1                     ew(2)*ev(2,2)*ev(3,2)+ew(3)*ev(2,3)*ev(3,3)
                  thnew(6,i1,i2,i3)=ew(1)*ev(3,1)*ev(3,1)+
     1                     ew(2)*ev(3,2)*ev(3,2)+ew(3)*ev(3,3)*ev(3,3)
               ELSE
C
C       valid tensor tensor, copy data
C
C                  DO k=1,6
C                     thnew(k,i1,i2,i3)=th(k,i1,i2,i3)
C                  END DO
               END IF
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize anisotropy index and direction of main anisotropy
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine projdt2(th,n1,n2,n3,thnew,ani,dir,det,eps)
C
C   th       -  observed diffusion tensor data
C   thnew    -  projected tensor data
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 th(6,n1,n2,n3),thnew(6,n1,n2,n3),ani(n1,n2,n3),
     1       dir(3,n1,n2,n3),det(n1,n2,n3),mew
      integer i1,i2,i3,ierr,k
      real*8 ew(3),ev(3,3),z1,z2,z3,z,eps,eps3
C  compute anisotropy index and direction of main anisotropy (nneded in statistical penalty)
      eps3=eps*eps*eps
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               call eigen3(th(1,i1,i2,i3),ew,ev,ierr)
               IF(ierr.ne.0) THEN
C
C       error in eigenvalue decomposition of tensor, set voxel inactive
C
                  thnew(1,i1,i2,i3)=eps
                  thnew(2,i1,i2,i3)=0.d0
                  thnew(3,i1,i2,i3)=0.d0
                  thnew(4,i1,i2,i3)=eps
                  thnew(5,i1,i2,i3)=0.d0
                  thnew(6,i1,i2,i3)=eps
C               ELSE IF(dmin1(ew(1),ew(2)).lt.dmax1(1.d-5,1d-5*ew(3))) 
C     1            THEN
               ELSE IF(min(ew(1),ew(2),ew(3)).lt.eps) THEN 
C
C       negative eigenvalue of tensor, project to space of positive definite tensors
C
C                  ew(1)=dmax1(ew(1),dmax1(1.d-5,1d-5*ew(3)))
C                  ew(2)=dmax1(ew(2),dmax1(1.d-5,1d-5*ew(3)))
                  ew(1)=max(eps,ew(1))
                  ew(2)=max(eps,ew(2))
                  ew(3)=max(eps,ew(3))
                  thnew(1,i1,i2,i3)=ew(1)*ev(1,1)*ev(1,1)+
     1                     ew(2)*ev(1,2)*ev(1,2)+ew(3)*ev(1,3)*ev(1,3)
                  thnew(2,i1,i2,i3)=ew(1)*ev(1,1)*ev(2,1)+
     1                     ew(2)*ev(1,2)*ev(2,2)+ew(3)*ev(1,3)*ev(2,3)
                  thnew(3,i1,i2,i3)=ew(1)*ev(1,1)*ev(3,1)+
     1                     ew(2)*ev(1,2)*ev(3,2)+ew(3)*ev(1,3)*ev(3,3)
                  thnew(4,i1,i2,i3)=ew(1)*ev(2,1)*ev(2,1)+
     1                     ew(2)*ev(2,2)*ev(2,2)+ew(3)*ev(2,3)*ev(2,3)
                  thnew(5,i1,i2,i3)=ew(1)*ev(2,1)*ev(3,1)+
     1                     ew(2)*ev(2,2)*ev(3,2)+ew(3)*ev(2,3)*ev(3,3)
                  thnew(6,i1,i2,i3)=ew(1)*ev(3,1)*ev(3,1)+
     1                     ew(2)*ev(3,2)*ev(3,2)+ew(3)*ev(3,3)*ev(3,3)
               ELSE
C
C       valid tensor tensor, copy data
C
                  DO k=1,6
                     thnew(k,i1,i2,i3)=th(k,i1,i2,i3)
                  END DO
               END IF
C
C       compute anisotropy index, direction corresponding to the first eigenvalue
C       and determinant
C
               mew=(ew(1)+ew(2)+ew(3))/3.d0
               z1=ew(1)-mew
               z2=ew(2)-mew
               z3=ew(3)-mew
               z=3.d0*(z1*z1+z2*z2+z3*z3)
               z1=ew(1)
               z2=ew(2)
               z3=ew(3)
               mew=2.d0*(z1*z1+z2*z2+z3*z3)
               if(mew.le.1d-20) mew=1.d0
               ani(i1,i2,i3)=sqrt(z/mew)
               DO k=1,3
                  dir(k,i1,i2,i3)=ev(k,3)
               END DO
               z=ew(1)*ew(2)*ew(3)
               IF(z.le.eps3) THEN
                  det(i1,i2,i3)=eps3
               ELSE
                  det(i1,i2,i3)=z
               END IF
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine smsigma(sigma2,n1,n2,n3,h,vext,sigma2h)
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 sigma2(n1,n2,n3),sigma2h(n1,n2,n3),vext(3),h
      integer i1,i2,i3,j1,j2,j3,ih1,ih2,ih3
      real*8 ssig,sw,w,h2,z1,z2,z3,z11,z22,z33
      h2=h*h
      ih1=h/vext(1)
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               sw=0.d0
               ssig=0.d0
               DO j1=i1-ih1,i1+ih1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  z1=(j1-i1)*vext(1)
                  z11=z1*z1
                  ih2=sqrt(h2-z11)/vext(2)
                  DO j2=i2-ih2,i2+ih2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     z2=(j2-i2)*vext(2)
                     z22=z11+z2*z2
                     ih3=sqrt(h2-z22)/vext(3)
                     DO j3=i3-ih3,i3+ih3
                        if(j3.lt.1.or.j3.gt.n3) CYCLE
                        z3=(j3-i3)*vext(3)
                        z33=z22+z3*z3                     
                        w=1.d0-z33/h2
                        sw=sw+w
                        ssig=ssig+w*sigma2(j1,j2,j3)
                     END DO
                  END DO
               END DO
               sigma2h(i1,i2,i3)=ssig/sw
            END DO
         END DO
      END DO
      RETURN
      END
