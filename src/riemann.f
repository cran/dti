CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   3D anisotropic smoothing of diffusion tensor data
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rawsdti(y,th,bi,ani,andir,det,bcov,sigma2,n1,n2,n3,h,
     1                  vext,rho,lambda,thnew,mask,wghts,yw,nw,rlm)
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
      integer n1,n2,n3,nw
      real*8 y(6,n1,n2,n3),th(6,n1,n2,n3),thnew(6,n1,n2,n3),h,rho,
     1       lambda,bi(n1,n2,n3),ani(n1,n2,n3),andir(3,n1,n2,n3),
     2       det(n1,n2,n3),bcov(6,6),sigma2(n1,n2,n3),vext(3),wghts(nw),
     3       yw(3,3,nw),thw(3,3),rem(3,3),rlm(3,3,nw)
      integer i1,j1,j1a,j1e,jj1,i2,j2,j2a,j2e,jj2,i3,j3,j3a,j3e,jj3,
     1        ierr,k,lw
      real*8 wij,adist,sw,swy(3,3),h3,thi(6),bii,sqrbii,ew(3),ev(3,3),
     1       mew,z1,z2,z3,dtidist2,sij,anii,deti,z,sew
      external adist,dtidist2
      logical aws,mask(n1,n2,n3) 
      aws=lambda.lt.1e20
      h3=h*h*h
C  now anisotropic smoothing 
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               lw=0
               if(mask(i1,i2,i3)) THEN
                  sw=0.d0
                  deti=exp(log(det(i1,i2,i3))/3)
                  bii=bi(i1,i2,i3)
                  sqrbii=sqrt(bii)*sigma2(i1,i2,i3)
                  thw(1,1)=th(1,i1,i2,i3)
                  thw(2,2)=th(4,i1,i2,i3)
                  thw(3,3)=th(6,i1,i2,i3)
                  thw(1,2)=th(2,i1,i2,i3)
                  thw(1,3)=th(3,i1,i2,i3)
                  thw(2,3)=th(5,i1,i2,i3)
                  DO k=1,6
                     thi(k)=th(k,i1,i2,i3)/deti
                  END DO
                  thi(1)=thi(1)+rho*sqrbii
                  thi(4)=thi(4)+rho*sqrbii
                  thi(6)=thi(6)+rho*sqrbii
C  this is scale invariant sice sqrbii scales with dsqrt(sigma2) (standard deviation)
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
                  anii=ani(i1,i2,i3)
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
                           if(wij.ge.h3) CYCLE
                           wij = (1.d0 - wij/h3)
                           IF(aws) THEN
                              sij=dtidist2(th(1,i1,i2,i3),
     1                             th(1,jj1,jj2,jj3),bcov)*bii/lambda
                              if(sij.gt.1.d0) CYCLE
                              wij=wij*(1.d0-sij)
                           END IF
                           lw=lw+1
                           IF(lw.gt.nw) THEN
                              call intpr("lw>nw",5,lw,1)
                              CYCLE
                           END IF
                           wij=wij/sigma2(jj1,jj2,jj3)
                           wghts(lw)=wij
                           yw(1,1,lw)=y(1,jj1,jj2,jj3)
                           yw(2,2,lw)=y(4,jj1,jj2,jj3)
                           yw(3,3,lw)=y(6,jj1,jj2,jj3)
                           yw(1,2,lw)=y(2,jj1,jj2,jj3)
                           yw(1,3,lw)=y(3,jj1,jj2,jj3)
                           yw(2,3,lw)=y(5,jj1,jj2,jj3)
                           yw(2,1,lw)=yw(1,2,lw)
                           yw(3,1,lw)=yw(1,3,lw)
                           yw(3,2,lw)=yw(2,3,lw)
                           sw=sw+wij
                        END DO
                     END DO
                  END DO
C
C    Compute the intrinsic mean of the Diffusion Tensors
C
                  IF(lw.gt.1) THEN
                     call rlogmap(thw,yw,lw,rlm)
                     DO j1=1,3
                        DO j2=j1,3
                           swy(j1,j2)=wghts(1)*rlm(j1,j2,1)
                        END DO
                     END DO
                     IF (lw.gt.1) THEN
                        DO k=2,lw
                           DO j1=1,3
                              DO j2=j1,3
                           swy(j1,j2)=swy(j1,j2)+wghts(k)*rlm(j1,j2,k)
                              END DO
                           END DO
                        END DO
                     END IF
                     DO j1=1,3
                        DO j2=j1,3
                           swy(j1,j2)=swy(j1,j2)/sw
                        END DO
                     END DO
                     swy(2,1)=swy(1,2)
                     swy(3,1)=swy(1,3)
                     swy(3,2)=swy(2,3)
                     call rexpmap(thw,swy,rem,ierr)
                     IF(ierr.eq.0) THEN
                        bi(i1,i2,i3)=sw
                        thnew(1,i1,i2,i3)=rem(1,1)
                        thnew(2,i1,i2,i3)=rem(1,2)
                        thnew(3,i1,i2,i3)=rem(1,3)
                        thnew(4,i1,i2,i3)=rem(2,2)
                        thnew(5,i1,i2,i3)=rem(2,3)
                        thnew(6,i1,i2,i3)=rem(3,3)
                     ELSE
                        DO k=1,6
                           thnew(k,i1,i2,i3)=th(k,i1,i2,i3)
                        END DO
                     END IF
                  ELSE
                     DO k=1,6
                        thnew(k,i1,i2,i3)=th(k,i1,i2,i3)
                     END DO
                  END IF
               ELSE
C   don't smooth, keep original data
                  bi(i1,i2,i3)=1.d0/sigma2(i1,i2,i3)
                  thnew(1,i1,i2,i3)=y(1,i1,i2,i3)
                  thnew(2,i1,i2,i3)=y(2,i1,i2,i3)
                  thnew(3,i1,i2,i3)=y(3,i1,i2,i3)
                  thnew(4,i1,i2,i3)=y(4,i1,i2,i3)
                  thnew(5,i1,i2,i3)=y(5,i1,i2,i3)
                  thnew(6,i1,i2,i3)=y(6,i1,i2,i3)
               END IF
               call eigen3(thnew(1,i1,i2,i3),ew,ev,ierr)
               if(ierr.ne.0) THEN
                  ani(i1,i2,i3)=0.d0
                  det(i1,i2,i3)=1
                  DO k=1,3
                     andir(k,i1,i2,i3)=0.d0
                  END DO
               ELSE
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
                  det(i1,i2,i3)=ew(1)*ew(2)*ew(3)
                  DO k=1,3
                     andir(k,i1,i2,i3)=ev(k,3)
                  END DO
               ENDIF
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute the Riemannian Exponential map Exp_{tensor}(tvect)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine rexpmap(tensor,tvect,rem,ierr)
      implicit logical (a-z)
      real*8 tensor(3,3),tvect(3,3),u(3,3),g(3,3),lam(3),y(3,3),
     1       v(3,3),sig(3),ginv(3,3),rtlam(3),w(3,3),rem(3,3),ts(3,3)
      integer i,j,k,ISUPPZ(6),iwork(50),ierr,nl
      real*8 work(104),vl,vu,z
      nl=3
C  get u and lambda
      DO i=1,3
         DO j=i,3
            ts(i,j)=tensor(i,j)
         END DO
      END DO
      call dsyevr('V','A','U',3,ts,3,vl,vu,1,3,1.d-10,nl,lam,
     1            u,3,ISUPPZ,work,104,iwork,50,ierr)
      if(ierr.ne.0) call intpr("ierr1",5,ierr,1)
      DO i=1,3
         rtlam(i)=sqrt(lam(i))
         DO j=1,3
            g(j,i)=u(j,i)*rtlam(i)
            ginv(i,j)=u(j,i)/rtlam(i)
         END DO
      END DO
C  get y
      DO i=1,3
         DO j=1,3
            z=0.d0
            DO k=1,3
               z=z+ginv(i,k)*tvect(k,j)
            END DO
            w(i,j)=z
         END DO
      END DO
      DO i=1,3
         DO j=i,3
            z=0.d0
            DO k=1,3
               z=z+w(i,k)*ginv(j,k)
            END DO
            y(i,j)=z
         END DO
      END DO
C     get v and sig
      call dsyevr('V','A','U',3,y,3,vl,vu,1,3,1.d-10,nl,sig,
     1            v,3,ISUPPZ,work,104,iwork,50,ierr)
C     get u=g%*%v 
      if(ierr.ne.0) call intpr("ierr2",5,ierr,1)
      DO i=1,3
         DO j=1,3
            z=0.d0
            DO k=1,3
               z=z+g(i,k)*v(k,j)
            END DO
            u(i,j)=z
         END DO
      END DO
C     get result in rem
      DO i=1,3
         if(sig(i).gt.20) THEN
            call dblepr("sig",3,sig,3)
            ierr=1
            RETURN
         ELSE
            sig(i)=exp(sig(i))
            ierr=0
         END IF
      END DO  
      DO i=1,3
         DO j=i,3
            z=0.d0
            DO k=1,3
               z=z+u(i,k)*u(j,k)*sig(k)
            END DO
            rem(i,j)=z
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute the Riemannian Log Map Log_{tensor}(tvect)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine rlogmap(tensor,tvect,n,rlm)
      implicit logical (a-z)
      integer n
      real*8 tensor(3,3),tvect(3,3,n),rlm(3,3,n)
      real*8 u(3,3),g(3,3),lam(3),y(3,3),ts(3,3),
     1       v(3,3),sig(3),ginv(3,3),rtlam(3),w(3,3)
      integer m,i,j,k,ISUPPZ(6),iwork(50),ierr,nl
      real*8 work(104),vl,vu,z
C  get u and lambda
      nl=3
      DO i=1,3
         DO j=i,3
            ts(i,j)=tensor(i,j)
         END DO
      END DO
      call dsyevr('V','A','U',3,ts,3,vl,vu,1,3,1.d-10,nl,lam,
     1            u,3,ISUPPZ,work,104,iwork,50,ierr)
      if(ierr.ne.0) call intpr("ierr3",5,ierr,1)
      if(nl.lt.3) call dblepr("rlp3nl",6,lam,3)
C  get g and ginv
      DO i=1,3
         rtlam(i)=sqrt(lam(i))
         DO j=1,3
            g(j,i)=u(j,i)*rtlam(i)
            ginv(i,j)=u(j,i)/rtlam(i)
         END DO
      END DO
      
C  thats the common part, now get results for all i = 1,n
C  get y
      DO m=1,n
         DO i=1,3
            DO j=1,3
               z=0.d0
               DO k=1,3
                  z=z+ginv(i,k)*tvect(k,j,m)
               END DO
               w(i,j)=z
            END DO
         END DO
         DO i=1,3
            DO j=i,3
               z=0.d0
               DO k=1,3
                  z=z+w(i,k)*ginv(j,k)
               END DO
               y(i,j)=z
            END DO
         END DO
C     get v and sig
         call dsyevr('V','A','U',3,y,3,vl,vu,1,3,1.d-10,nl,sig,
     1            v,3,ISUPPZ,work,104,iwork,50,ierr)
         if(ierr.ne.0) call intpr("ierr4",5,ierr,1)
         if(nl.lt.3) call dblepr("rlp4nl",6,sig,3)

C     get u=g%*%v 
         DO i=1,3
            DO j=1,3
               z=0.d0
               DO k=1,3
                  z=z+g(i,k)*v(k,j)
               END DO
               u(i,j)=z
            END DO
         END DO
C     get result in rlm
         DO i=1,3
            sig(i)=log(sig(i))
         END DO  
         DO i=1,3
            DO j=i,3
               z=0.d0
               DO k=1,3
                  z=z+u(i,k)*u(j,k)*sig(k)
               END DO
               rlm(i,j,m)=z
            END DO
         END DO
      END DO
      RETURN
      END
