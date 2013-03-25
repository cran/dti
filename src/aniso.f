CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute all eigenvalues (lambda) of a 3x3 matrix and the corresponding EV (theta)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine eigen3(y,lambda,theta,ierr)
      implicit logical (a-z)
      integer ierr
      real*8 y(6),a(3,3),lambda(3),theta(3,3)
      integer i,j,l,ISUPPZ(6),lwork,iwork(50),liwork,n,m
      real*8 work(104),vl,vu,eps
      n=3
      m=3
      eps=1.d-50
      l=1
      DO i=1,3
         DO j=i,3
            a(i,j)=y(l)
            l=l+1
         END DO
      END DO
      lwork=104
      liwork=50
      call dsyevr('V','A','U',n,a,n,vl,vu,1,n,eps,m,lambda,
     1            theta,n,ISUPPZ,work,lwork,iwork,liwork,ierr)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute all eigenvalues (lambda) of a 3x3 matrix 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine eigen30(y,lambda,ierr)
      implicit logical (a-z)
      integer ierr
      real*8 y(6),a(3,3),lambda(3),theta(3,3)
      integer i,j,l,ISUPPZ(6),lwork,iwork(50),liwork,n,m
      real*8 work(104),vl,vu,eps
      n=3
      m=3
      eps=1.d-50
      l=1
      DO i=1,3
         DO j=i,3
            a(i,j)=y(l)
            l=l+1
         END DO
      END DO
      lwork=104
      liwork=50
      call dsyevr('N','A','U',n,a,n,vl,vu,1,n,eps,m,lambda,
     1            theta,n,ISUPPZ,work,lwork,iwork,liwork,ierr)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute anisotropic distance
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function adist(a,x,y,z,vext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer x,y,z
      real*8 a(6),xx,yy,zz,vext(3)
      xx=x*vext(1)
      yy=y*vext(2)
      zz=z*vext(3)      
      adist=a(1)*xx*xx+a(4)*yy*yy+a(6)*zz*zz+
     1               2.d0*(a(2)*xx*yy+a(3)*xx*zz+a(5)*yy*zz)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute statistical penalty for dti
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function dtidist2(th1,th2,Bcov)
      implicit logical (a-z)
      real*8 th1(6),th2(6),Bcov(6,6),z,zd,zd2
      integer i,j
      z=0.d0
      DO i=1,6
         zd=th1(i)-th2(i)
         z=z+zd*zd*Bcov(i,i)
         if(i.lt.6) THEN
            DO j=i+1,6
               zd2=th1(j)-th2(j)
               z=z+2.d0*zd*zd2*Bcov(i,j)
            END DO
         END IF
      END DO
      dtidist2=z
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute range of x for anisotropic neighborhood
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rangex(a,h,ia,ie,vext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer ia,ie
      real*8 a(6),h,vext(3)
      real*8 p1,p2,p3,p4,p5,p6,z,s,t,p55,p66,p44,h1
      p1=a(1)
      p2=a(4)
      p3=a(6)
      p4=a(2)
      p5=a(3)
      p6=a(5)
      p44=p4*p4
      p55=p5*p5
      p66=p6*p6
      s=p2*p3-p6*p6
      t=p5*p4/p2/p3
      z=(p1-p44/p2-p55/p3+2.d0*p6*t)*s-p66*p55/p3-p66*p44/p2+
     1                                             2.d0*p66*p6*t
      if(abs(z).le.1e-40) THEN
         z=1
      END IF
      z=s/z
      if(z.le.0) THEN
         z=0.d0
      ELSE 
         z=sqrt(z)
      END IF
      h1=h/vext(1)
      ia=-z*h1
      ie=z*h1
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute range of x for anisotropic neighborhood
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rangey(a,ix,h,ja,je,vext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer ja,je,ix
      real*8 a(6),h,vext(3)
      real*8 p1,p2,p3,p4,p5,p6,z,s,t,p55,p66,p44,x,h2s
      x=ix/h*vext(1)
      p1=a(1)
      p2=a(4)
      p3=a(6)
      p4=a(2)
      p5=a(3)
      p6=a(5)
      p44=p4*p4
      p55=p5*p5
      p66=p6*p6
      s=p2*p3-p6*p6
      if(s.le.1e-10) THEN
         z=1
      END IF
      t=-(p4*p3 - p6*p5)* x
      z=(p44*p3*p3-2.d0*p4*p3*p6*p5+p66*p55-p1*p3*s+p55*s)*x*x+p3*s
      if(z.le.0) THEN
         z=0.d0
      ELSE 
         z=sqrt(z)
      END IF
      h2s=h/vext(2)/s
      ja=(t-z)*h2s
      je=(t+z)*h2s
      z=z*h2s
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute range of x for anisotropic neighborhood
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rangez(a,ix,jy,h,ka,ke,vext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer ka,ke,ix,jy
      real*8 a(6),h,vext(3)
      real*8 p3,z,t,x,y,h3
      x=ix/h*vext(1)
      y=jy/h*vext(2)
      h3=h/vext(3)
      p3=a(6)
      z=(a(5)*y + a(3)*x)
      t= - z/p3
      z=(-z*t - a(4)*y*y - a(1)*x*x - 2.d0*a(2)*x*y + 1.d0)/p3
      if(z.le.0) THEN
         z=0.d0
      ELSE 
         z=sqrt(z)
      END IF
      ka=(t-z)*h3
      ke=(t+z)*h3
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute FA for a slice  (used in plot.r only)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dti2Dfa(D,n,mask,fa,md,adir)
      implicit logical (a-z)
      integer n
      logical mask(n)
      real*8 D(6,n),fa(n),md(n),adir(3,n)
      integer i,ierr
      real*8 lambda(3),evec(3,3),trc,d1,d2,d3,a1,a2,a3,dd
      DO i=1,n
         if(mask(i)) THEN
            call eigen3(D(1,i),lambda,evec,ierr)
            a1=lambda(1)
            a2=lambda(2)
            a3=lambda(3)
            trc=(a1+a2+a3)/3.d0
            adir(1,i)=evec(1,3)
            adir(2,i)=evec(2,3)
            adir(3,i)=evec(3,3)
            md(i)=trc
            d1=a1-trc
            d2=a2-trc
            d3=a3-trc
            dd=a1*a1+a2*a2+a3*a3
            IF(dd.gt.1.d-12) THEN
               fa(i)=sqrt(1.5d0*(d1*d1+d2*d2+d3*d3)/dd)
            ELSE
               fa(i)=0.d0
            ENDIF
         ELSE
            md(i)=0.d0
            fa(i)=0.d0
            adir(1,i)=1.d0
            adir(2,i)=0.d0
            adir(3,i)=0.d0
         END IF
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute GA for a slice  (used in plot.r only)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dti2Dga(D,n,mask,ga,md,adir)
      implicit logical (a-z)
      integer n
      logical mask(n)
      real*8 D(6,n),ga(n),md(n),adir(3,n)
      integer i,ierr
      real*8 lambda(3),evec(3,3),trc,d1,d2,d3,a1,a2,a3,dd
      DO i=1,n
         if(mask(i)) THEN
            call eigen3(D(1,i),lambda,evec,ierr)
            a1=log(lambda(1))
            a2=log(lambda(2))
            a3=log(lambda(3))
            trc=(a1+a2+a3)/3.d0
            adir(1,i)=evec(1,3)
            adir(2,i)=evec(2,3)
            adir(3,i)=evec(3,3)
            md(i)=trc
            d1=a1-trc
            d2=a2-trc
            d3=a3-trc
            dd=a1*a1+a2*a2+a3*a3
            IF(dd.gt.1.d-12) THEN
               ga(i)=sqrt(1.5d0*(d1*d1+d2*d2+d3*d3)/dd)
            ELSE
               ga(i)=0.d0
            ENDIF
         ELSE
            md(i)=0.d0
            ga(i)=0.d0
            adir(1,i)=1.d0
            adir(2,i)=0.d0
            adir(3,i)=0.d0
         END IF
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute DTI-Indices for a volume
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dtiind3D(D,n,fa,ga,md,adir,bary)
      implicit logical (a-z)
      integer n
      real*8 D(6,n),fa(n),md(n),adir(3,n),bary(3,n),ga(n)
      integer i,ierr
      real*8 lambda(3),evec(9),trc,d1,d2,d3,a1,a2,a3,dd
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(D,n,fa,ga,md,adir,bary)
C$OMP& PRIVATE(i,ierr,lambda,evec,trc,d1,d2,d3,a1,a2,a3,dd)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         call eigen3(D(1,i),lambda,evec,ierr)
         a1=lambda(1)
         a2=lambda(2)
         a3=lambda(3)
         trc=(a1+a2+a3)/3.d0
         md(i)=trc
         d1=a1-trc
         d2=a2-trc
         d3=a3-trc
         dd=a1*a1+a2*a2+a3*a3
         IF(dd.gt.1.d-12) THEN
            fa(i)=sqrt(1.5d0*(d1*d1+d2*d2+d3*d3)/dd)
            bary(1,i)=(a3-a2)/trc/3.d0
            bary(2,i)=2.d0*(a2-a1)/trc/3.d0
            bary(3,i)=a1/trc
         ELSE
            fa(i)=0.d0
            bary(1,i)=0.d0
            bary(2,i)=0.d0
            bary(3,i)=1.d0
         ENDIF
         d1=log(a1)
         d2=log(a2)
         d3=log(a3)
         dd=(d1+d2+d3)/3.d0
         d1=d1-dd
         d2=d2-dd
         d3=d3-dd
         ga(i)=sqrt(d1*d1+d2*d2+d3*d3)
         adir(1,i)=evec(7)
         adir(2,i)=evec(8)
         adir(3,i)=evec(9)
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(fa,ga,md,adir,bary)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute subset of DTI-Indices for a volume
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dtieigen(D,n,fa,ev,adir)
      implicit logical (a-z)
      integer n
      real*8 D(6,n),fa(n),ev(3,n),adir(6,n)
      integer i,ierr
      real*8 lambda(3),evec(9),trc,d1,d2,d3,a1,a2,a3,dd,fai
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(D,n,fa,ev,adir)
C$OMP& PRIVATE(i,ierr,lambda,evec,trc,d1,d2,d3,a1,a2,a3,dd,fai)
C$OMP DO SCHEDULE(STATIC)
      DO i=1,n
         call eigen3(D(1,i),lambda,evec,ierr)
         a1=lambda(1)
         a2=lambda(2)
         a3=lambda(3)
         trc=(a1+a2+a3)/3.d0
         d1=a1-trc
         d2=a2-trc
         d3=a3-trc
         dd=a1*a1+a2*a2+a3*a3
         IF(dd.gt.1.d-12) THEN
            fai=sqrt(1.5d0*(d1*d1+d2*d2+d3*d3)/dd)
         ELSE
            fai=0.d0
         ENDIF
         adir(1,i)=evec(7)
         adir(2,i)=evec(8)
         adir(3,i)=evec(9)
         adir(4,i)=evec(4)
         adir(5,i)=evec(5)
         adir(6,i)=evec(6)
         ev(1,i)=a3
         ev(2,i)=a2
         ev(3,i)=a1
         fa(i)=fai
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(fa,adir,ev)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute DTI-Indices for a volume
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dti3Dall(D,n,fa,ga,md,adir,ev)
      implicit logical (a-z)
      integer n
      real*8 D(6,n),fa(n),md(n),adir(3,n),ev(3,n),ga(n)
      integer i,ierr
      real*8 evec(3,3),trc,d1,d2,d3,a1,a2,a3,dd
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(D,n,mask,fa,ga,md,adir,ev)
C$OMP& PRIVATE(i,ierr,evec,trc,d1,d2,d3,a1,a2,a3,dd)
C$OMP DO SCHEDULE(STATIC)
      DO i=1,n
         call eigen3(D(1,i),ev(1,i),evec,ierr)
         a1=ev(1,i)
         a2=ev(2,i)
         a3=ev(3,i)
         trc=(a1+a2+a3)/3.d0
         adir(1,i)=evec(1,3)
         adir(2,i)=evec(2,3)
         adir(3,i)=evec(3,3)
         md(i)=trc
         d1=a1-trc
         d2=a2-trc
         d3=a3-trc
         dd=a1*a1+a2*a2+a3*a3
         IF(dd.gt.1.d-12) THEN
            fa(i)=sqrt(1.5d0*(d1*d1+d2*d2+d3*d3)/dd)
         ELSE
            fa(i)=0.d0
            ev(1,i)=0.d0
            ev(2,i)=0.d0
            ev(3,i)=0.d0
         ENDIF
         d1=log(a1)
         d2=log(a2)
         d3=log(a3)
         dd=(d1+d2+d3)/3.d0
         d1=d1-dd
         d2=d2-dd
         d3=d3-dd
         ga(i)=sqrt(d1*d1+d2*d2+d3*d3)
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(fa,ga,md,adir,ev)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute DTI-eigenvalues for a volume
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dti3Dev(D,n,ev)
      implicit logical (a-z)
      integer n
      real*8 D(6,n),ev(3,n)
      integer i,ierr
      DO i=1,n
            call eigen30(D(1,i),ev(1,i),ierr)
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute DTI-eigenvectors for a volume 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dti3Dand(D,n,andir)
      implicit logical (a-z)
      integer n
      real*8 D(6,n),andir(3,n)
      integer i,ierr
      real*8 lambda(3),evec(3,3)
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(D,n,andir)
C$OMP& PRIVATE(i,lambda,evec,ierr)
C$OMP DO SCHEDULE(STATIC)
      DO i=1,n
         call eigen3(D(1,i),lambda,evec,ierr)
         andir(1,i)=evec(1,3)
         andir(2,i)=evec(2,3)
         andir(3,i)=evec(3,3)
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(andir)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute all(!) DTI-eigenvectors for a volume 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dti3DevAll( D, n, andir, evalues)
      implicit logical (a-z)
      integer n
      real*8 D( 6, n), andir( 9, n), evalues( 3, n)
      integer i, ierr
      real*8 lambda( 3), evec( 3, 3)
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(D,n,andir,evalues)
C$OMP& PRIVATE(i,lambda,evec,ierr)
C$OMP DO SCHEDULE(STATIC)
      DO i=1,n
         call eigen3(D(1,i),lambda,evec,ierr)
         andir( 1, i) = evec( 1, 3)
         andir( 2, i) = evec( 2, 3)
         andir( 3, i) = evec( 3, 3)
         andir( 4, i) = evec( 1, 2)
         andir( 5, i) = evec( 2, 2)
         andir( 6, i) = evec( 3, 2)
         andir( 7, i) = evec( 1, 1)
         andir( 8, i) = evec( 2, 1)
         andir( 9, i) = evec( 3, 1)
         evalues( 1, i) = lambda( 3)
         evalues( 2, i) = lambda( 2)
         evalues( 3, i) = lambda( 1)
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(andir)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Regularize tensors bu setting eigenvalues to zero
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dti3Dreg(D,n)
      implicit logical (a-z)
      integer n
      real*8 D(6,n)
      integer i,ierr
      real*8 lam(3),evec(3,3)
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(D,n)
C$OMP& PRIVATE(i,lam,evec,ierr)
C$OMP DO SCHEDULE(STATIC)
      DO i=1,n
         call eigen3(D(1,i),lam,evec,ierr)
         lam(1)=max(0.d0,lam(1))
         lam(2)=max(0.d0,lam(2))
         D(1,i) = lam(1)*evec(1,1)*evec(1,1)+
     1          lam(2)*evec(1,2)*evec(1,2)+lam(3)*evec(1,3)*evec(1,3)
         D(2,i) = lam(1)*evec(1,1)*evec(2,1)+
     1          lam(2)*evec(1,2)*evec(2,2)+lam(3)*evec(1,3)*evec(2,3)
         D(3,i) = lam(1)*evec(1,1)*evec(3,1)+
     1          lam(2)*evec(1,2)*evec(3,2)+lam(3)*evec(1,3)*evec(3,3)
         D(4,i) = lam(1)*evec(2,1)*evec(2,1)+
     1          lam(2)*evec(2,2)*evec(2,2)+lam(3)*evec(2,3)*evec(2,3)
         D(5,i) = lam(1)*evec(2,1)*evec(3,1)+
     1          lam(2)*evec(2,2)*evec(3,2)+lam(3)*evec(2,3)*evec(3,3)
         D(6,i) = lam(1)*evec(3,1)*evec(3,1)+
     1          lam(2)*evec(3,2)*evec(3,2)+lam(3)*evec(3,3)*evec(3,3)
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(D)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   determine sum of location weights for a given geometry a(6) and given 
C   bandwidth
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Algorithmus zur Nullstellenbestimmung einer monotonen Funktion auf(0,\infty)
      subroutine gethani(x,y,value,a,vext,eps,bw)
      implicit logical(a-z)
      real*8 x,y,value,a(6),vext(3),eps,bw
      real*8 fw1,fw2,fw3,z
      real*8 sofw3D
      external sofw3D
      if(x.ge.y) RETURN
      fw1=sofw3D(a,x,vext)
      fw2=sofw3D(a,y,vext)
      DO WHILE(fw1.gt.value)
         x=x*x/y
         fw1=sofw3D(a,x,vext)
      END DO
      DO WHILE(fw2.le.value)
         y=y*y/x
         fw2=sofw3D(a,y,vext)
      END DO
      DO WHILE(min(fw2/value,value/fw1).gt.1.d0+eps)
         z=x+(value-fw1)/(fw2-fw1)*(y-x)
         fw3=sofw3D(a,z,vext)
         if(fw3.le.value) THEN
            x=z
            fw1=fw3
         ENDIF
         if(fw3.ge.value) THEN
            y=z
            fw2=fw3
         ENDIF
               call rchkusr()
      END DO
      if(fw2/value.gt.value/fw1) THEN
          bw=x+(value-fw1)/(fw2-fw1)*(y-x)
      ELSE
          bw=y-(fw2-value)/(fw2-fw1)*(y-x)
      ENDIF
      RETURN
      END  
      subroutine getvofh(a,bw,vext,vol)
      implicit logical(a-z)
      real*8 a(6),bw,vext(3),vol,sofw3D
      vol=sofw3D(a,bw,vext)
      RETURN
      END
C  compute sum of weights for anisotropic smoothing      
      real*8 function sofw3D(a,bw,vext)
      implicit logical(a-z)
      real*8 a(6),bw,vext(3)
      integer ia1,ie1,ia2,ie2,ia3,ie3,i1,i2,i3
      real*8 wij,sw,h2,adist
      external adist
      h2=bw*bw
      sw=1.d0
      call rangex(a,bw,ia1,ie1,vext)
      DO i1=1,ie1
         call rangey(a,i1,bw,ia2,ie2,vext)
         DO i2=ia2,ie2
            call rangez(a,i1,i2,bw,ia3,ie3,vext)
            DO i3=ia3,ie3
               wij=max(0.d0,1.d0-adist(a,i1,i2,i3,vext)/h2)
               sw=sw+2.d0*wij
            END DO
         END DO
      END DO
C  now case i1=0
      call rangey(a,0,bw,ia2,ie2,vext)
      DO i2=1,ie2
         call rangez(a,0,i2,bw,ia3,ie3,vext)
         DO i3=ia3,ie3
            wij=max(0.d0,1.d0-adist(a,0,i2,i3,vext)/h2)
            sw=sw+2.d0*wij
         END DO
      END DO
C  now case i1=i2=0
      call rangez(a,0,0,bw,ia3,ie3,vext)
         DO i3=1,ie3
            wij=max(0.d0,1.d0-adist(a,0,0,i3,vext)/h2)
            sw=sw+2.d0*wij
         END DO
      sofw3D=sw
      RETURN
      END
