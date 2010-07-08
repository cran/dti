C
C __________________________________________________________________
C
      subroutine mfunpl0(par,siq,grad,m,lpar,ngrad,pen,z,w,erg)
C
C   model without isotropic compartment 
C   same as mfunpl but with unconstrained least squares for weights
C
C   code is restricted to m<=6
C
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad,m),erg,pen
      integer i,j,i3,mode,r
      real*8 th,w(ngrad),sw,sth,z1,p0,p1,d1,d2,d3,work(1000),s(6)
      th = max(par(1),-5.d0)
      DO i = 1,m
         i3=2*i
         p0=par(i3)
         p1=par(i3+1)
         sth = sin(p0)
         d1 = sth*cos(p1)
         d2 = sth*sin(p1)
         d3 = cos(p0)
         DO j = 1,ngrad
            z1 = d1*grad(1,j)+d2*grad(2,j)+d3*grad(3,j)
            z(j,i) = exp(-th*z1*z1)
         END DO
      END DO
C  
C    siq will be replaced, need to copy it if C-version of optim is used
C
      call dcopy(ngrad,siq,1,w,1)
      call dgelss(ngrad,m,1,z,ngrad,w,ngrad,s,-1.d0,r,work,1000,mode)
      IF(mode.ne.0) THEN
         call intpr("mode",4,mode,1)
         erg = 1d20
      ELSE
         sw=0.d0
C penalize for extreme th values
         if(th.gt.8.d0) sw=sw+th-8.d0
C penalize for negative weights
         if(th.lt.1.d-2) sw=sw-1.d2*th+1.d0
         DO i=1,m
            if(w(i).lt.0.d0) sw=sw-pen*w(i)
         END DO
         DO i=m+1,ngrad
            sw=sw+w(i)**2
         END DO
         erg=sw
      END IF
      call rchkusr()
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine mfunpl0g(par,s,g,m,lpar,n,d,z,v,w,dkgj,dkgj2,
     1                    ddkdphig,ddkdetag,dddphi,
     2                    dddeta,dvdth,dvdphi,dvdeta,dzdpars,
     3                    dwdpars,zs,work1,
     3                    work2,scopy,pen,dfdpar)
      implicit logical (a-z)
      integer m,n,lpar
      real*8 par(lpar),s(n),g(3,n),d(3,m),z(n,m),v(m,m),dkgj(n,m),
     1       w(m),dkgj2(n,m),ddkdphig(n,m),ddkdetag(n,m),dddphi(3,m),
     2       dddeta(3,m),dvdth(m,m),dvdphi(m,m,m),dvdeta(m,m,m),
     3       dzdpars(n,m,lpar),dwdpars(m,lpar),zs(n,m),dfdpar(lpar),pen
      integer i,j,k,l,i3,ind(10),mode
      real*8 th,sw,sphi,cphi,seta,ceta,z1,z2,p0,p1,
     1       work(1000),work1(n,m),work2(n,m),scopy(n)
      real*8 ddot
      external ddot
      th = max(-5.d0,par(1))
      sw = 0
C
C  get d, dkgj, dkgj2, dddphi, dddeta, z, ddkdphig, ddkdetag
C
      DO i = 1,m
         i3=2*i
         p0=par(i3)
         p1=par(i3+1)
         sphi = sin(p0)
         cphi = cos(p0)
         seta = sin(p1)
         ceta = cos(p1)
         d(1,i) = sphi*ceta
         d(2,i) = sphi*seta
         d(3,i) = cphi
         dddphi(1,i) = cphi*ceta
         dddphi(2,i) = cphi*seta
         dddphi(3,i) = -sphi
         dddeta(1,i) = -sphi*seta
         dddeta(2,i) = sphi*ceta
         dddeta(3,i) = 0.d0
         DO j = 1,n
            z1 = ddot(3,d(1,i),1,g(1,j),1)
            dkgj(j,i) = z1
            z2 = z1*z1
            dkgj2(j,i) = z2
            z(j,i) = exp(-th*z2)
            zs(j,i) = z(j,i)*s(j)
            ddkdphig(j,i) = ddot(3,dddphi(1,i),1,g(1,j),1)
            ddkdetag(j,i) = ddot(3,dddeta(1,i),1,g(1,j),1)
         END DO
      END DO
C  
C   we now have d, dkgj,dddphi, dddeta, z, ddkdphig, ddkdetag
C   next w
C
      call dcopy(n,s,1,scopy,1)
      call dcopy(n*m,z,1,work1,1)
      call dgelss(n,m,1,work1,n,scopy,n,w,-1.d0,r,work,1000,mode)
      IF(mode.gt.1) THEN
         call intpr("mode",4,mode,1)
      END IF
      call dcopy(m,scopy,1,w,1)
C
C   thats weights in w now V, dVdth, dVdphi, dVdeta, dzdpars
C
C   use work1, work2 and scopy for intermediate results
      DO k=1,m
         DO j=1,n
            work1(j,k) = dkgj(j,k)*ddkdphig(j,k)
            work2(j,k) = dkgj(j,k)*ddkdetag(j,k)
         END DO
      END DO
C initialize unneeded elements
      DO k=1,m
         DO l=1,m
            DO j=1,n
               dzdpars(j,k,1+l)=0.d0
            END DO
         END DO
      END DO
      DO k=1,m
         DO j=1,n
            dzdpars(j,k,1)=-dkgj2(j,k)*z(j,k)
            dzdpars(j,k,1+k)=-2.d0*th*work1(j,k)*z(j,k)
            dzdpars(j,k,1+m+k)=-2.d0*th*work2(j,k)*z(j,k)
         END DO
      END DO
      DO k=1,m
         v(k,k)=ddot(n,z(1,k),1,z(1,k),1)
         z1 = ddot(n,dzdpars(1,k,1),1,z(1,k),1)
         dVdth(k,k) = 2.d0*z1
         DO l=k+1,m
            z2=z1+ddot(n,dzdpars(1,l,1),1,z(1,l),1)
            dVdth(k,l) = z2
            dVdth(l,k) = z2
         END DO
         DO i=1,m
            DO l=1,m
               dVdphi(i,l,k) = 0.d0
               dVdeta(i,l,k) = 0.d0
            END DO
         END DO
         DO i=1,m
            dVdphi(i,k,k) = dVdphi(i,k,k)+
     1                      ddot(n,dzdpars(1,k,1+k),1,z(1,i),1)
            dVdphi(k,i,k) = dVdphi(i,k,k)+
     1                      ddot(n,dzdpars(1,k,1+k),1,z(1,i),1)
            dVdeta(i,k,k) = dVdeta(i,k,k)+
     1                      ddot(n,dzdpars(1,k,1+m+k),1,z(1,i),1)
            dVdeta(k,i,k) = dVdeta(i,k,k)+
     1                      ddot(n,dzdpars(1,k,1+m+k),1,z(1,i),1)
         END DO
      END DO
C
C   thats V, dVdth, dVdphi, dVdeta now fill dwdpars)
C
      DO k=1,m
         dwdpars(k,1) = ddot(n,dzdpars(1,k,1),1,s,1) -
     1                         ddot(m,dVdth(1,k),1,w,1)
         DO l=1,m
            dwdpars(k,1+l) = ddot(n,dzdpars(1,k,1+l),1,s,1)  - 
     1                              ddot(m,dVdphi(1,k,l),1,w,1)
            dwdpars(k,1+m+l) = ddot(n,dzdpars(1,k,1+m+l),1,s,1)  -
     1                              ddot(m,dVdeta(1,k,l),1,w,1)
         END DO
      END DO
C
C   thats  dzdpars  now compute  dw/dpar in dwdpars
C
      call dsysv("U",m,lpar,v,m,ind,dwdpars,m,work,1000,mode)
      IF(mode.ne.0) THEN
         call intpr("mode2",5,mode,1)
      END IF
C
C   now we have dw/dpar in dzdpars next residuals in scopy
C
      DO j=1,n
         z1 = s(j)
         DO k=1,m
            z1 = z1 - w(k)*z(j,k)
         END DO
         scopy(j) = z1
      END DO
C
C   now we have residuals in scopy compute gradient of f
C   use work for intermediate results
      DO j = 1,n
         z1 = 0.d0
         DO k=1,m
            z1=z1+w(k)*dzdpars(j,k,1)+dwdpars(k,1)*z(j,k)
         END DO
         work(j)=z1
      END DO
      dfdpar(1)=-2.d0*ddot(n,work,1,scopy,1)
      if(th.gt.8) dfdpar(1)=dfdpar(1)+1.d0
      if(th.lt.1d-2) dfdpar(1)=dfdpar(1)-1.d2
      DO k=1,m
         if(w(k).lt.0.d0) dfdpar(1)=dfdpar(1)-pen*dwdpars(k,1)
      END DO
C 
C    dzdpars contains dw/dpars 
C
      DO i = 1,m
         DO j = 1,n
            z1=w(i)*dzdpars(j,i,1+i)
            DO k=1,m
               z1=z1+dwdpars(k,1+i)*z(j,k)
            END DO
            work(j)=z1
         END DO
         dfdpar(1+2*i-1)=-2.d0*ddot(n,work,1,scopy,1)
         DO k=1,m
            if(w(k).lt.0) dfdpar(1+2*i-1)=dfdpar(1+2*i-1)-
     1                                    pen*dwdpars(k,1+i)
         END DO
         DO j = 1,n
            z1=w(i)*dzdpars(j,i,1+m+i)
            DO k=1,m
               z1=z1+dwdpars(k,1+m+i)*z(j,k)
            END DO
            work(j)=z1
         END DO
         dfdpar(1+2*i)=-2.d0*ddot(n,work,1,scopy,1)
         DO k=1,m
            if(w(k).lt.0) dfdpar(1+2*i)=dfdpar(1+2*i)-
     1                                pen*dwdpars(k,1+m+i)
         END DO
      END DO
C
C   thats derivative with respect to theta
C
      call rchkusr()
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine getsiin2(si,ngrad,n1,n2,n3,m,
     1         egrad,isample,ntry,sms,z,siind,mval,ns,mask)
C
C  compute diagnostics for initial estimates in siind
C  siind(1,i1,i2,i3) will contain the model order 
C  
C  si     - array of si-values
C  m      - model order
C  maxc   - maximum of cos(angle between directions)
C  exgrad - exp(-theta1*dgrad^2) 
C  isample - guesses for gradient directions
C  ntry   - number of guesses
C  sms    - copies of si
C  z      - array for design matrix corresponding to guesses
C  siind  - array of indices (output)
C  ns     - m+1
C  mask   - mask
C
C  restricted to ngrad<=1000 and m <=10
C
      implicit logical (a-z)
      integer n1,n2,n3,ngrad,ns,siind(ns,n1,n2,n3),m,ntry,
     1       isample(m,ntry)
      real*8 si(ngrad,n1,n2,n3),sms(ngrad),
     1       egrad(ngrad,ngrad),z(ngrad,ns),mval(n1,n2,n3)
      logical mask(n1,n2,n3)
      integer i1,i2,i3,k,ibest,mode,ind(10),l
      real*8 w(1000),krit,work1(1000),work2(10),erg
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(mask(i1,i2,i3)) THEN
C  now search for minima of sms (or weighted sms
                  ibest=1
                  krit=1e10
                  DO k=1,ntry
                     call dcopy(ngrad,si(1,i1,i2,i3),1,sms,1)
                     DO l=1,m
                  call dcopy(ngrad,egrad(1,isample(l,k)),1,z(1,l),1)
                     END DO
        call nnls(z,ngrad,ngrad,m,sms,w,erg,work2,work1,ind,mode)
                     IF(mode.gt.1) THEN
                        call intpr("mode",4,mode,1)
                        call intpr("isample",7,isample(1,k),m)
                     ELSE 
                        IF(erg.lt.krit) THEN
                           krit=erg
                           ibest=k
                        END IF  
                     END IF
                  END DO
                  siind(1,i1,i2,i3)=m
                  DO l=1,m
                     siind(l+1,i1,i2,i3)=isample(l,ibest)
                  END DO
                  mval(i1,i2,i3)=krit
               ELSE
                  siind(1,i1,i2,i3)=-1
                  mval(i1,i2,i3)=0
               END IF
            END DO
         END DO
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine getsiind(si,ngrad,n1,n2,n3,m,maxc,dgrad,
     1                    sms,siind,ns,mask)
C
C  compute diagnostics for initial estimates in siind
C  siind(1,i1,i2,i3) will contain the model order 
C  
C  si     - array of si-values
C  m      - model order
C  maxc   - maximum of cos(angle between directions)
C  dgrad  - matrix of pairwise cos of gradient directions
C  sms    - copies of si
C  siind  - array of indices (output)
C
      implicit logical (a-z)
      integer n1,n2,n3,ngrad,ns,siind(ns,n1,n2,n3),m
      real*8 si(ngrad,n1,n2,n3),sms(ngrad),dgrad(ngrad,ngrad),maxc
      logical mask(n1,n2,n3)
      integer i,i1,i2,i3,imin,k
      real*8 zmin,w
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
C               if(.not.mask(i1,i2,i3)) CYCLE
               if(mask(i1,i2,i3)) THEN
C  now search for minima of sms (or weighted sms
               call dcopy(ngrad,si(1,i1,i2,i3),1,sms,1)
               siind(1,i1,i2,i3)=m
               DO k=1,m
                  imin=1
                  zmin=sms(1)
                  DO i=2,ngrad
                     if(sms(i).lt.zmin) THEN
                        zmin=sms(i)
                        imin=i
                     END IF
                  END DO
                  siind(k+1,i1,i2,i3)=imin
                  IF(k.lt.m) THEN
                     DO i=1,ngrad
                        w=dgrad(i,imin)
                        if(w.gt.maxc) THEN
                           sms(i)=1d40
                        ELSE
                           sms(i)=sms(i)/(1.d0-w*w)
                        END IF
                     END DO
                  END IF
               END DO
               ELSE
                  siind(1,i1,i2,i3)=-1
               END IF
            END DO
         END DO
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine getev0(si,ngrad,n1,n2,n3,lev)
      implicit logical (a-z)
      integer n1,n2,n3,ngrad
      real*8 si(ngrad,n1,n2,n3),lev(2,n1,n2,n3)
      integer i1,i2,i3,j
      real*8 z,z1
      z=ngrad
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               lev(1,i1,i2,i3)=0
               z1=si(1,i1,i2,i3)
               DO j=2,ngrad
                  z1=z1+si(j,i1,i2,i3)
               END DO
               lev(2,i1,i2,i3)=-log(z1/z)
            END DO
         END DO
      END DO
      RETURN
      END
 