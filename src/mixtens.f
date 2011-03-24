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
      integer i,j,i3,mode,jpvt(6),rank
      real*8 th,w(ngrad),sw,sth,z1,p0,p1,d1,d2,d3,work(25)
      th = par(1)
      th = max(th,-5.d0)
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
         jpvt(i)=0
      END DO
C  
C    siq will be replaced, need to copy it if C-version of optim is used
C
      call dcopy(ngrad,siq,1,w,1)
      call dgelsy(ngrad,m,1,z,ngrad,w,ngrad,jpvt,1d-6,rank,work,25,
     1            mode)
C  1d-6  is a limit for condition number 
      IF(mode.ne.0) THEN
         call intpr("mode",4,mode,1)
      ELSE
         sw=0.d0
C penalize for extreme th values
         if(th.gt.1.d1) sw=sw+th-1.d1
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
      subroutine mfunpl0h(par,siq,grad,m,lpar,ngrad,z,w,b,
     1                    work1,erg)
C
C   model without isotropic compartment using Larsson-Hansson code
C   same as mfunpl but with unconstrained least squares for weights
C
C   code is restricted to m<=6
C
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad,m),erg,
     1       b(ngrad),work1(ngrad)
      integer i,j,i3,mode,ind(1000)
      real*8 th,w(ngrad),sth,z1,p0,p1,d1,d2,d3,work2(10)
      th = par(1)
      th = max(th,-5.d0)
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
      call dcopy(ngrad,siq,1,b,1)
      call nnls(z,ngrad,ngrad,m,b,w,erg,work2,work1,ind,mode)
      IF(mode.eq.2) erg = 1d10
      call rchkusr()
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine mfunpl1(par,siq,grad,mp1,lpar,ngrad,pen,z,w,erg)
C
C   model with isotropic compartment 
C
C   code is restricted to m<=6
C
      implicit logical (a-z)
      integer m,mp1,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad,mp1),erg,pen
      integer i,j,i3,mode,jpvt(6),rank
      real*8 th,w(ngrad),sw,sth,z1,p0,p1,d1,d2,d3,work(36),eth
      m = mp1-1
      th = par(1)
      th = max(th,-5.d0)
      eth = exp(-th)
      DO j = 1,ngrad
         z(j,1) = eth
      END DO
      DO i = 1,m
C maximal m-1 components
         i3=2*i
         p0=par(i3)
         p1=par(i3+1)
         sth = sin(p0)
         d1 = sth*cos(p1)
         d2 = sth*sin(p1)
         d3 = cos(p0)
         DO j = 1,ngrad
            z1 = d1*grad(1,j)+d2*grad(2,j)+d3*grad(3,j)
            z(j,i+1) = exp(-th*z1*z1)
         END DO
      END DO
C  
C    siq will be replaced, need to copy it if C-version of optim is used
C
      DO i=1,mp1
         jpvt(i)=0
      END DO
      call dcopy(ngrad,siq,1,w,1)
      call dgelsy(ngrad,mp1,1,z,ngrad,w,ngrad,jpvt,1d-6,rank,work,36,
     1            mode)
C  1d-6  is a limit for condition number 
      IF(mode.ne.0) THEN
         call intpr("mode",4,mode,1)
      ELSE
         sw=0.d0
C penalize for extreme th values
         if(th.gt.1.d1) sw=sw+th-1.d1
C penalize for negative weights
         if(th.lt.1.d-2) sw=sw-1.d2*th+1.d0
         DO i=1,mp1
            if(w(i).lt.0.d0) sw=sw-pen*w(i)
         END DO
         DO i=mp1+1,ngrad
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
      subroutine mfunpl0w(par,w,siq,grad,m,lpar,ngrad,z,erg)
C
C   model without isotropic compartment 
C   same as mfunpl but with unconstrained least squares for weights
C
C   code is restricted to m<=6
C
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad,m),w(m)
      integer i,j,i3
      real*8 th,sth,z1,p0,p1,d1,d2,d3,rss,res,erg
      th = par(1)
      th = max(th,-5.d0)
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
      rss =0.d0
      DO j=1,ngrad
         res=siq(j)
         DO i=1,m
            res=res-w(i)*z(j,i)
         END DO
         rss=rss+res*res
      END DO
      erg=rss
      call rchkusr()
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine mfunpl1w(par,w,siq,grad,mp1,lpar,ngrad,z,erg)
C
C   model with isotropic compartment 
C
C   code is restricted to m<=6
C
      implicit logical (a-z)
      integer m,mp1,lpar,ngrad
      real*8 par(lpar),w(mp1),siq(ngrad),grad(3,ngrad),z(ngrad,mp1),
     1       erg
      integer i,j,i3
      real*8 th,sth,z1,p0,p1,d1,d2,d3,eth,rss,res
      m = mp1-1
      th = par(1)
      th = max(th,-5.d0)
      eth = exp(-th)
      DO j = 1,ngrad
         z(j,1) = eth
      END DO
      DO i = 1,m
C maximal m-1 components
         i3=2*i
         p0=par(i3)
         p1=par(i3+1)
         sth = sin(p0)
         d1 = sth*cos(p1)
         d2 = sth*sin(p1)
         d3 = cos(p0)
         DO j = 1,ngrad
            z1 = d1*grad(1,j)+d2*grad(2,j)+d3*grad(3,j)
            z(j,i+1) = exp(-th*z1*z1)
         END DO
      END DO
      rss =0.d0
      DO j=1,ngrad
         res=siq(j)
         DO i=1,mp1
            res=res-w(i)*z(j,i)
         END DO
         rss=rss+res*res
      END DO
      erg=rss
      call rchkusr()
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine mfpl0gn(par,siq,grad,m,lpar,ngrad,pen,eps,z,w,
     1                   para,parb,dfdpar)
C
C   model without isotropic compartment 
C   same as mfunpl but with unconstrained least squares for weights
C
C   code is restricted to m<=6
C
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad,m),pen,
     1       dfdpar(lpar),para(lpar),parb(lpar),eps
      real*8 erga,ergb,deltai
      integer i
      deltai=0.5d0/eps
      DO i=1,lpar
         call dcopy(lpar,par,1,para,1)
         call dcopy(lpar,par,1,parb,1)
         para(i)=para(i)-eps
         parb(i)=parb(i)+eps
         call mfunpl0(para,siq,grad,m,lpar,ngrad,pen,z,w,erga)
         call mfunpl0(parb,siq,grad,m,lpar,ngrad,pen,z,w,ergb)
         if(max(ergb,erga).lt.1d10) THEN
            dfdpar(i)=(ergb-erga)*deltai
         ELSE
            dfdpar(i)=0.d0
         ENDIF
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine mfpl0hgn(par,siq,grad,m,lpar,ngrad,eps,z,w,b,
     1                    work1,para,parb,dfdpar)
C
C   model without isotropic compartment 
C   same as mfunpl but with unconstrained least squares for weights
C
C   code is restricted to m<=6
C
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad,m),b(lpar),
     1       dfdpar(lpar),para(lpar),parb(lpar),eps,work1(ngrad)
      real*8 erga,ergb,deltai
      integer i
      deltai=0.5d0/eps
      DO i=1,lpar
         call dcopy(lpar,par,1,para,1)
         call dcopy(lpar,par,1,parb,1)
         para(i)=para(i)-eps
         parb(i)=parb(i)+eps
         call mfunpl0h(para,siq,grad,m,lpar,ngrad,z,w,b,work1,erga)
         call mfunpl0h(parb,siq,grad,m,lpar,ngrad,z,w,b,work1,ergb)
         if(max(ergb,erga).lt.1d10) THEN
            dfdpar(i)=(ergb-erga)*deltai
         ELSE
            dfdpar(i)=0.d0
         ENDIF
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine mfunpl0g(par,s,g,m,lpar,n,d,z,v,w,dkgj,dkgj2,
     1                    ddkdphig,ddkdetag,dddphi,
     2                    dddeta,dvdth,dvdphi,dvdeta,dzdpars,
     3                    dwdpars,dwdpars2,zs,work1,
     3                    work2,scopy,pen,dfdpar)
      implicit logical (a-z)
      integer m,n,lpar
      real*8 par(lpar),s(n),g(3,n),d(3,m),z(n,m),v(m,m),dkgj(n,m),
     1       w(n),dkgj2(n,m),ddkdphig(n,m),ddkdetag(n,m),dddphi(3,m),
     2       dddeta(3,m),dvdth(m,m),dvdphi(m,m,m),dvdeta(m,m,m),
     3       dzdpars(n,m,3),dwdpars(m,lpar),dwdpars2(m,lpar),
     4       zs(n,m),dfdpar(lpar),pen,rcond,ferr(11),berr(11)
      integer i,j,k,l,i3,ind(5),mode,jpvt(5),rank
      real*8 th,sw,sphi,cphi,seta,ceta,z1,z2,p0,p1,
     1       work(250),work1(n,m),work2(n,m),scopy(n),m2th,af(25)
C  length of work needs to be larger or equal max(m*m,ngrad)
      real*8 ddot
      external ddot
      th = par(1)
      th = max(-5.d0,th)
      m2th = -2.d0*th
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
      DO i=1,m
         jpvt(i)=0
      END DO
C  
C   we now have d, dkgj,dddphi, dddeta, z, ddkdphig, ddkdetag
C   next w
C
      call dcopy(n,s,1,w,1)
      call dcopy(n*m,z,1,work1,1)
      call dgelsy(n,m,1,work1,n,w,n,jpvt,1d-6,rank,work,25,mode)
C  1d-6  is a limit for condition number 
      IF(mode.gt.1) THEN
         call intpr("mode",4,mode,1)
      END IF
C
C   thats weights in w now V, dVdth, dVdphi, dVdeta, dzdpars
C
C   use work1, work2 and scopy for intermediate results
C   componentwise products
      call dcprod0(dkgj,ddkdphig,n*m,work1)
      call dcprod0(dkgj,ddkdetag,n*m,work2)
C initialize unneeded elements
      call zerofill(dzdpars,m*3*n)
C compute componentswise -dkgj2*z
      call dcprod(dkgj2,z,-1.d0,n*m,dzdpars(1,1,1))
C compute componentswise m2th*work1*z
      call dcprod(work1,z,m2th,n*m,dzdpars(1,1,2))
C compute componentswise m2th*work2*z
      call dcprod(work2,z,m2th,n*m,dzdpars(1,1,3))
      call zerofill(dVdphi,m*m*m)
      call zerofill(dVdeta,m*m*m)
      DO k=1,m
         v(k,k)=ddot(n,z(1,k),1,z(1,k),1)
         z1 = ddot(n,dzdpars(1,k,1),1,z(1,k),1)
         dVdth(k,k) = 2.d0*z1
         DO l=k+1,m
            z2=ddot(n,z(1,k),1,z(1,l),1)
            v(l,k)=z2
            v(k,l)=z2
            z2=ddot(n,dzdpars(1,l,1),1,z(1,k),1)+
     1         ddot(n,dzdpars(1,k,1),1,z(1,l),1)
            dVdth(k,l) = z2
            dVdth(l,k) = z2
         END DO
         DO i=1,m
            z2 = ddot(n,dzdpars(1,k,2),1,z(1,i),1)
            dVdphi(i,k,k) = dVdphi(i,k,k) + z2
            dVdphi(k,i,k) = dVdphi(k,i,k) + z2
            z2 = ddot(n,dzdpars(1,k,3),1,z(1,i),1)
            dVdeta(i,k,k) = dVdeta(i,k,k) + z2
            dVdeta(k,i,k) = dVdeta(k,i,k) + z2
         END DO
      END DO
C
C   thats V, dVdth, dVdphi, dVdeta now fill dwdpars)
C
      DO k=1,m
         dwdpars(k,1) = ddot(n,dzdpars(1,k,1),1,s,1) -
     1                         ddot(m,dVdth(1,k),1,w,1)
         DO l=1,m
            dwdpars(k,1+l) = -ddot(m,dVdphi(1,k,l),1,w,1)    
            dwdpars(k,1+m+l) = -ddot(m,dVdeta(1,k,l),1,w,1)    
         END DO
         dwdpars(k,1+k) = dwdpars(k,1+k) + 
     1                    ddot(n,dzdpars(1,k,2),1,s,1)
         dwdpars(k,1+m+k) = dwdpars(k,1+m+k) + 
     1                    ddot(n,dzdpars(1,k,3),1,s,1)
      END DO
C
C   thats  dzdpars  now compute  dw/dpar in dwdpars
C
      call dcopy(m*lpar,dwdpars,1,dwdpars2,1)
      call dsysvx("N","U",m,lpar,v,m,af,m,ind,dwdpars2,m,dwdpars,
     1            m,rcond,ferr,berr,work,25,jpvt,mode)
      IF(mode.ne.0.or.rcond.lt.1d-8) THEN
C   solving the linear system fails due to renk deficiency of v
C compute numerical gradients instead
        call mfpl0gn(par,s,g,m,lpar,n,pen,1d-6,z,w,work1,work2,dfdpar)
         RETURN
      END IF
C
C   now we have dw/dpar in dzdpars next residuals in scopy
C
      DO j=1,n
         z1 = s(j)
         z2 = 0.d0
         DO k=1,m
            z1 = z1 - w(k)*z(j,k)
            z2=z2+w(k)*dzdpars(j,k,1)+dwdpars(k,1)*z(j,k)
         END DO
         scopy(j) = z1
         work(j)=z2
      END DO
C
C   now we have residuals in scopy compute gradient of f
C   use work for intermediate results
      dfdpar(1)=-2.d0*ddot(n,work,1,scopy,1)
      if(th.gt.1.d1) dfdpar(1)=dfdpar(1)+1.d0
      if(th.lt.1d-2) dfdpar(1)=dfdpar(1)-1.d2
      DO k=1,m
         if(w(k).lt.0.d0) dfdpar(1)=dfdpar(1)-pen*dwdpars(k,1)
      END DO
C    thats derivative with respect to theta
C 
C    dzdpars contains dw/dpars 
C
      DO i = 1,m
         DO j = 1,n
            z1=w(i)*dzdpars(j,i,2)
            DO k=1,m
               z1=z1+dwdpars(k,1+i)*z(j,k)
            END DO
            work(j)=z1
         END DO
         dfdpar(2*i)=-2.d0*ddot(n,work,1,scopy,1)
         DO k=1,m
            if(w(k).lt.0) dfdpar(2*i)=dfdpar(2*i)-
     1                                    pen*dwdpars(k,1+i)
         END DO
C    thats derivative with respect to phi
         DO j = 1,n
            z1=w(i)*dzdpars(j,i,3)
            DO k=1,m
               z1=z1+dwdpars(k,1+m+i)*z(j,k)
            END DO
            work(j)=z1
         END DO
         dfdpar(2*i+1)=-2.d0*ddot(n,work,1,scopy,1)
         DO k=1,m
            if(w(k).lt.0) dfdpar(1+2*i)=dfdpar(1+2*i)-
     1                                pen*dwdpars(k,1+m+i)
         END DO
      END DO
C
C   thats derivative with respect to eta
C
      call rchkusr()
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine mfpl1gn(par,siq,grad,m,lpar,ngrad,pen,eps,z,w,
     1                   para,parb,dfdpar)
C
C   model with isotropic compartment ( m = #comp +1 )
C   same as mfunpl but with unconstrained least squares for weights
C
C   code is restricted to m<=6
C
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad,m),pen,
     1       dfdpar(lpar),para(lpar),parb(lpar),eps
      real*8 erga,ergb,deltai
      integer i
      deltai=0.5d0/eps
      DO i=1,lpar
         call dcopy(lpar,par,1,para,1)
         call dcopy(lpar,par,1,parb,1)
         para(i)=para(i)-eps
         parb(i)=parb(i)+eps
         call mfunpl1(para,siq,grad,m,lpar,ngrad,pen,z,w,erga)
         call mfunpl1(parb,siq,grad,m,lpar,ngrad,pen,z,w,ergb)
         if(max(ergb,erga).lt.1d10) THEN
            dfdpar(i)=(ergb-erga)*deltai
         ELSE
            dfdpar(i)=0.d0
         ENDIF
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine mfunpl1g(par,s,g,m,mp1,lpar,n,d,z,v,w,dkgj,dkgj2,
     1                    ddkdphig,ddkdetag,dddphi,dddeta,dvdth,
     2                    dvdphi,dvdeta,dzdpars,dwdpars,dwdpars2,
     3                    zs,work1,work2,scopy,pen,dfdpar)
      implicit logical (a-z)
      integer m,mp1,n,lpar
      real*8 par(lpar),s(n),g(3,n),d(3,m),z(n,mp1),v(mp1,mp1),
     1       dkgj(n,m),w(n),dkgj2(n,m),ddkdphig(n,m),ddkdetag(n,m),
     2       dddphi(3,m),dddeta(3,m),dvdth(mp1,mp1),dvdphi(mp1,mp1,m),
     3       dvdeta(mp1,mp1,m),dzdpars(n,mp1,3),dwdpars(mp1,lpar),
     4       dwdpars2(mp1,lpar),zs(n,mp1),dfdpar(lpar),pen
      integer i,j,k,l,i3,ind(6),mode,ip1,kp1,lp1,jpvt(6),rank
      real*8 th,sw,sphi,cphi,seta,ceta,z1,z2,p0,p1,work(250),rcond,
     1       work1(n,mp1),work2(n,mp1),scopy(n),zji,m2th,eth,af(36),
     2       ferr(12),berr(12)
      real*8 ddot
      external ddot
      th = par(1)
      th = max(-5.d0,th)
      eth=exp(-th)
      m2th = -2.d0*th
      sw = 0
C
C  get d, dkgj, dkgj2, dddphi, dddeta, z, ddkdphig, ddkdetag
C
      DO i = 1,m
         ip1=i+1
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
            zji=exp(-th*z2)
            z(j,ip1) = zji
            zs(j,ip1) = zji*s(j)
            ddkdphig(j,i) = ddot(3,dddphi(1,i),1,g(1,j),1)
            ddkdetag(j,i) = ddot(3,dddeta(1,i),1,g(1,j),1)
         END DO
      END DO
C  isotrope part in z
      DO i=1,n
         z(i,1) = eth
         zs(i,1) = eth*s(i)
      END DO
C  
C   we now have d, dkgj,dddphi, dddeta, z, ddkdphig, ddkdetag
C   next w
C
      call dcopy(n,s,1,w,1)
      call dcopy(n*mp1,z,1,work1,1)
      call dgelsy(n,mp1,1,work1,n,w,n,jpvt,1d-6,rank,work,36,mode)
      IF(mode.gt.1) THEN
         call intpr("mode",4,mode,1)
      END IF
C
C   thats weights in w now V, dVdth, dVdphi, dVdeta, dzdpars
C
C   use work1, work2 and scopy for intermediate results
C   componentwise products
      call dcprod0(dkgj,ddkdphig,n*m,work1)
      call dcprod0(dkgj,ddkdetag,n*m,work2)
C initialize unneeded elements
      call zerofill(dzdpars,mp1*3*n)
C  derivatives of z(,1) are zero (isotrop part) dzdpars(.,1,.)=0
      DO j=1,n
         dzdpars(j,1,1)=-eth
      END DO
C compute componentswise -dkgj2*z(,+1)
      call dcprod(dkgj2,z(1,2),-1.d0,n*m,dzdpars(1,2,1))
C compute componentswise m2th*work1*z(,+1)
      call dcprod(work1,z(1,2),m2th,n*m,dzdpars(1,2,2))
C compute componentswise m2th*work2*z(,+1)
      call dcprod(work2,z(1,2),m2th,n*m,dzdpars(1,2,3))
C      DO k=1,m
C         kp1=k+1
C         DO j=1,n
C            dzdpars(j,kp1,1)=-dkgj2(j,k)*z(j,kp1)
C            dzdpars(j,kp1,kp1)=m2th*work1(j,k)*z(j,kp1)
C            dzdpars(j,kp1,m+kp1)=m2th*work2(j,k)*z(j,kp1)
C         END DO
C      END DO
C  derivatives of v with respect to theta
      DO k=1,mp1
         v(k,k)=ddot(n,z(1,k),1,z(1,k),1)
         z1 = ddot(n,dzdpars(1,k,1),1,z(1,k),1)
         dVdth(k,k) = 2.d0*z1
         DO l=k+1,mp1
            z2=ddot(n,z(1,k),1,z(1,l),1)
            v(l,k)=z2
            v(k,l)=z2
            z2=ddot(n,dzdpars(1,l,1),1,z(1,k),1)+
     1         ddot(n,dzdpars(1,k,1),1,z(1,l),1)
            dVdth(k,l) = z2
            dVdth(l,k) = z2
         END DO
      END DO
C  derivatives of v with respect to theta 
      call zerofill(dVdphi,m*mp1*mp1)
      call zerofill(dVdeta,m*mp1*mp1)
      DO k=1,m
         kp1=k+1
         DO i=1,m
            ip1=i+1
            z2 = ddot(n,dzdpars(1,kp1,2),1,z(1,ip1),1)
            dVdphi(ip1,kp1,k) = dVdphi(ip1,kp1,k)+ z2
            dVdphi(kp1,ip1,k) = dVdphi(kp1,ip1,k)+ z2
            z2 = ddot(n,dzdpars(1,kp1,3),1,z(1,ip1),1)
            dVdeta(ip1,kp1,k) = dVdeta(ip1,kp1,k)+ z2
            dVdeta(kp1,ip1,k) = dVdeta(kp1,ip1,k)+ z2
         END DO
      END DO
C
C   thats V, dVdth, dVdphi, dVdeta now fill dwdpars)
C
      DO k=1,mp1
         dwdpars(k,1) = ddot(n,dzdpars(1,k,1),1,s,1) -
     1                         ddot(mp1,dVdth(1,k),1,w,1)
         kp1=k+1
         DO l=1,m
            lp1=l+1
            dwdpars(k,lp1) = - ddot(mp1,dVdphi(1,k,l),1,w,1)
            dwdpars(k,m+lp1) = - ddot(mp1,dVdeta(1,k,l),1,w,1)
         END DO
         IF(k.ne.1) THEN
         dwdpars(k,k) = dwdpars(k,k) +
     1                    ddot(n,dzdpars(1,k,2),1,s,1)
         dwdpars(k,m+k) = dwdpars(k,m+k) +
     1                    ddot(n,dzdpars(1,k,3),1,s,1)
         END IF
      END DO
C
C   thats  dzdpars  now compute  dw/dpar in dwdpars
C
      call dcopy(mp1*lpar,dwdpars,1,dwdpars2,1)
      call dsysvx("N","U",mp1,lpar,v,mp1,af,mp1,ind,dwdpars2,mp1,
     1            dwdpars,mp1,rcond,ferr,berr,work,36,jpvt,mode)
      IF(mode.ne.0.or.rcond.lt.1d-8) THEN
C   solving the linear system fails due to renk deficiency of v
C compute numerical gradients instead
      call mfpl1gn(par,s,g,mp1,lpar,n,pen,1d-6,z,w,work1,work2,dfdpar)
         RETURN
      END IF
C
C   now we have dw/dpar in dzdpars next residuals in scopy
C
      DO j=1,n
         z1 = s(j)
         z2=0.d0
         DO k=1,mp1
            z1 = z1 - w(k)*z(j,k)
            z2=z2+w(k)*dzdpars(j,k,1)+dwdpars(k,1)*z(j,k)
         END DO
         scopy(j) = z1
         work(j)=z2
      END DO
C
C   now we have residuals in scopy compute gradient of f
C   use work for intermediate results
C   z(j,1)=0 und dzdpars(j,1,1) = 0
      dfdpar(1)=-2.d0*ddot(n,work,1,scopy,1)
      if(th.gt.1.d1) dfdpar(1)=dfdpar(1)+1.d0
      if(th.lt.1d-2) dfdpar(1)=dfdpar(1)-1.d2
      DO k=1,mp1
         if(w(k).lt.0.d0) dfdpar(1)=dfdpar(1)-pen*dwdpars(k,1)
      END DO
C    thats derivative with respect to theta
C 
C    dzdpars contains dw/dpars  dzdpars(j,1,.)=0 
C
      DO i = 1,m
         ip1=i+1
         DO j = 1,n
            z1=w(ip1)*dzdpars(j,ip1,2)
            DO k=1,mp1
               z1=z1+dwdpars(k,ip1)*z(j,k)
            END DO
            work(j)=z1
         END DO
         dfdpar(2*i)=-2.d0*ddot(n,work,1,scopy,1)
         DO k=1,mp1
            if(w(k).lt.0) dfdpar(2*i)=dfdpar(2*i)-
     1                                    pen*dwdpars(k,ip1)
         END DO
C    thats derivative with respect to phi
         DO j = 1,n
            z1=w(ip1)*dzdpars(j,ip1,3)
            DO k=1,mp1
               z1=z1+dwdpars(k,m+ip1)*z(j,k)
            END DO
            work(j)=z1
         END DO
         dfdpar(2*i+1)=-2.d0*ddot(n,work,1,scopy,1)
         DO k=1,mp1
            if(w(k).lt.0) dfdpar(1+2*i)=dfdpar(1+2*i)-
     1                                pen*dwdpars(k,m+ip1)
         END DO
      END DO
C
C   thats derivative with respect to eta
C
      z1=w(1)
      DO i=2,mp1
         z1=min(z1,w(i))
      END DO
      call rchkusr()
      RETURN
      END 
C
C __________________________________________________________________
C
      subroutine getsii30(si,vsi,ngrad,nvox,m,dgrad,nv,th,
     1     nth,indth,egrad,isample,ntry,sms,z,siind,mval,ns,mask)
C
C  compute diagnostics for initial estimates in siind
C  siind(1,i1,i2,i3) will contain the model order 
C  
C  si     - array of si-values
C  m      - model order
C  maxc   - maximum of cos(angle between directions)
C  th     - theta1
C  egrad - exp(-theta1*dgrad^2) 
C  isample - guesses for gradient directions
C  ntry   - number of guesses
C  sms    - copies of si
C  z      - array for design matrix corresponding to guesses
C  siind  - array of indices (output)
C  ns     - m+1
C  mask   - mask
C  mval   - aktual best risk
C
C  restricted to ngrad<=1000 and m <=10
C
      implicit logical (a-z)
      integer nvox,ngrad,ns,siind(ns,nvox),m,ntry,nth,nv,
     1       isample(m,ntry),indth(nvox)
      real*8 si(ngrad,nvox),sms(ngrad),dgrad(ngrad,nv),
     1       th(nvox),egrad(ngrad,nv),z(ngrad,ns),mval(nvox),
     2       vsi(nvox)
      logical mask(nvox)
      integer i,k,ibest,mode,ind(10),l,j,ii,iw,wind(5),nwi(5)
      real*8 w(1000),krit,work1(1000),work2(10),erg,thj,msi,m2si,
     1       z1,dng
      dng=ngrad
      iw=m
      DO i=1,m
         wind(i)=i
         nwi(i)=i
      END DO
      ibest=1
      DO i=1,nvox
         msi=0.d0
         m2si=0.d0
         z1=vsi(i)
         mval(i)=sqrt(dng*z1)
         if(.not.mask(i)) THEN
            siind(1,i)=-1
            mval(i)=0
         END IF
      END DO
      call rchkusr()
      DO j=1,nth
         thj=th(j)
         DO k=1,ngrad
            DO l=1,nv
               z1 = dgrad(k,l)
               egrad(k,l)=dexp(-thj*z1*z1)
            END DO
         END DO
         DO i=1,nvox
            if(mask(i)) THEN
               IF(j.ne.indth(i)) CYCLE
C  now search for minima of sms (or weighted sms
               ibest=0
               krit=mval(i)
               DO k=1,ntry
                  call dcopy(ngrad,si(1,i),1,sms,1)
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
                        iw=0
                        DO ii=1,m
                           if(w(ii).gt.1.d-12) THEN
                              iw=iw+1
                              wind(iw)=ii
                           ELSE
                              nwi(ii-iw)=ii
C   nonactive directions
                           END IF 
                        END DO
                     END IF  
                  END IF
               END DO
               if(ibest.gt.0) THEN
                  siind(1,i)=iw
                  siind(2,i)=j
                  IF (iw.gt.1) THEN
                     DO l=1,iw
                        siind(l+2,i)=isample(wind(l),ibest)
                     END DO
                  END IF
                  IF (iw.lt.m) THEN
                     DO l=1,m-iw
                        siind(m-l+3,i)=isample(nwi(l),ibest)
                     END DO
                  END IF
                  mval(i)=krit
               END IF
            END IF
            call rchkusr()
         END DO
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine getsii31(si,vsi,ngrad,nvox,m,dgrad,nv,iandir,th,
     1     nth,indth,egrad,isample,ntry,sms,z,siind,mval,ns,mask,
     2     dgradv,maxc)
C
C  compute diagnostics for initial estimates in siind
C  siind(1,i1,i2,i3) will contain the model order 
C  
C  si     - array of si-values
C  m      - model order
C  maxc   - maximum of cos(angle between directions)
C  th     - theta1
C  egrad - exp(-theta1*dgrad^2) 
C  isample - guesses for gradient directions
C  ntry   - number of guesses
C  sms    - copies of si
C  z      - array for design matrix corresponding to guesses
C  siind  - array of indices (output)
C  ns     - m+1
C  mask   - mask
C  mval   - aktual best risk
C
C  restricted to ngrad<=1000 and m <=10
C
      implicit logical (a-z)
      integer nvox,ngrad,ns,siind(ns,nvox),m,ntry,nth,nv,
     1       isample(1),indth(nvox),iandir(nvox)
      real*8 si(ngrad,nvox),sms(ngrad),dgrad(ngrad,nv),
     1       th(nvox),egrad(ngrad,nv),z(ngrad,ns),mval(nvox),
     2       vsi(nvox),dgradv(nv,nv),maxc
      logical mask(nvox),skip
      integer i,k,mode,ind(10),l,j,ii,iw,wind(5),nwi(5),mis,
     1        is(5),isbest(5),ntry0,iii
      real*8 w(1000),krit,work1(1000),work2(10),erg,thj,msi,m2si,
     1       z1,dng
      dng=ngrad
      mis=m-1
      ntry0=ntry
      if(m.eq.1) ntry0=1
      iw=m
      DO i=1,m
         wind(i)=i
         nwi(i)=i
         isbest(i)=i
      END DO
      DO i=1,nvox
         msi=0.d0
         m2si=0.d0
         z1=vsi(i)
         mval(i)=sqrt(dng*z1)
         if(.not.mask(i)) THEN
            siind(1,i)=-1
            mval(i)=0
         END IF
      END DO
      call rchkusr()
      DO j=1,nth
         thj=th(j)
         DO k=1,ngrad
            DO l=1,nv
               z1 = dgrad(k,l)
               egrad(k,l)=dexp(-thj*z1*z1)
            END DO
         END DO
         DO i=1,nvox
            if(mask(i)) THEN
               IF(j.ne.indth(i)) CYCLE
C  now search for minima of sms (or weighted sms)
               krit=mval(i)
               DO k=1,ntry0
                  IF(m.gt.1) THEN
                     skip=.FALSE.
                     DO l=1,m-1
                        iii=mis*(k-1)+l
                   if(dgradv(isample(mis*(k-1)+l),iandir(i)).gt.maxc)
     1                    skip=.TRUE.
                     END DO
                     IF(skip) CYCLE
                  END IF
                  call dcopy(ngrad,si(1,i),1,sms,1)
                  IF(m.gt.1) THEN
                     DO l=1,m-1
                        is(l)=isample(mis*(k-1)+l)
                        call dcopy(ngrad,egrad(1,is(l)),1,z(1,l),1)
                     END DO
                  END IF
                  is(m)=iandir(i)
                  call dcopy(ngrad,egrad(1,is(m)),1,z(1,m),1)
        call nnls(z,ngrad,ngrad,m,sms,w,erg,work2,work1,ind,mode)
                  IF(mode.gt.1) THEN
                     call intpr("mode",4,mode,1)
                     call intpr("isample",7,is,m)
                  ELSE 
                     IF(erg.lt.krit) THEN
                        krit=erg
                        iw=0
                        DO ii=1,m
                           isbest(ii)=is(ii)
                           if(w(ii).gt.1.d-12) THEN
                              iw=iw+1
                              wind(iw)=ii
                           ELSE
                              nwi(ii-iw)=ii
C   nonactive directions
                           END IF 
                        END DO
                     END IF  
                  END IF
               END DO
               siind(1,i)=iw
               siind(2,i)=j
               IF (iw.gt.1) THEN
                  DO l=1,iw
                     siind(l+2,i)=isbest(wind(l))
                  END DO
               END IF
               IF (iw.lt.m) THEN
                  DO l=1,m-iw
                     siind(m-l+3,i)=isbest(nwi(l))
                  END DO
               END IF
               mval(i)=krit
            END IF
            call rchkusr()
         END DO
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine getsi30i(si,vsi,ngrad,nvox,m,dgrad,nv,th,
     1     nth,indth,egrad,isample,ntry,sms,z,siind,mval,ns,mask)
C
C  compute diagnostics for initial estimates in siind
C  siind(1,i1,i2,i3) will contain the model order 
C  
C  si     - array of si-values
C  m      - model order
C  maxc   - maximum of cos(angle between directions)
C  th     - theta1
C  egrad - exp(-theta1*dgrad^2) 
C  isample - guesses for gradient directions
C  ntry   - number of guesses
C  sms    - copies of si
C  z      - array for design matrix corresponding to guesses
C  siind  - array of indices (output)
C  ns     - m+1
C  mask   - mask
C  mval   - aktual best risk
C
C  restricted to ngrad<=1000 and m <=10
C
      implicit logical (a-z)
      integer nvox,ngrad,ns,siind(ns,nvox),m,ntry,nth,nv,
     1       isample(m,ntry),indth(nvox)
      real*8 si(ngrad,nvox),sms(ngrad),dgrad(ngrad,nv),
     1       th(nvox),egrad(ngrad,nv),z(ngrad,1),mval(nvox),
     2       vsi(nvox)
      logical mask(nvox)
      integer i,k,ibest,mode,ind(12),l,j,ii,iw,wind(5),nwi(5),mp1
      real*8 w(1000),krit,work1(1000),work2(10),erg,thj,msi,m2si,
     1       z1,dng,ethj
      dng=ngrad
      iw=m
      mp1=m+1
      DO i=1,m
         wind(i)=i
         nwi(i)=i
      END DO
      ibest=1
      DO i=1,nvox
         msi=0.d0
         m2si=0.d0
         z1=vsi(i)
         mval(i)=sqrt(dng*z1)
         if(.not.mask(i)) THEN
            siind(1,i)=-1
            mval(i)=0
         END IF
      END DO
      call rchkusr()
      DO j=1,nth
         thj=th(j)
         ethj=dexp(-thj)
         DO k=1,ngrad
            DO l=1,nv
               z1 = dgrad(k,l)
               egrad(k,l)=dexp(-thj*z1*z1)
            END DO
         END DO
         DO i=1,nvox
            if(mask(i)) THEN
               IF(j.ne.indth(i)) CYCLE
C  now search for minima of sms (or weighted sms
               ibest=0
               krit=mval(i)
               DO k=1,ntry
                  DO l=1,ngrad
                     z(l,1)=ethj
                  END DO
                  call dcopy(ngrad,si(1,i),1,sms,1)
                  DO l=1,m
                  call dcopy(ngrad,egrad(1,isample(l,k)),1,z(1,l+1),1)
                  END DO
           call nnls(z,ngrad,ngrad,mp1,sms,w,erg,work2,work1,ind,mode)
                  IF(mode.gt.1) THEN
                     call intpr("mode",4,mode,1)
                     call intpr("isample",7,isample(1,k),m)
                  ELSE 
                     IF(erg.lt.krit) THEN
                        krit=erg
                        ibest=k
                        iw=0
                        DO ii=1,m
                           if(w(ii).gt.1.d-12) THEN
                              iw=iw+1
                              wind(iw)=ii
                           ELSE
                              nwi(ii-iw)=ii
C   nonactive directions
                           END IF 
                        END DO
                     END IF  
                  END IF
               END DO
               if(ibest.gt.0) THEN
                  siind(1,i)=iw
                  siind(2,i)=j
                  IF (iw.gt.1) THEN
                     DO l=1,iw
                        siind(l+2,i)=isample(wind(l),ibest)
                     END DO
                  END IF
                  IF (iw.lt.m) THEN
                     DO l=1,m-iw
                        siind(m-l+3,i)=isample(nwi(l),ibest)
                     END DO
                  END IF
                  mval(i)=krit
               END IF
            END IF
            call rchkusr()
         END DO
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine getsi31i(si,vsi,ngrad,nvox,m,dgrad,nv,iandir,th,
     1     nth,indth,egrad,isample,ntry,sms,z,siind,mval,ns,mask,
     2     dgradv,maxc)
C
C  compute diagnostics for initial estimates in siind
C  siind(1,i1,i2,i3) will contain the model order 
C  
C  si     - array of si-values
C  m      - model order
C  maxc   - maximum of cos(angle between directions)
C  th     - theta1
C  egrad - exp(-theta1*dgrad^2) 
C  isample - guesses for gradient directions
C  ntry   - number of guesses
C  sms    - copies of si
C  z      - array for design matrix corresponding to guesses
C  siind  - array of indices (output)
C  ns     - m+1
C  mask   - mask
C  mval   - aktual best risk
C
C  restricted to ngrad<=1000 and m <=10
C
      implicit logical (a-z)
      integer nvox,ngrad,ns,siind(ns,nvox),m,ntry,nth,nv,
     1       isample(1),indth(nvox),iandir(nvox)
      real*8 si(ngrad,nvox),sms(ngrad),dgrad(ngrad,nv),
     1       th(nvox),egrad(ngrad,nv),z(ngrad,1),mval(nvox),
     2       vsi(nvox),dgradv(nv,nv),maxc
      logical mask(nvox),skip
      integer i,k,mode,ind(12),l,j,ii,iw,wind(5),nwi(5),mis,
     1        is(5),isbest(5),ntry0,iii,mp1
      real*8 w(1000),krit,work1(1000),work2(10),erg,thj,msi,m2si,
     1       z1,dng,ethj
      dng=ngrad
      mis=m-1
      ntry0=ntry
      if(m.eq.1) ntry0=1
      iw=m
      mp1=m+1
      DO i=1,m
         wind(i)=i
         nwi(i)=i
         isbest(i)=i
      END DO
      DO i=1,nvox
         msi=0.d0
         m2si=0.d0
         z1=vsi(i)
         mval(i)=sqrt(dng*z1)
         if(.not.mask(i)) THEN
            siind(1,i)=-1
            mval(i)=0
         END IF
      END DO
      call rchkusr()
      DO j=1,nth
         thj=th(j)
         ethj=dexp(-thj)
         DO k=1,ngrad
            DO l=1,nv
               z1 = dgrad(k,l)
               egrad(k,l)=dexp(-thj*z1*z1)
            END DO
         END DO
         DO i=1,nvox
            if(mask(i)) THEN
               IF(j.ne.indth(i)) CYCLE
C  now search for minima of sms (or weighted sms)
               krit=mval(i)
               DO k=1,ntry0
                  IF(m.gt.1) THEN
                     skip=.FALSE.
                     DO l=1,m-1
                        iii=mis*(k-1)+l
                   if(dgradv(isample(mis*(k-1)+l),iandir(i)).gt.maxc)
     1                    skip=.TRUE.
                     END DO
                     IF(skip) CYCLE
                  END IF
                  call dcopy(ngrad,si(1,i),1,sms,1)
                  DO l=1,ngrad
                     z(l,1)=ethj
                  END DO
                  IF(m.gt.1) THEN
                     DO l=1,m-1
                        is(l)=isample(mis*(k-1)+l)
                        call dcopy(ngrad,egrad(1,is(l)),1,z(1,l+1),1)
                     END DO
                  END IF
                  is(m)=iandir(i)
                  call dcopy(ngrad,egrad(1,is(m)),1,z(1,mp1),1)
        call nnls(z,ngrad,ngrad,m,sms,w,erg,work2,work1,ind,mode)
                  IF(mode.gt.1) THEN
                     call intpr("mode",4,mode,1)
                     call intpr("isample",7,is,m)
                  ELSE 
                     IF(erg.lt.krit) THEN
                        krit=erg
                        iw=0
                        DO ii=1,m
                           isbest(ii)=is(ii)
                           if(w(ii).gt.1.d-12) THEN
                              iw=iw+1
                              wind(iw)=ii
                           ELSE
                              nwi(ii-iw)=ii
C   nonactive directions
                           END IF 
                        END DO
                     END IF  
                  END IF
               END DO
               siind(1,i)=iw
               siind(2,i)=j
               IF (iw.gt.1) THEN
                  DO l=1,iw
                     siind(l+2,i)=isbest(wind(l))
                  END DO
               END IF
               IF (iw.lt.m) THEN
                  DO l=1,m-iw
                     siind(m-l+3,i)=isbest(nwi(l))
                  END DO
               END IF
               mval(i)=krit
            END IF
            call rchkusr()
         END DO
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine sweeps0(si,s0,n1,n2,n3,ng0,ng1,level,siq,ms0,vsi,
     1                   mask)
C
C   calculate mean s0 value
C   generate mask
C   sweep s0 from si to generate  siq
C   calculate variance of siq
C
      integer n1,n2,n3,ng0,ng1,si(n1,n2,n3,ng1),s0(n1,n2,n3,ng0),
     1        level
      real*8 siq(n1,n2,n3,ng1),ms0(n1,n2,n3),vsi(n1,n2,n3)
      logical mask(n1,n2,n3),maskk
      integer i1,i2,i3,k
      real*8 s,z,z2,thresh,cv,s0mean
      thresh = level*ng0
      cv=ng1*(ng1-1)
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               z=0.d0
               DO k=1,ng0
                  z=z+s0(i1,i2,i3,k)
               END DO
               s0mean = z/ng0
               ms0(i1,i2,i3) = s0mean
               maskk = z.ge.thresh
               IF(maskk) THEN
                  z=0.d0
                  z2=0.d0
                  DO k=1,ng1
                     s=si(i1,i2,i3,k)/s0mean
                     if(s.gt.0.99d0) s=0.99d0
                     z=z+s
                     z2=z2+s*s
                     siq(i1,i2,i3,k)=s
                  END DO
                  vsi(i1,i2,i3)=(ng1*z2-z)/cv
                  if(vsi(i1,i2,i3).lt.1d-8) THEN
                     maskk = .FALSE.
                     vsi(i1,i2,i3)=0.d0
                  END IF
               ELSE
                  vsi(i1,i2,i3)=0.d0
                  DO k=1,ng1
                     siq(i1,i2,i3,k)=1.d0
                  END DO
               END IF
               mask(i1,i2,i3) = maskk
            END DO
         END DO
      call rchkusr()
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine iandir(vico,nvico,andir,nvox,landir,iandi)
      implicit logical(a-z)
      integer nvico,nvox,iandi(2,nvox)
      real*8 vico(3,nvico),andir(3,2,nvox)
      logical landir(nvox)
      integer i,j,jmax
      real*8 z,zmax,scprod3
      external scprod3
      DO i=1,nvox
         if(landir(i)) THEN
            zmax = scprod3(vico(1,1),andir(1,1,i))
            jmax = 1
            DO j=2,nvico
               z = scprod3(vico(1,j),andir(1,1,i))
               if(z.gt.zmax) THEN
                  zmax=z
                  jmax=j
               END IF
            END DO
            iandi(1,i)=jmax
         END IF
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      real*8 function scprod3(a,b)
      implicit logical (a-z) 
      real*8 a(3),b(3)
      scprod3=abs(a(1)*b(1)+a(2)*b(2)+a(3)*b(3))
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine selisamp(isample,nguess,maxcomp,dgrad,ndg,ind,maxc)
      implicit logical (a-z)
      integer nguess,maxcomp,ndg,isample(maxcomp,nguess)
      real*8  dgrad(ndg,ndg),maxc
      logical ind(nguess)
      integer i,j,k
      DO i=1,nguess
         ind(i)=.TRUE.
         DO j=1,maxcomp-1
            DO k=j+1,maxcomp
               IF(dgrad(isample(j,i),isample(k,i)).gt.maxc) THEN
                  ind(i)=.FALSE.
                  goto 1
               END IF
            END DO
         END DO
1        continue
      END DO
      RETURN
      END
C
C _________________________________________________________________
C
      subroutine zerofill(a,n)
      implicit logical(a-z)
      integer i,n
      real*8 a(n),ZERO
      PARAMETER   ( ZERO = 0.0D+0 )
      DO i=1,n
         a(i) = ZERO 
      END DO
      RETURN
      END
C
C _________________________________________________________________
C
      subroutine dcprod0(a,b,n,c)
C
C   compute component wise product of a and b in c
C
      implicit logical(a-z)
      integer i,n
      real*8 a(n),b(n),c(n)
      DO i=1,n
         c(i) = a(i)*b(i)
      END DO
      RETURN
      END
C
C _________________________________________________________________
C
      subroutine dcprod(a,b,alpha,n,c)
C
C   compute component wise product of a and b in c
C
      implicit logical(a-z)
      integer i,n
      real*8 a(n),b(n),c(n),alpha
      DO i=1,n
         c(i) = a(i)*b(i)*alpha
      END DO
      RETURN
      END
     
