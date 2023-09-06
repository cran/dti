C
C __________________________________________________________________
C
      subroutine mfunpl(par,w,siq,g,m,lpar,n,z,erg)
C
C   model without isotropic compartment
C   with explicit parametrization of weights
C
C   code is restricted to m<=6
C
      implicit none
      integer m,lpar,n
      double precision par(lpar),w(m),siq(n),g(3,n),z(n,m),erg
      integer i,j
      double precision th,sw,sth,z1,p0,p1,d1,d2,d3
      th = par(1)
      DO i = 1,m
         p0=par(2*i)
         p1=par(2*i+1)
         sth = sin(p0)
         d1 = sth*cos(p1)
         d2 = sth*sin(p1)
         d3 = cos(p0)
         DO j = 1,n
            z1 = d1*g(1,j)+d2*g(2,j)+d3*g(3,j)
            z(j,i) = exp(-th*z1*z1)
         END DO
      END DO
      sw=0.d0
      DO j=1,n
         z1 = siq(j)
         DO i=1,m
            z1=z1-w(i)*z(j,i)
         END DO
         sw=sw+z1*z1
      END DO
      erg=sw
      call rchkusr()
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine gmfunpl(par,w,siq,g,m,lpar,n,z,res,resd,dkgj,dkgj2,
     1                 ddkdphig,ddkdetag,dzdpars,work1,work2,dfdparw)
C
C   analytical gradients for c(par,w)
C   model without isotropic compartment
C   with explicit parametrization of weights
C
C   code is restricted to m<=6
C
      implicit none
      integer m,lpar,n
      double precision par(lpar),w(m),siq(n),g(3,n),z(n,m),res(n),
     1       resd(n),dkgj(n,m),dkgj2(n,m),ddkdphig(n,m),ddkdetag(n,m),
     2       dzdpars(n,m,3),work1(n,m),work2(n,m),dfdparw(*)
      integer i,j
      double precision th,m2th,sphi,cphi,seta,ceta,z1,z2,p0,p1,d1,d2,
     1       dphi1,dphi2
      double precision ddot
      external ddot
      th = par(1)
      m2th = -2.d0*th
      DO i = 1,m
         p0=par(2*i)
         p1=par(2*i+1)
         sphi = sin(p0)
         cphi = cos(p0)
         seta = sin(p1)
         ceta = cos(p1)
         d1 = sphi*ceta
         d2 = sphi*seta
         dphi1 = cphi*ceta
         dphi2 = cphi*seta
         DO j = 1,n
            z1 = d1*g(1,j)+d2*g(2,j)+cphi*g(3,j)
            dkgj(j,i) = z1
            z2 = z1*z1
            dkgj2(j,i) = z2
            z(j,i) = exp(-th*z2)
            ddkdphig(j,i) = dphi1*g(1,j)+dphi2*g(2,j)-sphi*g(3,j)
            ddkdetag(j,i) = d1*g(2,j)-d2*g(1,j)
         END DO
      END DO
      call dcprod0(dkgj,ddkdphig,n*m,work1)
      call dcprod0(dkgj,ddkdetag,n*m,work2)
C initialize unneeded elements
C      call zerofill(dzdpars,m*3*n)
C compute componentswise -dkgj2*z
      call dcprod(dkgj2,z,-1.d0,n*m,dzdpars(1,1,1))
C compute componentswise m2th*work1*z
      call dcprod(work1,z,m2th,n*m,dzdpars(1,1,2))
C compute componentswise m2th*work2*z
      call dcprod(work2,z,m2th,n*m,dzdpars(1,1,3))
C   Compute residuals
      DO j=1,n
         z1 = siq(j)
         z2 = 0.d0
         DO i=1,m
            z1 = z1 - w(i)*z(j,i)
            z2 = z2 + w(i)*dzdpars(j,i,1)
         END DO
         res(j)=z1
C         resd(j)=ddot(m,w,n,dzdpars(j,1,1),n)
         resd(j)=z2
      END DO
C Now compute gradient
      dfdparw(1)=-2.d0*ddot(n,res,1,resd,1)
      DO i=1,m
         dfdparw(2*i)=-2.d0*w(i)*ddot(n,res,1,dzdpars(1,i,2),1)
         dfdparw(2*i+1)=-2.d0*w(i)*ddot(n,res,1,dzdpars(1,i,3),1)
         dfdparw(2*m+1+i)=-2.d0*ddot(n,res,1,z(1,i),1)
      END DO
      call rchkusr()
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine mfunpli(par,w,siq,g,m,lpar,n,z,erg)
C
C   model with isotropic compartment
C   with explicit parametrization of weights in w(m+1)
C
C   code is restricted to m<=6
C
      implicit none
      integer m,lpar,n
      double precision par(lpar),w(*),siq(n),g(3,n),z(n,m),erg
      integer i,j
      double precision th,sw,sth,z1,p0,p1,d1,d2,d3
      th = par(1)
      DO i = 1,m
         p0=par(2*i)
         p1=par(2*i+1)
         sth = sin(p0)
         d1 = sth*cos(p1)
         d2 = sth*sin(p1)
         d3 = cos(p0)
         DO j = 1,n
            z1 = d1*g(1,j)+d2*g(2,j)+d3*g(3,j)
            z(j,i) = exp(-th*z1*z1)
         END DO
      END DO
      sw=0.d0
      DO j=1,n
         z1 = siq(j)-w(1)
         DO i=1,m
            z1=z1-w(i+1)*z(j,i)
         END DO
         sw=sw+z1*z1
      END DO
      erg=sw
      call rchkusr()
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine gmfunpli(par,w,siq,g,m,lpar,n,z,res,resd,dkgj,
     1           dkgj2,ddkdphig,ddkdetag,dzdpars,work1,work2,dfdparw)
C
C   analytical gradients for c(par,w)
C   model with isotropic compartment
C   with explicit parametrization of weights
C
C   code is restricted to m<=6
C
      implicit none
      integer m,lpar,n
      double precision par(lpar),w(m+1),siq(n),g(3,n),z(n,m),res(n),
     1       resd(n),dkgj(n,m),dkgj2(n,m),ddkdphig(n,m),ddkdetag(n,m),
     3       dzdpars(n,m,3),work1(n,m),work2(n,m),dfdparw(*)
      integer i,j
      double precision th,m2th,sphi,cphi,seta,ceta,z1,z2,p0,p1,sres,
     1       d1,d2,dphi1,dphi2
      double precision ddot
      external ddot
      th = par(1)
      m2th = -2.d0*th
      DO i = 1,m
         p0=par(2*i)
         p1=par(2*i+1)
         sphi = sin(p0)
         cphi = cos(p0)
         seta = sin(p1)
         ceta = cos(p1)
         d1 = sphi*ceta
         d2 = sphi*seta
         dphi1 = cphi*ceta
         dphi2 = cphi*seta
         DO j = 1,n
            z1 = d1*g(1,j)+d2*g(2,j)+cphi*g(3,j)
            dkgj(j,i) = z1
            z2 = z1*z1
            dkgj2(j,i) = z2
            z(j,i) = exp(-th*z2)
            ddkdphig(j,i) = dphi1*g(1,j)+dphi2*g(2,j)-sphi*g(3,j)
            ddkdetag(j,i) = d1*g(2,j)-d2*g(1,j)
         END DO
      END DO
      call dcprod0(dkgj,ddkdphig,n*m,work1)
      call dcprod0(dkgj,ddkdetag,n*m,work2)
C initialize unneeded elements
C      call zerofill(dzdpars,m*3*n)
C compute componentswise -dkgj2*z
      call dcprod(dkgj2,z,-1.d0,n*m,dzdpars(1,1,1))
C compute componentswise m2th*work1*z
      call dcprod(work1,z,m2th,n*m,dzdpars(1,1,2))
C compute componentswise m2th*work2*z
      call dcprod(work2,z,m2th,n*m,dzdpars(1,1,3))
C   Compute residuals
      sres=0.d0
      DO j=1,n
         z1 = siq(j)-w(1)
         z2 = 0.d0
         DO i=1,m
            z1 = z1 - w(i+1)*z(j,i)
            z2 = z2 + w(i+1)*dzdpars(j,i,1)
         END DO
         res(j) = z1
         resd(j) = z2
         sres=sres+z1
      END DO
C Now compute gradient
      dfdparw(1)=-2.d0*ddot(n,res,1,resd,1)
      dfdparw(2*m+2)=-2.d0*sres
      DO i=1,m
         dfdparw(2*i)=-2.d0*w(i+1)*ddot(n,res,1,dzdpars(1,i,2),1)
         dfdparw(2*i+1)=-2.d0*w(i+1)*ddot(n,res,1,dzdpars(1,i,3),1)
         dfdparw(2*m+2+i)=-2.d0*ddot(n,res,1,z(1,i),1)
      END DO
      call rchkusr()
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine mfunpl0(par,siq,g,m,lpar,n,pen,z,w,erg)
C
C   model without isotropic compartment
C   same as mfunpl but with unconstrained least squares for weights
C
C   code is restricted to m<=6
C
      implicit none
      integer m,lpar,n
      double precision par(lpar),siq(n),g(3,n),z(n,m),erg,pen
      integer i,j,mode,jpvt(6),rank
      double precision th,w(n),sw,sth,z1,p0,p1,d1,d2,d3,work(25)
      th = par(1)
      th = max(th,-5.d0)
      DO i = 1,m
         p0=par(2*i)
         p1=par(2*i+1)
         sth = sin(p0)
         d1 = sth*cos(p1)
         d2 = sth*sin(p1)
         d3 = cos(p0)
         DO j = 1,n
            z1 = d1*g(1,j)+d2*g(2,j)+d3*g(3,j)
            z(j,i) = exp(-th*z1*z1)
         END DO
         jpvt(i)=0
      END DO
C
C    siq will be replaced, need to copy it if C-version of optim is used
C
      call dcopy(n,siq,1,w,1)
      call dgelsy(n,m,1,z,n,w,n,jpvt,1d-8,rank,work,25,
     1            mode)
C  1d-6  is a limit for condition number
      IF(mode.eq.0) THEN
         sw=0.d0
C penalize for extreme th values
         if(th.gt.1.d1) sw=sw+th-1.d1
C penalize for negative weights
         if(th.lt.1.d-2) sw=sw-1.d2*th+1.d0
         DO i=1,m
            if(w(i).lt.0.d0) sw=sw-pen*w(i)
         END DO
         DO i=m+1,n
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
      subroutine mfunpl0w(par,w,siq,g,m,lpar,n,z,erg)
C
C   model without isotropic compartment
C   same as mfunpl but with unconstrained least squares for weights
C
C   code is restricted to m<=6
C
      implicit none
      integer m,lpar,n
      double precision par(lpar),siq(n),g(3,n),z(n,m),w(m)
      integer i,j
      double precision th,sth,z1,p0,p1,d1,d2,d3,rss,res,erg
      th = par(1)
      th = max(th,-5.d0)
      DO i = 1,m
         p0=par(2*i)
         p1=par(2*i+1)
         sth = sin(p0)
         d1 = sth*cos(p1)
         d2 = sth*sin(p1)
         d3 = cos(p0)
         DO j = 1,n
            z1 = d1*g(1,j)+d2*g(2,j)+d3*g(3,j)
            z(j,i) = exp(-th*z1*z1)
         END DO
      END DO
      rss =0.d0
      DO j=1,n
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
      subroutine mfpl0gn(par,siq,g,m,lpar,n,pen,eps,z,w,para,parb,
     1                   dfdpar)
C
C   model without isotropic compartment
C   same as mfunpl but with unconstrained least squares for weights
C
C   code is restricted to m<=6
C
      implicit none
      integer m,lpar,n
      double precision par(lpar),siq(n),g(3,n),z(n,m),pen,dfdpar(lpar),
     1       w(n),para(lpar),parb(lpar),eps
      double precision erga,ergb,deltai
      integer i
      deltai=0.5d0/eps
      DO i=1,lpar
         call dcopy(lpar,par,1,para,1)
         call dcopy(lpar,par,1,parb,1)
         para(i)=para(i)-eps
         parb(i)=parb(i)+eps
         call mfunpl0(para,siq,g,m,lpar,n,pen,z,w,erga)
         call mfunpl0(parb,siq,g,m,lpar,n,pen,z,w,ergb)
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
      subroutine mfunpl0g(par,s,g,m,lpar,n,z,v,w,dkgj,dkgj2,ddkdphig,
     1                  ddkdetag,dvdth,dvdphi,dvdeta,dzdpars,dwdpars,
     2                  dwdpars2,zs,work1,work2,scopy,pen,dfdpar)
      implicit none
      integer m,n,lpar
      double precision par(lpar),s(n),g(3,n),z(n,m),v(m,m),dkgj(n,m),
     1       w(n),dkgj2(n,m),ddkdphig(n,m),ddkdetag(n,m),
     2       dvdth(m,m),dvdphi(m,m,m),dvdeta(m,m,m),
     3       dzdpars(n,m,3),dwdpars(m,lpar),dwdpars2(m,lpar),
     4       zs(n,m),dfdpar(lpar),pen,rcond,ferr(11),berr(11)
      integer i,j,k,l,ind(5),mode,jpvt(5),rank
      double precision th,sphi,cphi,seta,ceta,z1,z2,p0,p1,d1,d2,dphi1,
     1       dphi2,work(250),work1(n,m),work2(n,m),scopy(n),m2th,af(25)
C  length of work needs to be larger or equal max(m*m,ngrad)
      double precision ddot
      external ddot
      th = par(1)
      th = max(-5.d0,th)
      m2th = -2.d0*th
C
C  get d, dkgj, dkgj2, dddphi, dddeta, z, ddkdphig, ddkdetag
C
      DO i = 1,m
         p0=par(2*i)
         p1=par(2*i+1)
         sphi = sin(p0)
         cphi = cos(p0)
         seta = sin(p1)
         ceta = cos(p1)
         d1 = sphi*ceta
         d2 = sphi*seta
C         d3 = cphi
         dphi1 = cphi*ceta
         dphi2 = cphi*seta
C         dphi3 = -sphi
C         deta1 = -sphi*seta
C         deta2 = sphi*ceta
         DO j = 1,n
            z1 = d1*g(1,j)+d2*g(2,j)+cphi*g(3,j)
            dkgj(j,i) = z1
            z2 = z1*z1
            dkgj2(j,i) = z2
            z(j,i) = exp(-th*z2)
            zs(j,i) = z(j,i)*s(j)
            ddkdphig(j,i) = dphi1*g(1,j)+dphi2*g(2,j)-sphi*g(3,j)
            ddkdetag(j,i) = d1*g(2,j)-d2*g(1,j)
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
      call dgelsy(n,m,1,work1,n,w,n,jpvt,1d-8,rank,work,25,mode)
C  1d-6  is a limit for condition number
C
C   thats weights in w now V, dVdth, dVdphi, dVdeta, dzdpars
C
C   use work1, work2 and scopy for intermediate results
C   componentwise products
      call dcprod0(dkgj,ddkdphig,n*m,work1)
      call dcprod0(dkgj,ddkdetag,n*m,work2)
C initialize unneeded elements
C      call zerofill(dzdpars,m*3*n)
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
        call mfpl0gn(par,s,g,m,lpar,n,pen,1d-8,z,w,work1,work2,dfdpar)
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
      implicit none
      integer nvox,ngrad,ns,siind(ns,nvox),m,ntry,nth,nv,
     1       isample(m,ntry),indth(nvox)
      double precision si(ngrad,nvox),sms(ngrad),dgrad(ngrad,nv),
     1       th(nth),egrad(ngrad,nv),z(ngrad,ns),mval(nvox),
     2       vsi(nvox)
      integer mask(nvox)
      integer i,k,ibest,mode,ind(10),l,j,ii,iw,wind(5),nwi(5)
      double precision w(1000),krit,work1(1000),work2(10),erg,thj,
     1       msi,m2si,z1,dng
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
         if(mask(i).eq.0) THEN
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
            if(mask(i).ne.0) THEN
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
                  IF(mode.le.1) THEN
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
                  IF (iw.ge.1) THEN
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
      subroutine pgtsii30(si,vsi,ngrad,nvox,m,dgrad,nv,th,
     1     nth,indth,egrad,isample,ntry,sms,z,siind,mval,ns)
C
C  compute diagnostics for initial estimates in siind
C  siind(1,i1,i2,i3) will contain the model order (mc parallel version)
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
      implicit none
      integer nvox,ngrad,ns,siind(ns,nvox),m,ntry,nth,nv,
     1       isample(m,ntry),indth(nvox)
      double precision si(ngrad,nvox),sms(ngrad),dgrad(ngrad,nv),
     1       th(nth),egrad(ngrad,nv),z(ngrad,ns),mval(nvox),
     2       vsi(nvox)
      integer i,k,ibest,mode,ind(10),l,j,ii,iw,wind(5),nwi(5)
      double precision w(1000),krit,work1(1000),work2(10),erg,thj,
     1       msi,m2si,z1,dng
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
               IF(mode.le.1) THEN
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
               IF (iw.ge.1) THEN
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
      implicit none
      integer nvox,ngrad,ns,siind(ns,nvox),m,ntry,nth,nv,
     1       isample(*),indth(nvox),iandir(nvox),mask(nvox)
      double precision si(ngrad,nvox),sms(ngrad),dgrad(ngrad,nv),
     1       th(nth),egrad(ngrad,nv),z(ngrad,ns),mval(nvox),
     2       vsi(nvox),dgradv(nv,nv),maxc
      logical skip
      integer i,k,mode,ind(10),l,j,ii,iw,wind(5),nwi(5),mis,
     1        is(5),isbest(5),ntry0,km1mis
      double precision w(1000),krit,work1(1000),work2(10),erg,thj,
     1       msi,m2si,z1,dng
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
         if(mask(i).eq.0) THEN
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
            if(mask(i).ne.0) THEN
               IF(j.ne.indth(i)) CYCLE
C  now search for minima of sms (or weighted sms)
               krit=mval(i)
               DO k=1,ntry0
                  km1mis=(k-1)*mis
                  IF(m.gt.1) THEN
                     skip=.FALSE.
                     DO l=1,m-1
                   if(dgradv(isample(km1mis+l),iandir(i)).gt.maxc)
     1                    skip=.TRUE.
                     END DO
                     IF(skip) CYCLE
                  END IF
                  call dcopy(ngrad,si(1,i),1,sms,1)
                  IF(m.gt.1) THEN
                     DO l=1,m-1
                        is(l)=isample(km1mis+l)
                        call dcopy(ngrad,egrad(1,is(l)),1,z(1,l),1)
                     END DO
                  END IF
                  is(m)=iandir(i)
                  call dcopy(ngrad,egrad(1,is(m)),1,z(1,m),1)
        call nnls(z,ngrad,ngrad,m,sms,w,erg,work2,work1,ind,mode)
                  IF(mode.le.1) THEN
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
               IF (iw.ge.1) THEN
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
      subroutine pgtsii31(si,vsi,ngrad,nvox,m,dgrad,nv,iandir,th,
     1     nth,indth,egrad,isample,ntry,sms,z,siind,mval,ns,
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
C  mval   - aktual best risk
C
C  restricted to ngrad<=1000 and m <=10
C
      implicit none
      integer nvox,ngrad,ns,siind(ns,nvox),m,ntry,nth,nv,
     1       isample(*),indth(nvox),iandir(nvox)
      double precision si(ngrad,nvox),sms(ngrad),dgrad(ngrad,nv),
     1       th(nth),egrad(ngrad,nv),z(ngrad,ns),mval(nvox),
     2       vsi(nvox),dgradv(nv,nv),maxc
      logical skip
      integer i,k,mode,ind(10),l,j,ii,iw,wind(5),nwi(5),mis,
     1        is(5),isbest(5),ntry0,km1mis
      double precision w(1000),krit,work1(1000),work2(10),erg,thj,
     1       msi,m2si,z1,dng
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
            IF(j.ne.indth(i)) CYCLE
C  now search for minima of sms (or weighted sms)
            krit=mval(i)
            DO k=1,ntry0
               km1mis=(k-1)*mis
               IF(m.gt.1) THEN
                  skip=.FALSE.
                  DO l=1,m-1
                   if(dgradv(isample(km1mis+l),iandir(i)).gt.maxc)
     1                    skip=.TRUE.
                  END DO
                  IF(skip) CYCLE
               END IF
               call dcopy(ngrad,si(1,i),1,sms,1)
               IF(m.gt.1) THEN
                  DO l=1,m-1
                     is(l)=isample(km1mis+l)
                     call dcopy(ngrad,egrad(1,is(l)),1,z(1,l),1)
                  END DO
               END IF
               is(m)=iandir(i)
               call dcopy(ngrad,egrad(1,is(m)),1,z(1,m),1)
        call nnls(z,ngrad,ngrad,m,sms,w,erg,work2,work1,ind,mode)
               IF(mode.le.1) THEN
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
            IF (iw.ge.1) THEN
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
            call rchkusr()
         END DO
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine iandir(vico,nvico,andir,nvox,landir,iandi)
      implicit none
      integer nvico,nvox,iandi(nvox)
      double precision vico(3,nvico),andir(3,2,nvox)
      integer landir(nvox)
      integer i,j,jmax
      double precision z,zmax,scprod3
      external scprod3
      DO i=1,nvox
         if(landir(i).ne.0) THEN
            zmax = scprod3(vico(1,1),andir(1,1,i))
            jmax = 1
            DO j=2,nvico
               z = scprod3(vico(1,j),andir(1,1,i))
               if(z.gt.zmax) THEN
                  zmax=z
                  jmax=j
               END IF
            END DO
            iandi(i)=jmax
         END IF
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      double precision function scprod3(a,b)
      implicit none
      double precision a(3),b(3)
      scprod3=abs(a(1)*b(1)+a(2)*b(2)+a(3)*b(3))
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine selisamp(isample,nguess,maxcomp,dgrad,ndg,ind,maxc)
      implicit none
      integer nguess,maxcomp,ndg,isample(maxcomp,nguess)
      double precision  dgrad(ndg,ndg),maxc
      integer ind(nguess)
      integer i,j,k
      DO i=1,nguess
         ind(i)=1
         DO j=1,maxcomp-1
            DO k=j+1,maxcomp
               IF(dgrad(isample(j,i),isample(k,i)).gt.maxc) THEN
                  ind(i)=0
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
      implicit none
      integer i,n
      double precision a(n),ZERO
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
      implicit none
      integer i,n
      double precision a(n),b(n),c(n)
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
      implicit none
      integer i,n
      double precision a(n),b(n),c(n),alpha
      DO i=1,n
         c(i) = a(i)*b(i)*alpha
      END DO
      RETURN
      END
