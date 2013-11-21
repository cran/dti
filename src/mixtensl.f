C
C  version for use with L-BFGS-B
C
      subroutine rskmixl2(par,npar,siq,g,b,ng,risk)
C
C    tensor-mixture models with lambda2 and fa as parameters
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit logical(a-z)
      integer npar,ng
      real*8 par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,risk
      integer i,ncomp
      real*8 resi,fval
      ncomp=(npar-2)/3
      lambda=par(3*ncomp+1)
      alpha=par(3*ncomp+2)
      risk=0.d0
      DO i=1,ng
         call fmixturl(par,npar-2,lambda,alpha,g(1,i),b(i),fval)
         resi=siq(i)-fval
         risk=risk+resi*resi
      END DO
      RETURN
      END
      subroutine rskmixl1(par,npar,siq,g,b,ng,alpha,risk)
C
C    tensor-mixture models with lambda2 as parameter, fa fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit logical(a-z)
      integer npar,ng
      real*8 par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,risk
      integer i,ncomp
      real*8 resi,fval
      ncomp=(npar-1)/3
      lambda=par(3*ncomp+1)
      risk=0.d0
      DO i=1,ng
         call fmixturl(par,npar-1,lambda,alpha,g(1,i),b(i),fval)
         resi=siq(i)-fval
         risk=risk+resi*resi
      END DO
      RETURN
      END
      subroutine rskmixl0(par,npar,siq,g,b,ng,lambda,alpha,risk)
C
C    tensor-mixture models with eigenvalues fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit logical(a-z)
      integer npar,ng
      real*8 par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,risk
      integer i
      real*8 resi,fval
      risk=0.d0
      DO i=1,ng
         call fmixturl(par,npar,lambda,alpha,g(1,i),b(i),fval)
         resi=siq(i)-fval
         risk=risk+resi*resi
      END DO
      RETURN
      END
C
C    Evaluations for tensor-mixture models
C
      subroutine fmixturl(par,npar,lambda,alpha,g,b,fval)
C
C  compute expected value
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda, c
C
      implicit logical (a-z)
      integer npar
      real*8 par(npar),g(3),fval,lambda,alpha
      integer ncomp,i
      real*8 w(5),phi(5),eta(5),b,wi,sw,blam,f1
      real*8 ddot3sq
      external ddot3sq
C  extract parameters
      ncomp=npar/3
      sw=1.d0
      DO i=1,ncomp
         wi=exp(par(3*i-2))
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
      blam=b*lambda
      f1=exp(-blam*(alpha+1.d0))
      DO i=1,ncomp
         f1=f1+w(i)*exp(-blam*(1.d0+alpha*ddot3sq(phi(i),eta(i),g)))
      END DO
      fval=f1/sw
      RETURN
      END
      subroutine drskml2(par,npar,siq,g,b,ng,drisk)
C
C    tensor-mixture models with lambda2 and fa as parameters
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit logical(a-z)
      integer npar,ng
      real*8 par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,drisk(npar)
      integer i,j,ncomp
      real*8 resi,fval,dval(15),drisk0(17),dlam,dalpha
C  calculate numerical gradients for comparison
      ncomp=(npar-2)/3
      lambda=par(3*ncomp+1)
      alpha=par(3*ncomp+2)
      DO j=1,npar
         drisk0(j)=0.d0   
      END DO
      DO i=1,ng
         call dfml2(par,npar,lambda,alpha,g(1,i),b(i),fval,dval,
     1                  dlam,dalpha)         
         resi=siq(i)-fval
         DO j=1,3*ncomp
            drisk0(j)=drisk0(j)-resi*dval(j)
         END DO
         drisk0(npar-1)=drisk0(npar-1)-resi*dlam
         drisk0(npar)=drisk0(npar)-resi*dalpha
      END DO
      DO j=1,npar
         drisk(j)=2.d0*drisk0(j)   
      END DO
      RETURN
      END
C
C  gradients
C
      subroutine drskml1(par,npar,siq,g,b,ng,alpha,drisk)
C
C    tensor-mixture models with lambda2 as parameter, fa fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit logical(a-z)
      integer npar,ng
      real*8 par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,drisk(npar)
      integer i,j,ncomp
      real*8 resi,fval,dval(15),drisk0(17),dlam
      ncomp=(npar-1)/3
      lambda=par(3*ncomp+1)
      DO j=1,npar
         drisk0(j)=0.d0   
      END DO
      DO i=1,ng
         call dfml1(par,npar,lambda,alpha,g(1,i),b(i),fval,dval,dlam)
         resi=siq(i)-fval
         DO j=1,3*ncomp
            drisk0(j)=drisk0(j)-resi*dval(j)
         END DO
         drisk0(npar)=drisk0(npar)-resi*dlam
      END DO
      DO j=1,npar
         drisk(j)=2.d0*drisk0(j)   
      END DO
      RETURN
      END
      subroutine drskml0(par,npar,siq,g,b,ng,lambda,alpha,drisk)
C
C    tensor-mixture models with eigenvalues fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit logical(a-z)
      integer npar,ng
      real*8 par(npar),g(3,ng),siq(ng),b(ng),lambda,alpha,drisk(npar)
      integer i,j,ncomp
      real*8 resi,fval,dval(15),drisk0(17)
      ncomp=npar/3
      DO j=1,npar
         drisk0(j)=0.d0   
      END DO
      DO i=1,ng
         call dfml0(par,npar,lambda,alpha,g(1,i),b(i),fval,dval)
         resi=siq(i)-fval
         DO j=1,3*ncomp
            drisk0(j)=drisk0(j)-resi*dval(j)
         END DO
      END DO
      DO j=1,npar
         drisk(j)=2.d0*drisk0(j)   
      END DO
      RETURN
      END
      subroutine dfml2(par,npar,lambda,alpha,g,b,fval,dval,dlam,
     1                  dalpha)
C
C  compute partial derivatives with respect to w_i, phi_i, eta_i, lambda and alpha
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda, alpha
C  fval contains function value
      implicit logical (a-z)
      integer npar
      real*8 par(npar),g(3),dval(npar),alpha,lambda,b,dlam,dalpha
      integer ncomp,i
      real*8 w(5),phi(5),eta(5),wi,sw,blam
      real*8 embclgd2(5),dgtd(3,5),fval,f1,f1l,f2l,f2c,adgtd(5)
C  extract parameters
      ncomp=(npar-2)/3
C need to have positive alpha to guarantee lambda1 >= lambda2
      sw=1.d0
      DO i=1,ncomp
         wi=exp(par(3*i-2))
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
C  precompute common terms
      blam=b*lambda
      f1=exp(-blam*(alpha+1.d0))
      fval = f1
      DO i=1,ncomp
         call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
         adgtd(i)=1.d0+alpha*dgtd(1,i)
         embclgd2(i)=exp(-blam*adgtd(i))
         fval=fval+w(i)*embclgd2(i)
      END DO
C  compute function value
      fval=fval/sw
C  compute derivatives
      f1l=-b*(1.d0+alpha)/sw*f1
      f2l=0.d0
      f2c=0.d0
      DO i=1,ncomp
         f2l=f2l+w(i)*embclgd2(i)*adgtd(i)
         f2c=f2c+w(i)*embclgd2(i)*dgtd(1,i)
         dval(3*i-1)=-alpha*blam*w(i)*embclgd2(i)*dgtd(2,i)/sw
C  derivative with respect to phi_i
         dval(3*i)=-alpha*blam*w(i)*embclgd2(i)*dgtd(3,i)/sw
C  derivative with respect to eta_i
         dval(3*i-2)=(embclgd2(i)-fval)/sw*w(i)
C  derivative with respect to eta_i
      END DO
      dlam=f1l-b/sw*f2l
C  derivative with respect to lambda
      dalpha=-blam/sw*(f1+f2c)
C  derivative with respect to c
      RETURN
      END
      subroutine dfml1(par,npar,lambda,alpha,g,b,fval,dval,dlam)
C
C  compute partial derivatives with respect to w_i, phi_i, eta_i and lambda
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda
C  fval contains function value
C
      implicit logical (a-z)
      integer npar
      real*8 par(npar),g(3),dval(npar),alpha,lambda,b,dlam
      integer ncomp,i
      real*8 w(5),phi(5),eta(5),wi,sw,blam
      real*8 embclgd2(5),dgtd(3,5),fval,f1,f1l,f2l,adgtd(5)
C  extract parameters
      ncomp=(npar-1)/3
C need to have positive alpha to guarantee lambda1 >= lambda2
      sw=1.d0
      DO i=1,ncomp
         wi=exp(par(3*i-2))
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
C  precompute common terms
      blam=b*lambda
      f1=exp(-blam*(alpha+1.d0))
      fval = f1
      DO i=1,ncomp
         call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
         adgtd(i)=1.d0+alpha*dgtd(1,i)
         embclgd2(i)=exp(-blam*adgtd(i))
         fval=fval+w(i)*embclgd2(i)
      END DO
C  compute function value
      fval=fval/sw
C  compute derivatives
      f1l=-b*(1.d0+alpha)/sw*f1
      f2l=0.d0
      DO i=1,ncomp
         f2l=f2l+w(i)*embclgd2(i)*adgtd(i)
         dval(3*i-1)=-alpha*blam*w(i)*embclgd2(i)*dgtd(2,i)/sw
C  derivative with respect to phi_i
         dval(3*i)=-alpha*blam*w(i)*embclgd2(i)*dgtd(3,i)/sw
C  derivative with respect to eta_i
         dval(3*i-2)=(embclgd2(i)-fval)/sw*w(i)
C  derivative with respect to eta_i
      END DO
      dlam=f1l-b/sw*f2l
C  derivative with respect to lambda
      RETURN
      END
      subroutine dfml0(par,npar,lambda,alpha,g,b,fval,dval)
C
C  compute partial derivatives with respect to w_i, phi_i, eta_i and lambda
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda
C  fval contains function value
C
      implicit logical (a-z)
      integer npar
      real*8 par(npar),g(3),dval(npar),alpha,lambda,b
      integer ncomp,i
      real*8 w(5),phi(5),eta(5),wi,sw,blam
      real*8 embclgd2(5),dgtd(3,5),fval,adgtd(5)
C  extract parameters
      ncomp=npar/3
C need to have positive alpha to guarantee lambda1 >= lambda2
      sw=1.d0
      DO i=1,ncomp
         wi=exp(par(3*i-2))
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
C  precompute common terms
      blam=b*lambda
      fval = exp(-blam*(alpha+1.d0))
      DO i=1,ncomp
         call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
         adgtd(i)=1.d0+alpha*dgtd(1,i)
         embclgd2(i)=exp(-blam*adgtd(i))
         fval=fval+w(i)*embclgd2(i)
      END DO
C  compute function value
      fval=fval/sw
C  compute derivatives
      DO i=1,ncomp
         dval(3*i-1)=-alpha*blam*w(i)*embclgd2(i)*dgtd(2,i)/sw
C  derivative with respect to phi_i
         dval(3*i)=-alpha*blam*w(i)*embclgd2(i)*dgtd(3,i)/sw
C  derivative with respect to eta_i
         dval(3*i-2)=(embclgd2(i)-fval)/sw*w(i)
C  derivative with respect to eta_i
      END DO
      RETURN
      END
      
      real*8 function ddot3sq(phi,eta,g)
C  compute (g^T n(phi,eta))^2 
      implicit logical (a-z)
      real*8 phi,eta,g(3),sphi,z
      sphi=dsin(phi)
      z=g(1)*dcos(eta)*sphi+g(2)*dsin(eta)*sphi+g(3)*dcos(phi)
      ddot3sq=z*z
      RETURN
      END
      subroutine dgtddphi(phi,eta,g, d)
C  compute (g^T n(phi,eta))^2 
C           d(g^T n(phi,eta))^2 / d phi
C           d(g^T n(phi,eta))^2 / d eta
C  as components of d
      implicit logical (a-z)
      real*8 phi,eta,g(3),sphi,cphi,ceta,seta,z,ze,zp,d(3)
      sphi=dsin(phi)
      cphi=dcos(phi)
      seta=dsin(eta)
      ceta=dcos(eta)
      z=g(1)*ceta*sphi+g(2)*seta*sphi+g(3)*cphi
      zp=g(1)*ceta*cphi+g(2)*seta*cphi-g(3)*sphi
      ze=g(2)*ceta*sphi-g(1)*seta*sphi
      d(1)=z*z
      d(2)=2.d0*z*zp
      d(3)=2.d0*z*ze
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine getsii(si,vsi,ngrad,nvox,m,dgrad,bv,nv,alpha,
     1        lambda,egrad,isample,ntry,sms,z0,z,siind,mval,ns)
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
      integer nvox,ngrad,ns,siind(ns,nvox),m,ntry,nv,
     1       isample(m,ntry)
      real*8 si(ngrad,nvox),sms(ngrad),dgrad(ngrad,nv),
     1       egrad(ngrad,nv),z(ngrad,ns),mval(nvox),
     2       vsi(nvox),bv(ngrad),alpha,lambda,z0(ngrad)
      integer i,k,ibest,mode,ind(10),l,ii,iw,wind(6),nwi(6)
      real*8 w(1000),krit,work1(1000),work2(12),erg,msi,m2si,
     1       z1,dng,albv,lbv
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
      DO k=1,ngrad
         lbv = lambda*bv(k)
         z0(k) = dexp(-lbv*(1.d0+alpha))
         DO l=1,nv
            albv = alpha*lbv
            z1 = dgrad(k,l)
            egrad(k,l)=dexp(-lbv-albv*z1*z1)
         END DO
      END DO
      DO i=1,nvox
C  now search for minima of sms (or weighted sms
         ibest=0
         krit=mval(i)
         DO k=1,ntry
            call dcopy(ngrad,si(1,i),1,sms,1)
            call dcopy(ngrad,z0,1,z(1,1),1)
            DO l=1,m
               call dcopy(ngrad,egrad(1,isample(l,k)),1,z(1,l+1),1)
            END DO
         if(i.eq.16) THEN
         END IF
         call nnls(z,ngrad,ngrad,m+1,sms,w,erg,work2,work1,ind,mode)
            IF(mode.gt.1) THEN
               call intpr("mode",4,mode,1)
               call intpr("isample",7,isample(1,k),m)
            ELSE 
               IF(erg.lt.krit) THEN
                  krit=erg
                  ibest=k
                  iw=0
                  DO ii=2,m+1
                     if(w(ii).gt.1.d-12) THEN
                        iw=iw+1
                        wind(iw)=ii-1
                     ELSE
                        nwi(ii-iw-1)=ii-1
C   nonactive directions
                     END IF 
                  END DO
               END IF  
            END IF
         END DO
         if(ibest.gt.0) THEN
            siind(1,i)=iw
            IF (iw.ge.1) THEN
               DO l=1,iw
                  siind(l+1,i)=isample(wind(l),ibest)
               END DO
            END IF
            IF (iw.lt.m) THEN
               DO l=1,m-iw
                  siind(m-l+2,i)=isample(nwi(l),ibest)
               END DO
            END IF
            mval(i)=krit
         END IF
      END DO
      RETURN
      END
      
