C
C  version for use with L-BFGS-B for general b-values
C
      subroutine rskmixb2(par,npar,si,g,b,ng,risk)
C
C    tensor-mixture models with lambda2 and fa as parameters
C     compute sum_i (si-f(par,g(,i),b(i))^2
C     differs from old version in using S_i instead of the E_i=S_i/S_0
C     parameter vector: w_1,phi_1,eta_1,...,w0,lambda,alpha
      implicit none
      integer npar,ng
      double precision par(npar),g(3,ng),si(ng),b(ng),lambda,alpha,
     1       risk,w0
      integer i
      double precision resi,fval
      lambda=par(npar-1)
      alpha=par(npar)
      w0=par(npar-2)
      risk=0.d0
      DO i=1,ng
         call fmixturb(par,npar-3,w0,lambda,alpha,g(1,i),b(i),fval)
         resi=si(i)-fval
         risk=risk+resi*resi
      END DO
      RETURN
      END
      subroutine rskmixb1(par,npar,si,g,b,ng,alpha,risk)
C
C    tensor-mixture models with lambda2 as parameter, fa fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit none
      integer npar,ng
      double precision par(npar),g(3,ng),si(ng),b(ng),lambda,alpha,
     1       risk,w0
      integer i
      double precision resi,fval
      lambda=par(npar)
      w0=par(npar-1)
      risk=0.d0
      DO i=1,ng
         call fmixturb(par,npar-2,w0,lambda,alpha,g(1,i),b(i),fval)
         resi=si(i)-fval
         risk=risk+resi*resi
      END DO
      RETURN
      END
      subroutine rskmixb0(par,npar,si,g,b,ng,lambda,alpha,risk)
C
C    tensor-mixture models with eigenvalues fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit none
      integer npar,ng
      double precision par(npar),g(3,ng),si(ng),b(ng),lambda,alpha,
     1       risk,w0
      integer i
      double precision resi,fval
      w0=par(npar)
      risk=0.d0
      DO i=1,ng
         call fmixturb(par,npar-1,w0,lambda,alpha,g(1,i),b(i),fval)
         resi=si(i)-fval
         risk=risk+resi*resi
      END DO
      RETURN
      END
C
C    Evaluations for tensor-mixture models
C
      subroutine fmixturb(par,npar,w0,lambda,alpha,g,b,fval)
C
C  compute expected value
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda, c
C
      implicit none
      integer npar
      double precision par(npar),g(3),fval,lambda,alpha,w0
      integer ncomp,i
      double precision w(5),phi(5),eta(5),b,blam,f1
      double precision ddot3sq
      external ddot3sq
C  extract parameters
      ncomp=npar/3
      DO i=1,ncomp
         w(i)=par(3*i-2)
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
      END DO
      blam=b*lambda
      f1=w0*exp(-blam*(alpha+1.d0))
      DO i=1,ncomp
         f1=f1+w(i)*exp(-blam*(1.d0+alpha*ddot3sq(phi(i),eta(i),g)))
      END DO
      fval=f1
      RETURN
      END
      subroutine drskmb2(par,npar,si,g,b,ng,drisk)
C
C    tensor-mixture models with lambda2 and fa as parameters
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit none
      integer npar,ng
      double precision par(npar),g(3,ng),si(ng),b(ng),lambda,alpha,
     1       drisk(npar)
      integer i,j
      double precision resi,fval,dval(15),drisk0(18),dlam,dalpha,w0,dw0
C  calculate numerical gradients for comparison
      w0=par(npar-2)
      lambda=par(npar-1)
      alpha=par(npar)
      DO j=1,npar
         drisk0(j)=0.d0
      END DO
      DO i=1,ng
         call dfmb2(par,npar,w0,lambda,alpha,g(1,i),b(i),fval,dval,
     1                  dw0,dlam,dalpha)
         resi=si(i)-fval
         DO j=1,npar-3
            drisk0(j)=drisk0(j)-resi*dval(j)
         END DO
         drisk0(npar-2)=drisk0(npar-2)-resi*dw0
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
      subroutine drskmb1(par,npar,si,g,b,ng,alpha,drisk)
C
C    tensor-mixture models with lambda2 as parameter, fa fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit none
      integer npar,ng
      double precision par(npar),g(3,ng),si(ng),b(ng),lambda,alpha,
     1       drisk(npar)
      integer i,j
      double precision resi,fval,dval(15),drisk0(17),dlam,w0,dw0
      w0=par(npar-1)
      lambda=par(npar)
      DO j=1,npar
         drisk0(j)=0.d0
      END DO
      DO i=1,ng
         call dfmb1(par,npar,w0,lambda,alpha,g(1,i),b(i),fval,dval,
     1                  dw0,dlam)
         resi=si(i)-fval
         DO j=1,npar-2
            drisk0(j)=drisk0(j)-resi*dval(j)
         END DO
         drisk0(npar-1)=drisk0(npar-1)-resi*dw0
         drisk0(npar)=drisk0(npar)-resi*dlam
      END DO
      DO j=1,npar
         drisk(j)=2.d0*drisk0(j)
      END DO
      RETURN
      END
      subroutine drskmb0(par,npar,si,g,b,ng,lambda,alpha,drisk)
C
C    tensor-mixture models with eigenvalues fixed
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit none
      integer npar,ng
      double precision par(npar),g(3,ng),si(ng),b(ng),lambda,alpha,
     1       drisk(npar)
      integer i,j
      double precision resi,fval,dval(15),drisk0(16),w0,dw0
      w0=par(npar)
      DO j=1,npar
         drisk0(j)=0.d0
      END DO
      DO i=1,ng
         call dfmb0(par,npar,w0,lambda,alpha,g(1,i),b(i),fval,dval,dw0)
         resi=si(i)-fval
         DO j=1,npar-1
            drisk0(j)=drisk0(j)-resi*dval(j)
         END DO
      END DO
      DO j=1,npar-1
         drisk(j)=2.d0*drisk0(j)
      END DO
      drisk0(npar)=drisk0(npar)-resi*dw0
      RETURN
      END
      subroutine dfmb2(par,npar,w0,lambda,alpha,g,b,fval,dval,dw0,
     1                  dlam,dalpha)
C
C  compute partial derivatives with respect to w_i, phi_i, eta_i, lambda and alpha
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda, alpha
C  fval contains function value
      implicit none
      integer npar
      double precision par(npar),g(3),dval(npar),w0,alpha,lambda,b,
     1       dw0,dlam,dalpha
      integer ncomp,i
      double precision w(5),phi(5),eta(5),blam
      double precision embclgd2(5),dgtd(3,5),fval,f1,f1l,f2l,f2c,
     1       adgtd(5)
C  extract parameters
      ncomp=(npar-3)/3
C need to have positive alpha to guarantee lambda1 >= lambda2
      DO i=1,ncomp
         w(i)=par(3*i-2)
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
      END DO
C  precompute common terms
      blam=b*lambda
      f1=w0*exp(-blam*(alpha+1.d0))
      fval = f1
      DO i=1,ncomp
C
C   dgt(1,i) contains   (g^T n(phi,eta))^2
C   dgt(2,i) contains  d(g^T n(phi,eta))^2 / d phi
C   dgt(2,i) contains  d(g^T n(phi,eta))^2 / d eta
C
         call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
         adgtd(i)=1.d0+alpha*dgtd(1,i)
         embclgd2(i)=exp(-blam*adgtd(i))
         fval=fval+w(i)*embclgd2(i)
      END DO
C  computed function value
C  now compute derivatives
      f1l=-b*(1.d0+alpha)*f1
      f2l=0.d0
      f2c=0.d0
      DO i=1,ncomp
         f2l=f2l+w(i)*embclgd2(i)*adgtd(i)
         f2c=f2c+w(i)*embclgd2(i)*dgtd(1,i)
         dval(3*i-1)=-alpha*blam*w(i)*embclgd2(i)*dgtd(2,i)
C  derivative with respect to phi_i
         dval(3*i)=-alpha*blam*w(i)*embclgd2(i)*dgtd(3,i)
C  derivative with respect to eta_i
         dval(3*i-2)=embclgd2(i)
C  derivative with respect to wi_i
      END DO
      dlam=f1l-b*f2l
C  derivative with respect to lambda
      dalpha=-blam*(f1+f2c)
C  derivative with respect to alpha
      dw0=exp(-blam*(alpha+1.d0))
C  derivative with respect to w0
      RETURN
      END
      subroutine dfmb1(par,npar,w0,lambda,alpha,g,b,fval,dval,dw0,
     1                     dlam)
C
C  compute partial derivatives with respect to w_i, phi_i, eta_i and lambda
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda
C  fval contains function value
C
      implicit none
      integer npar
      double precision par(npar),g(3),dval(npar),alpha,lambda,b,dlam,
     1       w0,dw0
      integer ncomp,i
      double precision w(5),phi(5),eta(5),blam
      double precision embclgd2(5),dgtd(3,5),fval,f1,f1l,f2l,adgtd(5)
C  extract parameters
      ncomp=(npar-1)/3
C need to have positive alpha to guarantee lambda1 >= lambda2
      DO i=1,ncomp
         w(i)=par(3*i-2)
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
      END DO
C  precompute common terms
      blam=b*lambda
      f1=w0*exp(-blam*(alpha+1.d0))
      fval = f1
      DO i=1,ncomp
C
C   dgt(1,i) contains   (g^T n(phi,eta))^2
C   dgt(2,i) contains  d(g^T n(phi,eta))^2 / d phi
C   dgt(2,i) contains  d(g^T n(phi,eta))^2 / d eta
C
         call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
         adgtd(i)=1.d0+alpha*dgtd(1,i)
         embclgd2(i)=exp(-blam*adgtd(i))
         fval=fval+w(i)*embclgd2(i)
      END DO
C  computed function value
C  now compute derivatives
      f1l=-b*(1.d0+alpha)*f1
      f2l=0.d0
      DO i=1,ncomp
         f2l=f2l+w(i)*embclgd2(i)*adgtd(i)
         dval(3*i-1)=-alpha*blam*w(i)*embclgd2(i)*dgtd(2,i)
C  derivative with respect to phi_i
         dval(3*i)=-alpha*blam*w(i)*embclgd2(i)*dgtd(3,i)
C  derivative with respect to eta_i
         dval(3*i-2)=embclgd2(i)
C  derivative with respect to w_i
      END DO
      dlam=f1l-b*f2l
C  derivative with respect to lambda
      dw0=exp(-blam*(alpha+1.d0))
C  derivative with respect to w0
      RETURN
      END
      subroutine dfmb0(par,npar,w0,lambda,alpha,g,b,fval,dval,dw0)
C
C  compute partial derivatives with respect to w_i, phi_i, eta_i and lambda
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda
C  fval contains function value
C
      implicit none
      integer npar
      double precision par(npar),g(3),dval(npar),alpha,lambda,b,dw0,w0
      integer ncomp,i
      double precision w(5),phi(5),eta(5),blam
      double precision embclgd2(5),dgtd(3,5),fval,adgtd(5)
C  extract parameters
      ncomp=npar/3
C need to have positive alpha to guarantee lambda1 >= lambda2
      DO i=1,ncomp
         w(i)=par(3*i-2)
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
      END DO
C  precompute common terms
      blam=b*lambda
      fval = w0*exp(-blam*(alpha+1.d0))
      DO i=1,ncomp
C
C   dgt(1,i) contains   (g^T n(phi,eta))^2
C   dgt(2,i) contains  d(g^T n(phi,eta))^2 / d phi
C   dgt(2,i) contains  d(g^T n(phi,eta))^2 / d eta
C
         call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
         adgtd(i)=1.d0+alpha*dgtd(1,i)
         embclgd2(i)=exp(-blam*adgtd(i))
         fval=fval+w(i)*embclgd2(i)
      END DO
C  computed function value
C  now compute derivatives
      DO i=1,ncomp
         dval(3*i-1)=-alpha*blam*w(i)*embclgd2(i)*dgtd(2,i)
C  derivative with respect to phi_i
         dval(3*i)=-alpha*blam*w(i)*embclgd2(i)*dgtd(3,i)
C  derivative with respect to eta_i
         dval(3*i-2)=embclgd2(i)
C  derivative with respect to eta_i
      END DO
      dw0=exp(-blam*(alpha+1.d0))
C  derivative with respect to w0
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
      implicit none
      integer nvox,ngrad,ns,siind(ns,nvox),m,ntry,nv,
     1       isample(m,ntry)
      double precision si(ngrad,nvox),sms(ngrad),dgrad(ngrad,nv),
     1       egrad(ngrad,nv),z(ngrad,ns),mval(nvox),
     2       vsi(nvox),bv(ngrad),alpha,lambda,z0(ngrad)
      integer i,k,ibest,mode,ind(10),l,ii,iw,wind(6),nwi(6)
      double precision w(1000),krit,work1(1000),work2(12),erg,msi,m2si,
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
          IF(mode.le.1) THEN
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
C
C
C __________________________________________________________________
C
      subroutine getsiibv(si,ngrad,nvox,m,dgrad,bv,nv,alpha,
     1        lambda,egrad,isample,ntry,sms,z0,z,siind,wi,mval,ns)
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
      integer nvox,ngrad,ns,siind(ns,nvox),m,ntry,nv,
     1       isample(m,ntry)
      double precision si(ngrad,nvox),sms(ngrad),dgrad(ngrad,nv),
     1       egrad(ngrad,nv),z(ngrad,ns),mval(nvox),
     2       wi(ns,nvox),bv(ngrad),alpha,lambda,z0(ngrad)
      integer i,k,ibest,mode,ind(10),l,ii,iw,wind(6),nwi(6)
      double precision w(1000),krit,work1(1000),work2(12),erg,msi,m2si,
     1       z1,dng,albv,lbv,w0
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
         krit=1d20
         DO k=1,ntry
            call dcopy(ngrad,si(1,i),1,sms,1)
            call dcopy(ngrad,z0,1,z(1,1),1)
            DO l=1,m
               call dcopy(ngrad,egrad(1,isample(l,k)),1,z(1,l+1),1)
            END DO
         call nnls(z,ngrad,ngrad,m+1,sms,w,erg,work2,work1,ind,mode)
C
C  thats the Hansson-Larsson code for linear regression with positivity
C          constraints
C
            IF(mode.le.1) THEN
               IF(erg.lt.krit) THEN
                 krit=erg
                 ibest=k
                 iw=0
                 w0 = w(1)
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
              wi(1,i)=w0
              DO l=1,iw
                 siind(l+1,i)=isample(wind(l),ibest)
                 wi(l+1,i)=w(wind(l))
              END DO
           END IF
           IF (iw.lt.m) THEN
              DO l=1,m-iw
                 siind(m-l+2,i)=isample(nwi(l),ibest)
                 wi(m-l+2,i)=1e-2
C use something small as initial value to be in the interior
              END DO
           END IF
           mval(i)=krit
        END IF
      END DO
      RETURN
      END
