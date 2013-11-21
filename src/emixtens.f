C
C  version for use with L-BFGS-B
C
      subroutine erskmxl2(par,npar,siq,g,b,ng,fval,vw,nvox,risk)
C
C    tensor-mixture models with lambda2 and fa as parameters
C     compute sum_j sum_i vw(j)(siq_j-f_j(par,g(,i),b(i))^2  #j over neighboring voxel
      implicit logical(a-z)
      integer npar,ng,nvox
      real*8 par(npar),g(3,ng),siq(ng,nvox),b(ng),lambda,alpha,risk,
     1       fval(nvox),vw(nvox)
      integer i,ncomp,j
      real*8 resi
      lambda=par(1)
      alpha=par(2)
      ncomp = (npar-2)/(nvox+2)
      risk=0.d0
      DO i=1,ng
         call efmixl(par(3),3*ncomp,npar-2,lambda,alpha,g(1,i),b(i),
     1               fval)
         Do j=1,nvox
            resi=siq(i,j)-fval(j)
            risk=risk+vw(j)*resi*resi
         END DO
      END DO
      RETURN
      END
C
C    Evaluations for tensor-mixture models
C
      subroutine efmixl(par,npar0,npar,lambda,alpha,g,b,fval)
C
C  compute expected value
C  par(npar0) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc
C  par((npar0+1):npar) contains weights w_{ij} for neighboring voxel j
C  fval has length 1+(npar-npar0)/ncomp and returns function values at 
C  voxel j gradient direction g
C
      implicit logical (a-z)
      integer npar,npar0
      real*8 par(npar),g(3),fval(*),lambda,alpha
      integer ncomp,i,nvox,j,nstart
      real*8 w(5),phi(5),eta(5),b,wi,sw,blam,f1,f10,expi(5)
      real*8 ddot3sq
      external ddot3sq
C  extract parameters
      ncomp=npar0/3
      sw=1.d0
      DO i=1,ncomp
         wi=par(3*i-2)
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
      blam=b*lambda
      f10=exp(-blam*(alpha+1.d0))
      f1=f10
      DO i=1,ncomp
         expi(i)=exp(-blam*(1.d0+alpha*ddot3sq(phi(i),eta(i),g)))
         f1=f1+w(i)*expi(i)
      END DO
      fval(1)=f1/sw
      nvox = (npar-npar0)/ncomp+1
      nstart = npar0
C  now compute function value for neighboring voxel keeping everyting except
C  weights
      DO j=2,nvox
         sw=1.d0
         DO i=1,ncomp
            wi=par(nstart+i)
            w(i)=wi
            sw=sw+wi
         END DO
         f1=f10
         DO i=1,ncomp
            f1=f1+w(i)*expi(i)
         END DO
         fval(j)=f1/sw
         nstart=nstart+ncomp
      END DO
      RETURN
      END
      subroutine edrskml2(par,npar,ncomp,siq,g,b,ng,fval,vw,nvox,
     1                    dval,drisk0,driskla,drisk,dlam,dalpha)
C
C    tensor-mixture models with lambda2 and fa as parameters
C     compute sum_i (siq-f(par,g(,i),b(i))^2
      implicit logical(a-z)
      integer npar,ng,nvox,ncomp
      real*8 par(npar),g(3,ng),siq(ng,nvox),b(ng),lambda,alpha,
     1       drisk(npar),fval(nvox),vw(nvox),drisk0(3,ncomp,nvox),
     2       driskla(2,nvox),dval(3,ncomp,nvox)
      integer i,j,npar0,k
      real*8 resi,dlam(nvox),dalpha(nvox)
C  calculate numerical gradients for comparison
C      ncomp = (npar-2)/(2*nvox+1)
      npar0 = 2+3*ncomp
      lambda=par(1)
      alpha=par(2)
      DO k=1,ncomp
         DO j=1,nvox
            drisk0(1,k,j)=0.d0   
            drisk0(2,k,j)=0.d0   
            drisk0(3,k,j)=0.d0   
         END DO
      END DO
      DO i=1,ng
         call edfml2(par(3),3*ncomp,npar-2,lambda,alpha,g(1,i),b(i),
     1               fval,dval,dlam,dalpha)  
         DO j=1,nvox       
            resi=siq(i,j)-fval(j)
            DO k=1,ncomp
               drisk0(1,k,j)=drisk0(1,k,j)-resi*dval(1,k,j)
               drisk0(2,k,j)=drisk0(2,k,j)-resi*dval(2,k,j)
               drisk0(3,k,j)=drisk0(3,k,j)-resi*dval(3,k,j)
            END DO
            driskla(1,j)=driskla(1,j)-resi*dlam(j)
            driskla(2,j)=driskla(2,j)-resi*dalpha(j)
         END DO
      END DO
      DO k=1,npar
         drisk(k) = 0.d0
      END DO
      DO j=1,nvox
         drisk(1) = drisk(1)+vw(j)*driskla(1,j)
         drisk(2) = drisk(2)+vw(j)*driskla(2,j)
         DO k=1,ncomp
            if(j.eq.1) drisk(3*k)=drisk(3*k)+vw(j)*drisk0(1,k,j)
            drisk(3*k+1)=drisk(3*k+1)+vw(j)*drisk0(2,k,j)
            drisk(3*k+2)=drisk(3*k+2)+vw(j)*drisk0(3,k,j)
            if(j.gt.1) THEN
               drisk(npar0+(j-2)*ncomp+k)=drisk(npar0+(j-2)*ncomp+k)+
     1              vw(j)*drisk0(1,k,j)
            END IF
         END DO
      END DO
      DO k=1,npar
         drisk(k)=2.d0*drisk(k) 
      END DO
      RETURN
      END
      subroutine edfml2(par,npar0,npar,lambda,alpha,g,b,fval,dval,dlam,
     1                  dalpha)
C
C  compute partial derivatives with respect to w_i, phi_i, eta_i, lambda and alpha
C  par(npar) containes  w_1, phi_1, eta_1, ..., w_nc, phi_nc, eta_nc, lambda, alpha
C  fval contains function value
      implicit logical (a-z)
      integer npar,npar0
      real*8 par(npar),g(3),dval(npar0,*),fval(*),alpha,lambda,b,
     1       dlam(*),dalpha(*)
      integer ncomp,i,nstart,nvox,j
      real*8 w(5),phi(5),eta(5),wi,sw,blam
      real*8 embclgd2(5),dgtd(3,5),f1,f10,f1l,f2l,f2c,adgtd(5)
C  extract parameters
      ncomp=npar0/3
C need to have positive alpha to guarantee lambda1 >= lambda2
      sw=1.d0
      DO i=1,ncomp
         wi=par(3*i-2)
         w(i)=wi
         phi(i)=par(3*i-1)
         eta(i)=par(3*i)
         sw=sw+wi
      END DO
C  precompute common terms
      blam=b*lambda
      f10=exp(-blam*(alpha+1.d0))
      f1 = f10
      DO i=1,ncomp
         call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
         adgtd(i)=1.d0+alpha*dgtd(1,i)
         embclgd2(i)=exp(-blam*adgtd(i))
         f1=f1+w(i)*embclgd2(i)
      END DO
C  compute function value
      fval(1)=f1/sw
C  compute derivatives
      f1l=-b*(1.d0+alpha)*f10
      f2l=0.d0
      f2c=0.d0
      DO i=1,ncomp
         f2l=f2l+w(i)*embclgd2(i)*adgtd(i)
         f2c=f2c+w(i)*embclgd2(i)*dgtd(1,i)
         dval(3*i-1,1)=-alpha*blam*w(i)*embclgd2(i)*dgtd(2,i)/sw
C  derivative with respect to phi_i
         dval(3*i,1)=-alpha*blam*w(i)*embclgd2(i)*dgtd(3,i)/sw
C  derivative with respect to eta_i
         dval(3*i-2,1)=(embclgd2(i)-fval(1))/sw
C  derivative with respect to w_i1
      END DO
      dlam(1)=(f1l-b*f2l)/sw
C  derivative with respect to lambda
      dalpha(1)=-blam/sw*(f10+f2c)
      nvox = (npar-npar0)/ncomp+1
      nstart = npar0
      DO j=2,nvox
         sw=1.d0
         DO i=1,ncomp
            wi=par(nstart+i)
            w(i)=wi
            sw=sw+wi
         END DO
         f1 = f10
         DO i=1,ncomp
            call dgtddphi(phi(i),eta(i),g, dgtd(1,i))
            adgtd(i)=1.d0+alpha*dgtd(1,i)
            embclgd2(i)=exp(-blam*adgtd(i))
            f1=f1+w(i)*embclgd2(i)
         END DO
C  compute function value
         fval(j)=f1/sw
         f2l=0.d0
         f2c=0.d0
         DO i=1,ncomp
            f2l=f2l+w(i)*embclgd2(i)*adgtd(i)
            f2c=f2c+w(i)*embclgd2(i)*dgtd(1,i)
            dval(3*i-1,j)=-alpha*blam*w(i)*embclgd2(i)*dgtd(2,i)/sw
C  derivative with respect to phi_i
            dval(3*i,j)=-alpha*blam*w(i)*embclgd2(i)*dgtd(3,i)/sw
C  derivative with respect to eta_i
            dval(3*i-2,j)=(embclgd2(i)-fval(j))/sw
C  derivative with respect to w_ij
         END DO
         dlam(j)=(f1l-b*f2l)/sw
C  derivative with respect to lambda
         dalpha(j)=-blam/sw*(f10+f2c)  
         nstart=nstart+ncomp
      END DO
C  derivative with respect to c
      RETURN
      END
C
C    construct index file for neighborhoods
C    in mask(n1,n2,n3), vn(3,nvox) - x,y,z coordinates 
C    out  index(nvox,n)  indices for neighboring voxel
C                        contains 0 if voxel outside mask
C
      subroutine nvindex(n1,n2,n3,mask,nvox,vn,ind,n)
      implicit logical (a-z)
      integer n1,n2,n3,nvox,n,ind(nvox,n),vn(3,nvox)
      logical mask(n1,n2,n3),valid
      integer i1,i2,i3,i,j1,j2,j3,im,j
      i=0
      im=0 
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               i=i+1
               if(mask(i1,i2,i3)) THEN
                  im=im+1
                  ind(1,im) = i
                  DO j=2,nvox
                     j1=i1+vn(1,j)
                     j2=i2+vn(2,j)
                     j3=i3+vn(3,j)
                     valid=j1.ge.1.and.j1.le.n1
                     valid=valid.and.j2.ge.1.and.j2.le.n2
                     valid=valid.and.j3.ge.1.and.j3.le.n3
                     if(valid) valid=valid.and.mask(j1,j2,j3)
                     if(valid) THEN
                        ind(j,im)=j1+(j2-1)*n1+(j3-1)*n1*n2
                     ELSE
                        ind(j,im)=0
                     ENDIF
                 END DO
               END IF
            END DO
         END DO
      END DO
      n=im
      RETURN
      END
      
      
      
       
