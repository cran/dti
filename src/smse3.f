      subroutine ghfse3i(i4,kstar,k456,ng,
     1                    kappa,vext,h,varred,n,dist)
C
C   compute bandwidth sequence for given kappa and gradients  
C   lkfse3i0 computes weighting schemes and corresponding variance reduction
C   k456,ng (input) contain auxiliary statistics(gradients)
C
      implicit logical (a-z)
      integer ng,i4,kstar,n,dist
      real*8 k456(3,ng,ng),vext(2),
     1       kappa,h(kstar),varred(kstar)
      logical getkappa
      integer k,n0,maxn
      real*8 hakt,hakt0,vr,ch,chk,vred,v0r
      ch=1.25d0
      getkappa=.FALSE.
      hakt=1.d0
C   initialize kappa
C   loop over steps
      call lkfse3i0(hakt,kappa/hakt,i4,k456,ng,
     1                   vext,vr,n,dist)
      chk=ch*vr
      maxn = 1
      DO k=1,kstar
         call lkfse3i0(hakt,kappa/hakt,i4,k456,ng,
     1                   vext,vr,n,dist)
C  search for new hakt   
         vred=vr/chk
         DO WHILE (vred.lt.1) 
            hakt=hakt*1.05
            call lkfse3i0(hakt,kappa/hakt,i4,k456,ng,
     1                   vext,vr,n,dist)
            vred=vr/chk
         END DO
         DO WHILE (vred.gt.1.01) 
            hakt0=hakt
            v0r=vr
            n0=n
            hakt=hakt/1.005
            call lkfse3i0(hakt,kappa/hakt,i4,k456,ng,
     1                   vext,vr,n,dist)
            vred=vr/chk
            If (vred.lt.1) THEN
               hakt=hakt0
               vr=v0r
               n=n0
            END IF
         END DO
         h(k) = hakt
         varred(k) = vr
         chk=chk*ch
         maxn=max(maxn,n)
         if(k.eq.kstar) THEN
            call lkfse3i0(h(k),kappa/hakt,i4,k456,ng,
     1                   vext,vr,n,dist)
         END IF
C  number of positive weights for last step in n
      END DO
      n=maxn
      RETURN
      END
      subroutine exppm6(p,ex)
      implicit logical (a-z)
      real*8 p,ex(3,3)
      ex(1,1)=dcos(p)
      ex(1,2)=dsin(p)
      ex(1,3)=0.d0
      ex(2,1)=-dsin(p)
      ex(2,2)=dcos(p)
      ex(2,3)=0.d0
      ex(3,1)=0.d0
      ex(3,2)=0.d0
      ex(3,3)=1.d0
      RETURN
      END
      subroutine exppm5(p,ex)
      implicit logical (a-z)
      real*8 p,ex(3,3)
      ex(1,1)=dcos(p)
      ex(1,2)=0.d0
      ex(1,3)=-dsin(p)
      ex(2,1)=0.d0
      ex(2,2)=1.d0
      ex(2,3)=0.d0
      ex(3,1)=dsin(p)
      ex(3,2)=0.d0
      ex(3,3)=dcos(p)
      RETURN
      END
      subroutine exppm4(p,b,ex)
      implicit logical (a-z)
      real*8 b,p,ex(3,3)
      real*8 D,sqr2,sb,pDs,D2,cpds,spds,spdsh
      sqr2 = dsqrt(2.d0)
      sb = dsin(b)
      D = dsqrt(3.d0-dcos(2.d0*b))
      pDs = p*D/sqr2
      cpds = dcos(pDs)
      spds = dsin(pDs)
      spdsh = dsin(pDs/2.d0)
      D2 = D*D
      ex(1,1)=2.d0*(1.d0+cpDs*sb*sb)/D2
      ex(1,2)=-sqr2/D*sb*spDs
      ex(1,3)=-4.d0/D2*sb*spdsh*spdsh
      ex(2,1)=-ex(1,2)
      ex(2,2)=cpDs
      ex(2,3)=sqr2/D*spDs
      ex(3,1)=ex(1,3)
      ex(3,2)=-ex(2,3)
      ex(3,3)=1.d0-2.d0/D2*(1.d0-cpDs)
      RETURN
      END
      subroutine k456krb(par,b,matm,erg)
C
C  Solve exponential equation for dicrepance parameters
C  compute ||\prod_{i=4}^6 exp(par[i] m_i) - matm||^2
C
      implicit logical (a-z)
      real*8 par(3),b,matm(3,3),erg
      integer i1,i2
      real*8 s,z,em4(3,3),em5(3,3),em6(3,3),am4(3,3),am5(3,3)
      call exppm4(par(1),b,em4)
      call exppm5(par(2),em5)
      call exppm6(par(3),em6)
      DO i1=1,3
         DO i2=1,3
            am5(i1,i2)=em5(i1,1)*em6(1,i2)+em5(i1,2)*em6(2,i2)+
     1                                     em5(i1,3)*em6(3,i2)
            END DO
         END DO
      DO i1=1,3
         DO i2=1,3
            am4(i1,i2)=em4(i1,1)*am5(1,i2)+em4(i1,2)*am5(2,i2)+
     1                                     em4(i1,3)*am5(3,i2)
         END DO
      END DO
      s=0.d0
      DO i1=1,3
         DO i2=1,3
            z=matm(i1,i2)-am4(i1,i2)
            s=s+z*z
         END DO
      END DO
      erg=s
      RETURN
      END
      subroutine abofg(g,n,bg)
C
C  compute spherical coordinates for gradient vectors
C
      implicit logical (a-z)
      integer n
      real*8 g(3,n),bg(2,n)
      integer i
      real*8 z,g1,g2,g3,beta,gamma,onemeps
      onemeps=1.d0-1d-10
      DO i=1,n
         g1=g(1,i)
         g2=g(2,i)
         g3=g(3,i)
C standardize to length 1
         z=g1*g1+g2*g2+g3*g3
         z=sqrt(z)
         g1=g1/z
         g2=g2/z
         g3=g3/z
         beta=asin(g1)
         if(abs(g1).lt.onemeps) THEN
            z=g3/cos(beta)
            if(abs(z).lt.onemeps) THEN
               gamma=acos(z)
            ELSE
               gamma=1.570796327d0-sign(1.570796327d0,z)
            END IF
         ELSE 
            gamma=0.d0
         END IF
         if(g2.gt.0.d0) gamma = -gamma
         bg(1,i)=beta
         bg(2,i)=gamma
      END DO
      RETURN
      END
      subroutine bgstats(g,n,bg,bghat)
      implicit logical (a-z)
      integer n
      real*8 g(3,n),bg(2,n),bghat(2,n,n)
      integer i1,i2
      real*8 dgamma,cb1,sb1,cb2,sb2,betah,cbh,z,gammah,cdg
C   first get sperical coordinates in bg
      call abofg(g,n,bg)
      DO i1=1,n
         sb1=sin(bg(1,i1))
         cb1=cos(bg(1,i1))
C    die normalen-vektoren der Gradienten
         DO i2=1,n
            dgamma=bg(2,i1)-bg(2,i2)
            cdg=cos(dgamma)
            if(abs(cdg).gt.1-1d-8) THEN
               bghat(1,i1,i2)=asin(sin(bg(1,i1)-
     1                             sign(1.d0,cdg)*bg(1,i2)))
               bghat(2,i1,i2)=0.d0
               CYCLE
            END IF
            sb2=sin(bg(1,i2))
            cb2=cos(bg(1,i2))
            z=sb1*cb2-cb1*sb2*cdg
            betah=asin(z)
            cbh=cos(betah)
            IF(abs(cbh).gt.1d-8) THEN
               z=cb1*sin(dgamma)/cbh
               z=sign(min(1.d0,abs(z)),z)
               gammah=asin(z)
            ELSE
               if(abs(cb1).gt.1d-6) THEN
C   this should not happen
                  call dblepr("cb1",3,cb1,1)
                  call dblepr("cbh",3,cbh,1)
               END IF
               gammah = dgamma*sign(1d0,cb1*cbh)
            END IF
            bghat(1,i1,i2)=betah
            bghat(2,i1,i2)=gammah
C    die spherischen Koordinaten der Gradientenpaare (Parameter der Rotationsmatix)
         END DO
      END DO  
      RETURN
      END
      subroutine lkfse3i(h,kappa,i4,k456,ng,
     1                   vext,ind,wght,n,dist)
      implicit logical (a-z)
      integer ng,n,ind(5,n),i4,dist
      real*8 h,kappa,k456(3,ng,ng),vext(2),wght(n)
      integer ih1,ih2,ih3,i,j1,j2,j3,j4
      real*8 h2,kap2,x1,x2,x3,k4,k5,k6,z,z1,vd2,vd3
      ih1 = max(1.d0,5.0d0*h)
      ih2 = max(1.d0,5.0d0*h/vext(1))
      ih3 = max(1.d0,5.0d0*h/vext(2))
      h2 = h*h
      kap2 = kappa*kappa
      vd2 = vext(1)*vext(1)
      vd3 = vext(2)*vext(2)
      i = 1
      z = 0.d0
      k5 = 0.d0
      k6 = 0.d0
C  just to prevent compiler warnings
      DO j4 = 1,ng
         k4 = k456(1,i4,j4)
         IF(dist.lt.3) THEN
            k5 = k456(2,i4,j4)
            k6 = k456(3,i4,j4)
         END IF
         SELECT CASE (dist)
            CASE (1) 
               z = (k4*k4+k5*k5+abs(k6))/kap2
            CASE (2) 
               z = (k4*k4+k5*k5+k6*k6)/kap2
            CASE (3) 
               z = k4*k4/kap2
            CASE DEFAULT 
               z = k4*k4/kap2
         END SELECT
         if(z.gt.h2) CYCLE
C   last three komponents already to large
         DO j1 = 0,ih1
            x1 = j1
            x1 = z + x1*x1
            if(x1.gt.h2) CYCLE
            DO j2 = -ih2,ih2
               x2 = j2
               x2 = x1 + vd2*x2*x2
               if(x2.gt.h2) CYCLE
               DO j3 = -ih3,ih3
                  x3 = j3
                  z1 = x2+vd3*x3*x3
                  if(z1.gt.h2) CYCLE
                  if(i.gt.n) THEN
                     call intpr("Exceeded max i",14,i,1)
                     call intpr("for i4",6,i4,1)
                     n = i-1
                     return
                  END IF
                  wght(i)= (1.d0-z1/h2)
                  ind(1,i) = j1
                  ind(2,i) = j2
                  ind(3,i) = j3
                  ind(4,i) = i4
                  ind(5,i) = j4
                  i = i+1
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      n = i-1
      RETURN
      END
      subroutine lkfse3i0(h,kappa,i4,k456,ng,vext,vred,n,dist)
      implicit logical (a-z)
      integer ng,n,i4,dist
      real*8 h,kappa,k456(3,ng,ng),vext(2),vred
      integer ih1,ih2,ih3,j1,j2,j3,j4,mj1,mj2,mj3,rad
      real*8 x1,x2,x3,k4,k5,k6,z,z1,
     1       sw,sw2,wght,anz,h2,kap2,vd2,vd3
      ih1 = max(1.d0,5.0d0*h)
      ih2 = max(1.d0,5.0d0*h/vext(1))
      ih3 = max(1.d0,5.0d0*h/vext(2))
      sw=0.d0
      sw2=0.d0
      mj1=0
      mj2=0
      mj3=0
      h2 = h*h
      kap2 = kappa*kappa
      vd2 = vext(1)*vext(1)
      vd3 = vext(2)*vext(2)
      n = 0
      z = 0.d0
      k5 = 0.d0
      k6 = 0.d0
C  just to prevent compiler warnings
      DO j4 = 1,ng
         k4 = k456(1,i4,j4)
         IF(dist.lt.3) THEN
            k5 = k456(2,i4,j4)
            k6 = k456(3,i4,j4)
         END IF
         SELECT CASE (dist)
            CASE (1) 
               z = (k4*k4+k5*k5+abs(k6))/kap2
            CASE (2) 
               z = (k4*k4+k5*k5+k6*k6)/kap2
            CASE (3) 
               z = k4*k4/kap2
            CASE DEFAULT 
               z = k4*k4/kap2
         END SELECT
         if(z.gt.h2) CYCLE
C   last three komponents already to large
         DO j1 = 0,ih1
            x1 = j1
            x1 = z + x1*x1
            if(x1.gt.h2) CYCLE
            DO j2 = -ih2,ih2
               x2 = j2
               x2 = x1 + vd2*x2*x2
               if(x2.gt.h2) CYCLE
               DO j3 = -ih3,ih3
                  x3 = j3
                  z1 = x2+vd2*x3*x3
                  if(z1.gt.h2) CYCLE
                  wght= (1.d0-z1/h2)
C   if j1>0  (-j1,-j2,-j3) gets the same weight, so count it twice
                  if(j1.eq.0) THEN
                     anz=1.d0
                  ELSE
                     anz=2.d0
                  ENDIF
                  sw=sw+anz*wght
                  wght=wght*wght
                  sw2=sw2+anz*wght
                  n=n+1
                  mj1=max(j1,mj1)
                  mj2=max(abs(j2),mj2)
                  mj3=max(abs(j3),mj3)
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      vred = sw*sw/sw2
      rad = max(mj1,max(mj2,mj3))
      if(rad.gt.2.d0*h) THEN
         call dblepr("h",1,h,1)
         call intpr("radius",6,rad,1)
      END IF
      RETURN
      END
      subroutine lkfulse3(h,kappa,k456,ng,vext,ind,wght,n,dist)
      implicit logical (a-z)
      integer ng,n,ind(5,n),dist
      real*8 h(ng),kappa(ng),k456(3,ng,ng),vext(2),wght(n)
      integer ns,ni,i
      ns = 0
      DO i = 1,ng
         ni = n-ns
         call lkfse3i(h(i),kappa(i),i,k456,ng,
     1                vext,ind(1,ns+1),wght(ns+1),ni,dist)
         ns = ns+ni
      END DO
      n = ns
      RETURN
      END      
      subroutine adsmse3p(y,th,ni,mask,n1,n2,n3,ngrad,lambda,ncoils,
     1                    ncores,ind,w,n,thn,ldf,sw,swy,model,minlev)
C   model=1 takes noncentral Chi-sq values in y
C   model=0 takes noncentral Chi values in y
C
C   perform adaptive smoothing on SE(3) 
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C
c   model=0  Chi^2-based KL-distance, y ~ Chi, th on same scale, smooth y
c   model=1  Chi^2-based KL-distance, y ~ Chi^2, th on same scale, smooth y
c   model=2  Gauss-based KL-distance, y ~ Chi, th on same scale, smooth y^2
      implicit logical (a-z)
      integer n1,n2,n3,ngrad,n,ind(5,n),ncoils,model,ncores
      logical mask(n1,n2,n3)
      real*8 y(n1,n2,n3,ngrad),thn(n1,n2,n3,ngrad),ni(n1,n2,n3,ngrad),
     1     ldf(n1,n2,n3,ngrad),th(n1,n2,n3,ngrad)
      real*8 lambda,w(n),sw(ngrad,ncores),swy(ngrad,ncores),
     2       lgfi,dgfi,fici,df,minlev
      integer iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr
      real*8 z,thi,nii,thj,ldfi,ldfj,yj
!$      integer omp_get_thread_num 
!$      external omp_get_thread_num
      real*8 kldisnc1
      external kldisnc1
      df=2.d0*ncoils
      nii=1.d0
      thrednr = 1
C just to prevent a compiler warning
C  precompute values of lgamma(corrected df/2) in each voxel
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,n3,ngrad,ncores,mask,y,thn,ni,ldf,th,w,sw,swy,
C$OMP& ind,ncoils,model,df,n,lambda,minlev)
C$OMP& FIRSTPRIVATE(nii)
C$OMP& PRIVATE(iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,z,thi,thj,
C$OMP& ldfi,ldfj,dgfi,fici,lgfi,yj)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         if(model.lt.2) THEN
            i1=mod(iind,n1)
            if(i1.eq.0) i1=n1
            i2=mod((iind-i1)/n1+1,n2)
            if(i2.eq.0) i2=n2
            i3=(iind-i1-(i2-1)*n1)/n1/n2+1
            if(.not.mask(i1,i2,i3)) CYCLE
C  not needed for Gauss approximation
            DO i4=1,ngrad
               thi=th(i1,i2,i3,i4)
               call lgstats(thi,df,model,ldfi)
               ldf(i1,i2,i3,i4)=ldfi
            END DO
         END IF
      END DO
C$OMP END DO 
C$OMP BARRIER
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1         
         if(.not.mask(i1,i2,i3)) CYCLE
         DO i4=1,ngrad
            sw(i4,thrednr)=0.d0
            swy(i4,thrednr)=0.d0
         END DO
         i4=0
         DO i=1,n
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               thi = th(i1,i2,i3,i4)
               IF(model.ne.2) THEN
                  ldfi=ldf(i1,i2,i3,i4)
                  call ncstats0(thi,ldfi,df,model,lgfi,dgfi,fici)
               END IF
               nii = ni(i1,i2,i3,i4)/lambda
           END IF
            j1=i1+ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE          
            j4=ind(5,i)
C adaptation 
            if(lambda.lt.1d10) THEN
               thj=th(j1,j2,j3,j4)
               if(model.ge.2) THEN
                  z=(thi-thj)
                  z=nii*z*z
               ELSE
                  ldfj=ldf(j1,j2,j3,j4)
                  z=nii*kldisnc1(lgfi,dgfi,fici,thj,ldfj,df,model)
               END IF
C  do not adapt on the sphere !!! 
            ELSE
            z=0.d0
            END IF
            if(z.ge.1.d0) CYCLE
            z=w(i)*min(1.d0,2.d0-2.d0*z)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            yj=y(j1,j2,j3,j4)
            if(model.eq.2) yj=yj*yj
            swy(i4,thrednr)=swy(i4,thrednr)+z*yj
         END DO
C  now opposite directions
         DO i=1,n
            if(ind(1,i).eq.0) CYCLE
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               thi = th(i1,i2,i3,i4)
               IF(model.ne.2) THEN
                  ldfi=ldf(i1,i2,i3,i4)
                  call ncstats0(thi,ldfi,df,model,lgfi,dgfi,fici)
               END IF
               nii = ni(i1,i2,i3,i4)/lambda
            END IF
C
C   handle case j1-i1 < 0 which is not contained in ind 
C   using axial symmetry
C
            j1=i1-ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2-ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3-ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE          
            j4=ind(5,i)
            if(lambda.lt.1d10) THEN
               thj=th(j1,j2,j3,j4)
               if(model.ge.2) THEN
                  z=(thi-thj)
                  z=nii*z*z
               ELSE
                  ldfj=ldf(j1,j2,j3,j4)
                  z=nii*kldisnc1(lgfi,dgfi,fici,thj,ldfj,df,model)
               END IF
C  do not adapt on the sphere !!! 
            ELSE
               z=0.d0
            END IF
            if(z.ge.1.d0) CYCLE
            z=w(i)*min(1.d0,2.d0-2.d0*z)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            yj=y(j1,j2,j3,j4)
            if(model.eq.2) yj=yj*yj
            swy(i4,thrednr)=swy(i4,thrednr)+z*yj
         END DO
         DO i4=1,ngrad
            z=swy(i4,thrednr)/sw(i4,thrednr)
            if(model.eq.2) z=sqrt(z)
            thn(i1,i2,i3,i4) = max(z,minlev)
            ni(i1,i2,i3,i4) = sw(i4,thrednr)
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni)
      RETURN
      END
      subroutine adsmse3m(y,th,ni,sthi,mask,ns,n1,n2,n3,ngrad,lambda,
     1             ncoils,ncores,ind,w,n,thn,sw,swy,si,thi)
C   
C  Multi-shell version (differs in dimension of th 
C  KL-distance based on all spheres and Gauss-approximation only
C
C   perform adaptive smoothing on SE(3) 
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C
c   model=2  Gauss-based KL-distance, y ~ Chi, th on same scale, smooth y^2
      implicit logical (a-z)
      integer ns,n1,n2,n3,ngrad,n,ind(5,n),ncoils,ncores
      logical mask(n1,n2,n3)
      real*8 y(n1,n2,n3,ngrad),thn(n1,n2,n3,ngrad),ni(n1,n2,n3,ngrad),
     1     th(ns,n1,n2,n3,ngrad),sthi(ns,n1,n2,n3,ngrad)
      real*8 lambda,w(n),sw(ngrad,ncores),swy(ngrad,ncores)
      integer iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,k
      real*8 sz,z,nii,thj,si(ns,ncores),thi(ns,ncores) 
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      real*8 kldisnc1
      external kldisnc1
      nii=1.d0
      thrednr = 1
C just to prevent a compiler warning
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,n3,ngrad,ncores,mask,y,thn,ni,th,sthi,w,sw,swy,
C$OMP& ind,ncoils,n,lambda,ns,thi,si)
C$OMP& FIRSTPRIVATE(a,b,nii)
C$OMP& PRIVATE(iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,z,thj,sz)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
!$       thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1 
         if(.not.mask(i1,i2,i3)) CYCLE
         DO i4=1,ngrad
            sw(i4,thrednr)=0.d0
            swy(i4,thrednr)=0.d0
         END DO
         i4=0
         DO i=1,n
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               DO k=1,ns
                  thi(k,thrednr) = th(k,i1,i2,i3,i4)
C   thast the approximated standard deviation
                  si(k,thrednr) = sthi(k,i1,i2,i3,i4)
               END DO
               nii = ni(i1,i2,i3,i4)/lambda
            END IF
            j1=i1+ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE          
            j4=ind(5,i)
C adaptation 
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k,thrednr)-th(k,j1,j2,j3,j4))/si(k,thrednr)
                  sz=sz+z*z
               END DO
               z=nii*sz
C  do not adapt on the sphere !!! 
            ELSE
               z=0.d0
            END IF
            if(z.ge.1.d0) CYCLE
            z=w(i)*min(1.d0,2.d0-2.d0*z)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            swy(i4,thrednr)=swy(i4,thrednr)+z*y(j1,j2,j3,j4)
         END DO
C  now opposite directions
         DO i=1,n
            if(ind(1,i).eq.0) CYCLE
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               DO k=1,ns
                  thi(k,thrednr) = th(k,i1,i2,i3,i4)
                  si(k,thrednr) = sthi(k,i1,i2,i3,i4)
               END DO
               nii = ni(i1,i2,i3,i4)/lambda
            END IF
C
C   handle case j1-i1 < 0 which is not contained in ind 
C   using axial symmetry
C
            j1=i1-ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2-ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3-ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE          
            j4=ind(5,i)
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k,thrednr)-th(k,j1,j2,j3,j4))/si(k,thrednr)
                  sz=sz+z*z
               END DO
               z=nii*sz
C  do not adapt on the sphere !!! 
            ELSE
               z=0.d0
            END IF
            if(z.ge.1.d0) CYCLE
            z=w(i)*min(1.d0,2.d0-2.d0*z)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            swy(i4,thrednr)=swy(i4,thrednr)+z*y(j1,j2,j3,j4)
         END DO
         DO i4=1,ngrad
            thn(i1,i2,i3,i4) = swy(i4,thrednr)/sw(i4,thrednr)
            ni(i1,i2,i3,i4) = sw(i4,thrednr)
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni)
      RETURN
      END
      subroutine adsmse3c(y,y0,th,ni,th0,ni0,mask,ns,n1,n2,n3,ngrad,
     1                lambda,ws0,ncores,ind,w,n,ind0,w0,
     2                n0,thn,nin,th0n,ni0n,sw,swy,thi,nii)
C   
C  Multi-shell version (differs in dimension of th 
C  KL-distance based on all spheres and Gauss-approximation only
C
C   perform adaptive smoothing on SE(3) multishell including s0
C   y  -  si images
C   y0 -  mean s0 image
C   th -  estimated/interpolated \E si on all shells (including 0 shell)
C   ni -  corresponding sum of weights
C   th0 -  estimated/interpolated \E s0 and mean_g(\E si) on all other shells  
C   ni0 -  corresponding sum of weights
C   mask - head mask
C   ns   - number of shells (including 0 shell)
C   n1,n2,n3,ngrad - dimensions, number of gradients (bv!=0)
C   lambda - skale parameter
C   ws0  - relative weight for information from s0 images (should be in [0,1])
C   ncoils - df/2 of \chi distributions
C   minlev - expectation of central chi distr. (needed in variance estimates)
C   ncores - number of cores
C   ind    - index vectors for si weighting schemes 
C   w    - corresponding weights
C   n    - number of si weights
C   ind0    - index vectors for s0 weighting schemes 
C   w0    - corresponding weights
C   n0    - number of s0 weights
C   thn   - new estimates of \E si
C   thn0   - new estimates of \E s0
C   nin   - new sum of weight for si
C   ni0n   - new sum of weight for s0
C...sw,swy,si,thi,nii - working areas
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i) for si images
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C
      implicit logical (a-z)
      integer ns,n1,n2,n3,ngrad,n,n0,ind(5,n),ind0(3,n0),ncores
      logical mask(n1,n2,n3)
      real*8 y(n1,n2,n3,ngrad),y0(n1,n2,n3),th(ns,n1,n2,n3,ngrad),
     1     ni(ns,n1,n2,n3,ngrad),th0(ns,n1,n2,n3),ni0(ns,n1,n2,n3),
     2     thn(n1,n2,n3,ngrad),th0n(n1,n2,n3),
     3     nin(n1,n2,n3,ngrad),ni0n(n1,n2,n3)
      real*8 w(n),w0(n0),lambda,thi(ns,ncores),ws0,
     4     sw(ngrad,ncores),swy(ngrad,ncores),nii(ns,ncores)
      integer iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,k
      real*8 sz,z,sw0,swy0
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      thrednr = 1

C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ns,n1,n2,n3,ngrad,n,n0,ind,ind0,ncoils,ncores,y,y0,
C$OMP&       th,ni,th0,ni0,w,w0,thn,th0n,nin,ni0n,thi,sw,swy,nii,
C$OMP&       lambda,mask,ws0)
C$OMP& PRIVATE(iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,k,sz,z,
C$OMP&       sw0,swy0)
C$OMP DO SCHEDULE(GUIDED)
C  First si - images
      DO iind=1,n1*n2*n3
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1         
         if(.not.mask(i1,i2,i3)) CYCLE
         DO i4=1,ngrad
            sw(i4,thrednr)=0.d0
            swy(i4,thrednr)=0.d0
         END DO
         i4=0
         DO i=1,n
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               DO k=1,ns
                  thi(k,thrednr) = th(k,i1,i2,i3,i4)
                  nii(k,thrednr) = ni(k,i1,i2,i3,i4)/lambda
               END DO
               nii(1,thrednr)=ws0*nii(1,thrednr)
            END IF
            j1=i1+ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE          
            j4=ind(5,i)
C adaptation 
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k,thrednr)-th(k,j1,j2,j3,j4))
                  sz=sz+nii(k,thrednr)*z*z
               END DO
C  do not adapt on the sphere !!! 
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w(i)*min(1.d0,2.d0-2.d0*sz)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            swy(i4,thrednr)=swy(i4,thrednr)+z*y(j1,j2,j3,j4)
         END DO
C  now opposite directions
         DO i=1,n
            if(ind(1,i).eq.0) CYCLE
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               DO k=1,ns
                  thi(k,thrednr) = th(k,i1,i2,i3,i4)
                  nii(k,thrednr) = ni(k,i1,i2,i3,i4)/lambda
               END DO
               nii(1,thrednr)=ws0*nii(1,thrednr)
C  first component corresponds to S0 image, ws0 is used to downweight its influence
C  when smoothing diffusion weighted data
            END IF
C
C   handle case j1-i1 < 0 which is not contained in ind 
C   using axial symmetry
C
            j1=i1-ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2-ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3-ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE          
            j4=ind(5,i)
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k,thrednr)-th(k,j1,j2,j3,j4))
                  sz=sz+nii(k,thrednr)*z*z
               END DO
C  do not adapt on the sphere !!! 
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w(i)*min(1.d0,2.d0-2.d0*sz)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            swy(i4,thrednr)=swy(i4,thrednr)+z*y(j1,j2,j3,j4)
         END DO
         DO i4=1,ngrad
            thn(i1,i2,i3,i4) = swy(i4,thrednr)/sw(i4,thrednr)
            nin(i1,i2,i3,i4) = sw(i4,thrednr)
         END DO
C    now the s0 image in iind
         sw0=0.d0
         swy0=0.d0
         DO k=1,ns
            thi(k,thrednr) = th0(k,i1,i2,i3)
            nii(k,thrednr) = ni0(k,i1,i2,i3)/lambda
         END DO
         DO i=1,n0
            j1=i1+ind0(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind0(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind0(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE          
C adaptation 
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k,thrednr)-th0(k,j1,j2,j3))
                  sz=sz+nii(k,thrednr)*z*z
               END DO
C  do not adapt on the sphere !!! 
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w0(i)*min(1.d0,2.d0-2.d0*sz)
            sw0=sw0+z
            swy0=swy0+z*y0(j1,j2,j3)
         END DO
C  now opposite directions
         DO i=1,n0
            if(ind0(1,i).eq.0) CYCLE
C
C   handle case j1-i1 < 0 which is not contained in ind 
C   using axial symmetry
C
            j1=i1-ind0(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2-ind0(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3-ind0(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE          
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k,thrednr)-th0(k,j1,j2,j3))
                  sz=sz+nii(k,thrednr)*z*z
               END DO
C  do not adapt on the sphere !!! 
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w0(i)*min(1.d0,2.d0-2.d0*sz)
            sw0=sw0+z
            swy0=swy0+z*y0(j1,j2,j3)
         END DO
         th0n(i1,i2,i3) = swy0/sw0
         ni0n(i1,i2,i3) = sw0
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,nin,th0n,ni0n)
      RETURN
      END
      subroutine adsmse3q(y,y0,th,ni,th0,ni0,mask,ns,nsp1,n1,n2,n3,
     1                ngrad,lambda,ncores,ind,w,n,
     2                ind0,w0,n0,thn,nin,th0n,ni0n,sw,swy,thi,nii)
C   
C  Multi-shell version (differs in dimension of th 
C  KL-distance based on all spheres and Gauss-approximation only
C
C   perform adaptive smoothing on SE(3) multishell including s0
C   y  -  si images
C   y0 -  mean s0 image
C   th -  estimated/interpolated \E si on all shells (including 0 shell)
C   ni -  corresponding sum of weights
C   th0 -  estimated/interpolated \E s0 and mean_g(\E si) on all other shells  
C   ni0 -  corresponding sum of weights
C   mask - head mask
C   ns   - number of shells 
C   nsp1   - number of shells (including 0 shell)
C   n1,n2,n3,ngrad - dimensions, number of gradients (bv!=0)
C   lambda - skale parameter
C   ncores - number of cores
C   ind    - index vectors for si weighting schemes 
C   w    - corresponding weights
C   n    - number of si weights
C   ind0    - index vectors for s0 weighting schemes 
C   w0    - corresponding weights
C   n0    - number of s0 weights
C   thn   - new estimates of \E si
C   thn0   - new estimates of \E s0
C   nin   - new sum of weight for si
C   ni0n   - new sum of weight for s0
C...sw,swy,si,thi,nii - working areas
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i) for si images
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C
      implicit logical (a-z)
      integer ns,nsp1,n1,n2,n3,ngrad,n,n0,ind(5,n),ind0(3,n0),ncores
      logical mask(n1,n2,n3)
      real*8 y(n1,n2,n3,ngrad),y0(n1,n2,n3),th(ns,n1,n2,n3,ngrad),
     1    ni(ns,n1,n2,n3,ngrad),th0(nsp1,n1,n2,n3),ni0(nsp1,n1,n2,n3),
     2    thn(n1,n2,n3,ngrad),th0n(n1,n2,n3),
     3    nin(n1,n2,n3,ngrad),ni0n(n1,n2,n3)
      real*8 w(n),w0(n0),lambda,thi(nsp1,ncores),s0i,
     1      sw(ngrad,ncores),swy(ngrad,ncores),nii(nsp1,ncores)
      integer iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,k
      real*8 sz,z,sw0,swy0
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      thrednr = 1
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ns,nsp1,n1,n2,n3,ngrad,n,n0,ind,ind0,ncores,y,y0,
C$OMP&       th,ni,th0,ni0,w,w0,thn,th0n,nin,ni0n,thi,sw,swy,nii,
C$OMP&       lambda,mask)
C$OMP& PRIVATE(iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,k,sz,z,
C$OMP&       sw0,swy0,s0i)
C$OMP DO SCHEDULE(GUIDED)
C  First si - images
      DO iind=1,n1*n2*n3
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1         
         if(.not.mask(i1,i2,i3)) CYCLE
         DO i4=1,ngrad
            sw(i4,thrednr)=0.d0
            swy(i4,thrednr)=0.d0
         END DO
         i4=0
         s0i = th0(1,i1,i2,i3)
         DO i=1,n
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               DO k=1,ns
                  thi(k,thrednr) = th(k,i1,i2,i3,i4)
                  nii(k,thrednr) = ni(k,i1,i2,i3,i4)/lambda*s0i*s0i
               END DO
            END IF
            j1=i1+ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE          
            j4=ind(5,i)
C adaptation 
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=thi(k,thrednr)-th(k,j1,j2,j3,j4)
                  sz=sz+nii(k,thrednr)*z*z
               END DO
C  do not adapt on the sphere !!! 
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w(i)*min(1.d0,2.d0-2.d0*sz)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            swy(i4,thrednr)=swy(i4,thrednr)+z*y(j1,j2,j3,j4)
         END DO
C  now opposite directions
         DO i=1,n
            if(ind(1,i).eq.0) CYCLE
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               DO k=1,ns
                  thi(k,thrednr) = th(k,i1,i2,i3,i4)
                  nii(k,thrednr) = ni(k,i1,i2,i3,i4)/lambda*s0i*s0i
               END DO
            END IF
C
C   handle case j1-i1 < 0 which is not contained in ind 
C   using axial symmetry
C
            j1=i1-ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2-ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3-ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE          
            j4=ind(5,i)
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=thi(k,thrednr)-th(k,j1,j2,j3,j4)
                  sz=sz+nii(k,thrednr)*z*z
               END DO
C  do not adapt on the sphere !!! 
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w(i)*min(1.d0,2.d0-2.d0*sz)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            swy(i4,thrednr)=swy(i4,thrednr)+z*y(j1,j2,j3,j4)
         END DO
         DO i4=1,ngrad
            thn(i1,i2,i3,i4) = swy(i4,thrednr)/sw(i4,thrednr)
            nin(i1,i2,i3,i4) = sw(i4,thrednr)
         END DO
C    now the s0 image in iind
         sw0=0.d0
         swy0=0.d0
         DO k=1,nsp1
            thi(k,thrednr) = th0(k,i1,i2,i3)
            nii(k,thrednr) = ni0(k,i1,i2,i3)/lambda
         END DO
         DO i=1,n0
            j1=i1+ind0(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind0(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind0(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE          
C adaptation 
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,nsp1
                  z=thi(k,thrednr)-th0(k,j1,j2,j3)
                  sz=sz+nii(k,thrednr)*z*z
               END DO
C  do not adapt on the sphere !!! 
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w0(i)*min(1.d0,2.d0-2.d0*sz)
            sw0=sw0+z
            swy0=swy0+z*y0(j1,j2,j3)
         END DO
C  now opposite directions
         DO i=1,n0
            if(ind0(1,i).eq.0) CYCLE
C
C   handle case j1-i1 < 0 which is not contained in ind 
C   using axial symmetry
C
            j1=i1-ind0(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2-ind0(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3-ind0(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE          
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,nsp1
                  z=thi(k,thrednr)-th0(k,j1,j2,j3)
                  sz=sz+nii(k,thrednr)*z*z
               END DO
C  do not adapt on the sphere !!! 
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w0(i)*min(1.d0,2.d0-2.d0*sz)
            sw0=sw0+z
            swy0=swy0+z*y0(j1,j2,j3)
         END DO
         th0n(i1,i2,i3) = swy0/sw0
         ni0n(i1,i2,i3) = sw0
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,nin,th0n,ni0n)
      RETURN
      END
      subroutine asmse30p(y,th,ni,mask,n1,n2,n3,lambda,ncoils,ind,w,
     1                    n,starts,nstarts,thn,ldf,swi,model,minlev)
C   model=1 takes noncentral Chi-sq values in y0
C   model=0 takes noncentral Chi values in y0
C   perform adaptive smoothing on R^3
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C
c   model=0  Chi^2-based KL-distance, y ~ Chi, th on same scale, smooth y
c   model=1  Chi^2-based KL-distance, y ~ Chi^2, th on same scale, smooth y
c   model=2  Gauss-based KL-distance, y ~ Chi, th on same scale, smooth y^2

      implicit logical (a-z)
      integer n1,n2,n3,n,ind(5,n),starts(1),nstarts,ncoils,model
      logical mask(n1,n2,n3)
      real*8 y(n1,n2,n3),th(n1,n2,n3),ni(1),thn(1),minlev,
     1       lambda,w(n),sw0,swy0,swi(nstarts),ldf(n1,n2,n3)
      integer i,i0,i1,i2,i3,j1,j2,j3,l1,l2,l3,nn
      real*8 z,ldfi,lgfi,dgfi,fici,sc,df,yj,maxswi
      real*8 kldisnc1
      external kldisnc1
      df=2.d0*ncoils
      nn=n1*n2*n3
C
C  compute location weights in swi(i0) as sum of location weights used for same 
C  relative voxel in R3 and coinciding gradients (sum over gradients)
C
      maxswi = 0.d0
      DO i0=1,nstarts
         swi(i0)=0.d0
         DO i=starts(i0)+1,starts(i0+1)
            swi(i0)=swi(i0)+w(i)
         END DO
         swi(i0)=swi(i0)/(starts(i0+1)-starts(i0))
      END DO
      maxswi=1.d0
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,n3,nn,ncoils,mask,y,th,ni,thn,ind,ldf,n,nstarts,
C$OMP& starts,w,swi,df,minlev)
C$OMP& FIRSTPRIVATE(model,lambda,maxswi)
C$OMP& PRIVATE(sw0,swy0,i,i0,i1,i2,i3,j1,j2,j3,l1,l2,l3,z,lgfi,dgfi,
C$OMP& ldfi,fici,sc,yj)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,nn
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1
         thn(i) = th(i1,i2,i3)
         if(.not.mask(i1,i2,i3)) CYCLE
         if(model.lt.2) THEN
            call lgstats(th(i1,i2,i3),df,model,ldfi)
            ldf(i1,i2,i3) = ldfi
         END IF
      END DO
C$OMP END DO
C$OMP BARRIER
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,nn
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1         
         if(.not.mask(i1,i2,i3)) CYCLE
         sw0=0.d0
         swy0=0.d0
         sc = lambda/ni(i)
         IF(model.lt.2) THEN
            call ncstats0(th(i1,i2,i3),ldf(i1,i2,i3),df,
     1                       model,lgfi,dgfi,fici)
         END IF
         DO i0=1,nstarts
            l1=ind(1,1+starts(i0))
            j1=i1+l1
            if(j1.le.0.or.j1.gt.n1) CYCLE
            l2=ind(2,1+starts(i0))
            j2=i2+l2
            if(j2.le.0.or.j2.gt.n2) CYCLE
            l3=ind(3,1+starts(i0))
            j3=i3+l3
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE
            if(model.ge.2) THEN
               z=th(i1,i2,i3)-th(j1,j2,j3)
               z=z*z
            ELSE
               z=kldisnc1(lgfi,dgfi,fici,th(j1,j2,j3),
     1                        ldf(j1,j2,j3),df,model)
            END IF
            if(z.ge.sc) CYCLE
            if(z.lt.sc) THEN
               z=swi(i0)*min(1.d0,2.d0-2.d0*z/sc)
               sw0=sw0+z
               yj=y(j1,j2,j3)
               if(model.eq.2) yj=yj*yj
               swy0=swy0+z*yj
            END IF
         END DO
         DO i0=1,nstarts
            l1=ind(1,1+starts(i0))
            if(l1.eq.0) CYCLE
C  this part for negative l1 only (opposite directions)
            j1=i1-l1
            if(j1.le.0.or.j1.gt.n1) CYCLE
            l2=ind(2,1+starts(i0))
            j2=i2-l2
            if(j2.le.0.or.j2.gt.n2) CYCLE
            l3=ind(3,1+starts(i0))
            j3=i3-l3
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(.not.mask(j1,j2,j3)) CYCLE
            if(model.ge.2) THEN
               z=th(i1,i2,i3)-th(j1,j2,j3)
               z=z*z
            ELSE
               z=kldisnc1(lgfi,dgfi,fici,th(j1,j2,j3),
     1         ldf(j1,j2,j3),df,model)
            END IF
            if(z.ge.sc) CYCLE
            if(z.lt.sc) THEN
               z=swi(i0)*min(1.d0,2.d0-2.d0*z/sc)
               sw0=sw0+z
               yj=y(j1,j2,j3)
               if(model.eq.2) yj=yj*yj
               swy0=swy0+z*yj
            END IF
         END DO
         z = swy0/sw0
         if(model.eq.2) z=sqrt(z)
         thn(i) = max(z,minlev)
         ni(i) = sw0/maxswi
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni)
      RETURN
      END
C
C   Kullback-leibler distance for noncentral-Chi-distributions
C   with parameters thi and thj and 2*nc degrees of freedom
C
C
C   variant with precomputed quantities
C
      subroutine lgstats(thi,df,model,lgfi)
      implicit logical (a-z)
      integer model
      real*8 thi,lgfi,df
      real*8 mu2i,z1,z2,fi
      real*8 lgammaf,digammaf
      external lgammaf, digammaf
      mu2i = thi
      if(model.eq.0) mu2i = mu2i*mu2i
      mu2i = max(0.d0,mu2i-df)
      z1 = df + mu2i
      z2 = z1 + mu2i
      fi = z1*z1/z2
      lgfi = lgammaf(fi/2.d0)
      RETURN
      END
      subroutine ncstats0(thi,lgfi0,df,model,lgfi,dgfi,fici)
      implicit logical (a-z)
      integer model
      real*8 thi,lgfi0,lgfi,dgfi,fici,df
      real*8 mu2i,z1,z2,dlci,fi,ci
      real*8 digammaf
      external digammaf
      mu2i = thi
      if(model.eq.0) mu2i = mu2i*mu2i
      mu2i = max(0.d0,mu2i-df)
      z1 = df + mu2i
      z2 = z1 + mu2i
      ci = z2/z1
      fi = z1/ci
      fici = fi*ci
      dlci = dlog(ci)
      dgfi = digammaf(0.5d0*fi)+dlci
      lgfi = lgfi0+0.5d0*(fi*dlci+fi-fi*dgfi)
      RETURN
      END
      real*8 function kldisnc1(lgfi,dgfi,fici,thj,lgfj,df,
     1                          model)
C    for smoothing noncentral Chi values
C    thi,thj  current estimates
C    df= 2 * number of coils
C    model = 0   smoothing of chi values
C    model = 1   smoothing of chi^2 values
      implicit logical (a-z)
      integer model
      real*8 lgfi,dgfi,fici,thj,df,lgfj
      real*8 mu2j,fj,cj,z1,z2
C  use Chi^2 instead of Chi for KL distance
      mu2j = thj
      if(model.eq.0) mu2j = mu2j*mu2j
      mu2j = max(0.d0,mu2j-df)
C  Approximation by Patnaik (1949)
      z1 = df + mu2j
      z2 = z1 + mu2j
      cj = z2/z1
      fj = z1/cj
      kldisnc1 = lgfj-lgfi+0.5d0*(fj*dlog(cj)+
     1                 fici/cj-fj*dgfi)
      RETURN
      END
      subroutine awsvchi2(y,th,ni,mask,n1,n2,n3,lambda,ncoils,
     1                    thn,th2,ni2,ldf,h,vext)
C   model=1 takes noncentral Chi-sq values in y0
C   model=0 takes noncentral Chi values in y0
C   perform adaptive smoothing on R^3
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C
      implicit logical (a-z)
      integer n1,n2,n3,ncoils
      logical mask(n1,n2,n3)
      real*8 y(n1,n2,n3),th(n1,n2,n3),ni(1),thn(1),
     1       th2(1),ni2(1),lambda,h,vext(2),
     2       ldf(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,cw1,cw2,cw3,i,n
      real*8 z,lgfi,dgfi,fici,df,kval,w,w0,h2,sw,sw2,swy,swy2,yj,
     1       z1,z2,z3
      real*8 kldisnc1
      external kldisnc1
      df=2.d0*ncoils
      h2=h*h
      cw1=h
      cw2=h/vext(1)
      cw3=h/vext(2)
      n = n1*n2*n3
C  precompute values of lgamma(corrected df/2) in each voxel
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,n3,ncoils,mask,y,th,ni,thn,th2,ni2,h,vext,ldf,n)
C$OMP& FIRSTPRIVATE(lambda,cw1,cw2,cw3,df,h2)
C$OMP& PRIVATE(i1,i2,i3,j1,j2,j3,i,z,lgfi,dgfi,fici,kval,w,w0,
C$OMP& sw,sw2,swy,swy2,yj,z1,z2,z3)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1         
         call lgstats(th(i1,i2,i3),df,1,ldf(i1,i2,i3))
      END DO
C$OMP END DO
C$OMP BARRIER
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1         
         if(.not.mask(i1,i2,i3)) CYCLE
         sw=0.d0
         swy=0.d0
         sw2=0.d0
         swy2=0.d0
         kval = lambda/ni(i)
         call ncstats0(th(i1,i2,i3),ldf(i1,i2,i3),df,
     1                       1,lgfi,dgfi,fici)
         DO j1=max(1,i1-cw1),min(i1+cw1,n1)
            z1=j1-i1
            z1=z1*z1
            DO j2=max(1,i2-cw2),min(i2+cw2,n2)
               z2=(j2-i2)*vext(1)
               z2=z2*z2
               DO j3=max(1,i3-cw3),min(i3+cw3,n3)
                  if(.not.mask(j1,j2,j3)) CYCLE
                  z3=(j3-i3)*vext(2)
                  w0=z1+z2+z3*z3
                  if(w0.ge.h2) CYCLE
                  w0=1.d0-w0/h2 
                  z=kldisnc1(lgfi,dgfi,fici,th(j1,j2,j3),
     1                       ldf(j1,j2,j3),df,1)
                  if(z.ge.kval) CYCLE
                  w=w0*min(1.d0,2.d0-2.d0*z/kval)
                  sw=sw+w
                  sw2=sw2+w*w
                  yj=y(j1,j2,j3)
                  swy=swy+w*yj
                  swy2=swy2+w*yj*yj
               END DO
            END DO
         END DO
         thn(i) = swy/sw
         th2(i) = swy2/sw
         ni(i) = sw
         ni2(i) = sw2
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni,ni2,th2)
      RETURN
      END
      subroutine paramw3(h,vext,ind,w,n)
C  compute a description of local weights 
C  h    - bandwidth
C  vext - vector (length 2) of relative voxel extensions
C  ind  - integer array dim (3,n) containing relative indices in xyz
C  w    - vector of corresponding weights
C  n    - number of positive weights (initial value 
C         (2 int(h)+1)*(2 int(h/vext(1))+1)*(2 int(h/vext(2))+1)
      integer n,ind(3,n)
      real*8 h,vext(2),w(n)
      integer i,i1,i2,i3,ih1,ih2,ih3
      real*8 hsq,z1,z2,z3
      hsq=h*h
      ih1 = h
      ih2 = h/vext(1)
      ih3 = h/vext(2)
      i=1
      DO i1=-ih1,ih1
         z1=i1*i1
         DO i2=-ih2,ih2
            z2=i2*vext(1)
            z2=z1+z2*z2
            IF(z2.ge.hsq) CYCLE
            DO i3=-ih3,ih3
               z3=i3*vext(2)
               z3=z2+z3*z3
               IF(z3.ge.hsq) CYCLE
               ind(1,i)=i1
               ind(2,i)=i2
               ind(3,i)=i3
               w(i)=1.d0-z3/hsq
               i=i+1
            END DO
         END DO
      END DO
      n=i-1
      RETURN
      END
      subroutine awsvchi(y,th,ni,mask,n1,n2,n3,ind,w,nw,lambda,
     1                    sigma,thn,sy)
C   Takes noncentral Chi values in y
C   perform adaptive smoothing on R^3
C   th containes previous estimates
C   ni containes previous sum of weights divided by variance of Chi(2,th/sigma)
C   mask - logical mask (use if mask==TRUE)
C   n1,n2,n3 - dimensions
C   ind  - integer array dim (3,n) containing relative indices in xyz
C   w    - vector of corresponding location weights
C   nw   - number of positive weights (initial value 
C   lambda   - kritical value for pairwise tests
C   sigma    - actual estimate of sigma
C   thn      - new estimate sum_j w_a(j) Y_j
C   th2      - sum_j w_a(j) Y_j^2
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C
      implicit logical (a-z)
      integer n1,n2,n3,nw,ind(3,nw)
      logical mask(n1,n2,n3)
      real*8 y(n1,n2,n3),th(n1,n2,n3),ni(1),thn(1),
     1       sy(1),lambda,w(nw),sigma
      integer i1,i2,i3,j1,j2,j3,i,j,n
      real*8 z,sw,sw2,swy,swy2,yj,thi,wj,kval,cw
      n = n1*n2*n3
C  precompute values of lgamma(corrected df/2) in each voxel
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,j,i1,i2,i3,j1,j2,j3,z,sw,swy,sw2,swy2,thi,kval,
C$OMP& wj,yj,cw)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1         
         if(.not.mask(i1,i2,i3)) CYCLE
         sw=0.d0
         swy=0.d0
         sw2=0.d0
         swy2=0.d0
         thi = th(i1,i2,i3)
C   thats the estimated standard deviation of y(i1,i2,i3)
         kval = lambda/ni(i)*sigma*sigma
         Do j=1,nw
            j1=i1+ind(1,j)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,j)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,j)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            wj=w(j)
            z=thi-th(j1,j2,j3)
            z=z*z
            if(z.ge.kval) CYCLE
            wj=wj*min(1.d0,2.d0-2.d0*z/kval)
            sw=sw+wj
            sw2=sw2+wj*wj
            yj=y(j1,j2,j3)
            swy=swy+wj*yj
            swy2=swy2+wj*yj*yj
         END DO
         thi = swy/sw
         z = swy2/sw
C  z2-thi^2  is an estimate of the variance of y(i) 
         cw = 1.d0-sw2/sw/sw
         IF(cw.gt.0.d0) THEN
            sy(i) = sqrt((z-thi*thi)/cw)
C  sy(i)  is an estimate of sigma corrected for 
C       simultaneously estimating the mean and for non-central chi-bias 
         ELSE
            sy(i) = 0.d0
C  case ni(i) = 1
         END IF
         thn(i) = thi
         ni(i) = sw
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni,sy)
      RETURN
      END
      subroutine awsadchi(y,th,ni,mask,n1,n2,n3,ind,w,nw,lambda,
     1                    sigma,wad,nthreds,thn,sy)
C   Takes noncentral Chi values in y
C   perform adaptive smoothing on R^3
C   th containes previous estimates
C   ni containes previous sum of weights
C   mask - logical mask (use if mask==TRUE)
C   n1,n2,n3 - dimensions
C   ind  - integer array dim (3,n) containing relative indices in xyz
C   w    - vector of corresponding location weights
C   nw   - number of positive weights (initial value 
C   lambda   - kritical value for pairwise tests
C   sigma    - actual estimate of sigma
C   thn      - new estimate sum_j w_a(j) Y_j
C   th2      - sum_j w_a(j) Y_j^2
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C
      implicit logical (a-z)
      integer n1,n2,n3,nw,ind(3,nw),nthreds
      logical mask(n1,n2,n3)
      real*8 y(n1,n2,n3),th(n1,n2,n3),ni(1),thn(1),
     1       sy(1),lambda,w(nw),sigma,wad(nw,nthreds)
      integer i1,i2,i3,j1,j2,j3,i,j,n,thrednr
      real*8 z,sw,sw2,swy,swy2,yj,thi,wj,kval,cw
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      n = n1*n2*n3
      thrednr = 1
C  precompute values of lgamma(corrected df/2) in each voxel
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,j,i1,i2,i3,j1,j2,j3,z,sw,swy,sw2,swy2,thi,kval,
C$OMP& wj,yj,cw,thrednr)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1         
         if(.not.mask(i1,i2,i3)) CYCLE
!$         thrednr = omp_get_thread_num()+1
         sw=0.d0
         swy=0.d0
         sw2=0.d0
         swy2=0.d0
         thi = th(i1,i2,i3)
C   thats the estimated standard deviation of y(i1,i2,i3)
         kval = lambda/ni(i)*sigma*sigma
         DO j=1,nw
            wad(j,thrednr)=0.d0
            j1=i1+ind(1,j)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,j)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,j)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            wj=w(j)
            z=thi-th(j1,j2,j3)
            z=z*z
            if(z.ge.kval) CYCLE
            wj=wj*min(1.d0,2.d0-2.d0*z/kval)
            wad(j,thrednr)=wj
            sw=sw+wj
            sw2=sw2+wj*wj
            yj=y(j1,j2,j3)
            swy=swy+wj*yj
         END DO
         thi = swy/sw
         DO j=1,nw
            wj=wad(j,thrednr)
            if(wj.le.1d-8) CYCLE
            j1=i1+ind(1,j)
            j2=i2+ind(2,j)
            j3=i3+ind(3,j)
C no need to test for grid coordinates since wj>0
            swy2=swy2+wj*abs(thi-y(j1,j2,j3))
         END DO
         z = swy2/sw/.8d0
C  z  is an estimate of the standard deviation of y(i) 
         cw = 1.d0-sw2/sw/sw
         IF(cw.gt.0.d0) THEN
            sy(i) = z/sqrt(cw)
C  sy(i)  is an estimate of sigma corrected for 
C       simultaneously estimating the mean and for non-central chi-bias 
         ELSE
            sy(i) = 0.d0
C  case ni(i) = 1
         END IF
         thn(i) = thi
         ni(i) = sw
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni,sy)
      RETURN
      END
      subroutine lkfuls0(h,vext,ind,wght,n)
      implicit logical (a-z)
      integer n,ind(3,n)
      real*8 h,vext(2),wght(n)
      integer ih1,ih2,ih3,i,j1,j2,j3
      real*8 h2,x1,x2,x3,z,z1,vd2,vd3
      vd2 = vext(1)
      vd3 = vext(2)
      ih1 = max(1.d0,5.0d0*h)
      ih2 = max(1.d0,5.0d0*h/vd2)
      ih3 = max(1.d0,5.0d0*h/vd3)
      h2 = h*h
      i = 1
      z = 0.d0
C  just to prevent compiler warnings
      DO j1 = 0,ih1
         x1 = j1
         DO j2 = -ih2,ih2
            x2 = vd2*j2
            DO j3 = -ih3,ih3
               x3 = vd3*j3
               z1 = z+x1*x1+x2*x2+x3*x3
               if(z1.gt.h2) CYCLE
               if(i.gt.n) THEN
                  call intpr("Exceeded max i",14,i,1)
                  n = i-1
                  return
               END IF
               wght(i)= (1.d0-z1/h2)
               ind(1,i) = j1
               ind(2,i) = j2
               ind(3,i) = j3
               i = i+1
            END DO
         END DO
      END DO
      n = i-1
      RETURN
      END
      subroutine ipolsp(theta,th0,ni,ni0,n,ng,gind,gw,nbv,nbvp1,
     1                    msth,msni)
C   interpolate values of theta on spheres where it was not observed
      implicit logical (a-z)
      integer n,ng,nbv,nbvp1,gind(3,nbv,ng)
      real*8 theta(n,ng),th0(n),ni(n,ng),ni0(n),gw(3,nbv,ng),
     1       msth(nbvp1,n,ng),msni(nbvp1,n,ng)
      integer i,j,k,i1,i2,i3,ip1
      real*8 w1,w2,w3
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,j,k,w1,w2,w3,i1,i2,i3,ip1)
C$OMP DO SCHEDULE(GUIDED)
      DO j=1,ng
         DO k=1,n
            msth(1,k,j) = th0(k)
            msni(1,k,j) = ni0(k)
            DO i=1,nbv
               ip1 = i+1
               i1 = gind(1,i,j)
               IF(i1.eq.j) THEN
C  same shell just copy
                  msth(ip1,k,j) = theta(k,j)
                  msni(ip1,k,j) = ni(k,j)
               ELSE
C  different shell need to interpolate
                  i2 = gind(2,i,j)
                  i3 = gind(3,i,j)
                  w1 = gw(1,i,j)
                  w2 = gw(2,i,j)
                  w3 = gw(3,i,j)
                  msth(ip1,k,j) = theta(k,i1)*w1+theta(k,i2)*w2+
     1                            theta(k,i3)*w3
                  msni(ip1,k,j) = 1.d0/
     1                           (w1/ni(k,i1)+w2/ni(k,i2)+w3/ni(k,i3))
               END IF
            END DO
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(msth,msni)
      RETURN
      END
