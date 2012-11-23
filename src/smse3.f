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
      vd2 = vext(1)
      vd3 = vext(2)
      i = 1
      z = 0.d0
C  just to prevent compiler warnings
      DO j4 = 1,ng
         k4 = k456(1,i4,j4)
         k5 = k456(2,i4,j4)
         k6 = k456(3,i4,j4)
         if(dist.eq.1) z = (k4*k4+k5*k5+abs(k6))/kap2
         if(dist.eq.2) z = (k4*k4+k5*k5+k6*k6)/kap2
         if(dist.eq.3) z = k4*k4/kap2
         if(z.gt.h2) CYCLE
C   last three komponents already to large
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
                     call intpr("for i4",6,i4,1)
                     n = i-1
                     return
                  END IF
                  wght(i)= (1.d0-z1/h2)
C   *cb
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
      subroutine lkfse3i0(h,kappa,i4,k456,ng,
     1                   vext,vred,n,dist)
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
      vd2 = vext(1)
      vd3 = vext(2)
      n = 0
      z = 0.d0
C  just to prevent compiler warnings
      DO j4 = 1,ng
         k4 = k456(1,i4,j4)
         k5 = k456(2,i4,j4)
         k6 = k456(3,i4,j4)
         if(dist.eq.1) z = (k4*k4+k5*k5+abs(k6))/kap2
         if(dist.eq.2) z = (k4*k4+k5*k5+k6*k6)/kap2
         if(dist.eq.3) z = k4*k4/kap2
         if(z.gt.h2) CYCLE
C   last three komponents already to large
         DO j1 = 0,ih1
            x1 = j1
            DO j2 = -ih2,ih2
               x2 = vd2*j2
               DO j3 = -ih3,ih3
                  x3 = vd3*j3
                  z1 = z+x1*x1+x2*x2+x3*x3
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
      subroutine adrsmse3(y,th,ni,mask,n1,n2,n3,ngrad,lambda,ncoils,
     1                    ind,w,n,thn,ldf,sigma,sw,swy,model)
C   model=1 takes noncentral Chi-sq values in y
C   model=0 takes noncentral Chi values in y
C
C   perform adaptive smoothing on SE(3) 
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C
      implicit logical (a-z)
      integer n1,n2,n3,ngrad,n,ind(5,n),ncoils,model
      logical mask(n1,n2,n3)
      real y(n1,n2,n3,ngrad),thn(n1,n2,n3,ngrad),ni(n1,n2,n3,ngrad),
     1     ldf(n1,n2,n3,ngrad),th(n1,n2,n3,ngrad)
      real*8 lambda,w(n),sw(ngrad),swy(ngrad),
     2       sigma,lgfi,dgfi,fici,df
      integer i,i1,i2,i3,i4,j1,j2,j3,j4
      real*8 z,thi,nii,thj,ldfi,ldfj
      real*8 kldisnc1
      external kldisnc1
      nii=1.d0
      df=2.d0*ncoils
C just to prevent a compiler warning
C  precompute values of lgamma(corrected df/2) in each voxel
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               DO i4=1,ngrad
                  thi=th(i1,i2,i3,i4)
                  call lgstats(thi,sigma,df,model,ldfi)
                  ldf(i1,i2,i3,i4)=ldfi
               END DO
            END DO
         END DO
      END DO
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(.not.mask(i1,i2,i3)) CYCLE
               DO i4=1,ngrad
                  sw(i4)=0.d0
                  swy(i4)=0.d0
               END DO
               i4=0
               DO i=1,n
                  if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
                     i4 = ind(4,i)
                     thi = th(i1,i2,i3,i4)
                     nii = ni(i1,i2,i3,i4)/lambda
                     ldfi=ldf(i1,i2,i3,i4)
                     call ncstats0(thi,ldfi,sigma,df,
     1                             model,lgfi,dgfi,fici)
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
                     ldfj=ldf(j1,j2,j3,j4)
                     z=nii*kldisnc1(lgfi,dgfi,fici,thj,
     1                              ldfj,sigma,df,model)
C  do not adapt on the sphere !!! 
                  ELSE
                     z=0.d0
                  END IF
                  if(z.ge.1.d0) CYCLE
                  z=w(i)*min(1.d0,2.d0-2.d0*z)
                  sw(i4)=sw(i4)+z
                  swy(i4)=swy(i4)+z*y(j1,j2,j3,j4)
               END DO
               DO i=1,n
                  if(ind(1,i).eq.0) CYCLE
                  if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
                     i4 = ind(4,i)
                     thi = th(i1,i2,i3,i4)
                     nii = ni(i1,i2,i3,i4)/lambda
                     ldfi=ldf(i1,i2,i3,i4)
                     call ncstats0(thi,ldfi,sigma,df,
     1                             model,lgfi,dgfi,fici)
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
                     ldfj=ldf(j1,j2,j3,j4)
                     z=nii*kldisnc1(lgfi,dgfi,fici,thj,
     1                              ldfj,sigma,df,model)
C  do not adapt on the sphere !!! 
                  ELSE
                     z=0.d0
                  END IF
                  if(z.ge.1.d0) CYCLE
                  z=w(i)*min(1.d0,2.d0-2.d0*z)
                  sw(i4)=sw(i4)+z
                  swy(i4)=swy(i4)+z*y(j1,j2,j3,j4)
               END DO
               DO i4=1,ngrad
                  thn(i1,i2,i3,i4) = swy(i4)/sw(i4)
                  ni(i1,i2,i3,i4) = sw(i4)
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine adsmse3p(y,th,ni,mask,n1,n2,n3,ngrad,lambda,ncoils,
     1                    ncores,ind,w,n,thn,ldf,sigma,sw,swy,model)
C   model=1 takes noncentral Chi-sq values in y
C   model=0 takes noncentral Chi values in y
C
C   perform adaptive smoothing on SE(3) 
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C
      implicit logical (a-z)
      integer n1,n2,n3,ngrad,n,ind(5,n),ncoils,model,ncores
      logical mask(n1,n2,n3)
      real y(n1,n2,n3,ngrad),thn(n1,n2,n3,ngrad),ni(n1,n2,n3,ngrad),
     1     ldf(n1,n2,n3,ngrad),th(n1,n2,n3,ngrad)
      real*8 lambda,w(n),sw(ngrad,ncores),swy(ngrad,ncores),
     2       sigma,lgfi,dgfi,fici,df
      integer iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr
      real*8 z,thi,nii,thj,ldfi,ldfj
      integer omp_get_thread_num
      real*8 kldisnc1
      external kldisnc1,omp_get_thread_num
      nii=1.d0
      df=2.d0*ncoils
C just to prevent a compiler warning
C  precompute values of lgamma(corrected df/2) in each voxel
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,n3,ngrad,ncores,mask,y,thn,ni,ldf,th,w,sw,swy,
C$OMP& ind,ncoils,model,df,n,lambda,sigma)
C$OMP& PRIVATE(iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,z,thi,nii,thj,
C$OMP& ldfi,ldfj,dgfi,fici,lgfi)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1
         if(.not.mask(i1,i2,i3)) CYCLE
C      DO i1=1,n1
C         DO i2=1,n2
C            DO i3=1,n3
               DO i4=1,ngrad
                  thi=th(i1,i2,i3,i4)
                  call lgstats(thi,sigma,df,model,ldfi)
                  ldf(i1,i2,i3,i4)=ldfi
               END DO
C            END DO
C         END DO
      END DO
C$OMP END DO 
C$OMP BARRIER
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1         
C      DO i1=1,n1
C         DO i2=1,n2
C            DO i3=1,n3
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
                     nii = ni(i1,i2,i3,i4)/lambda
                     ldfi=ldf(i1,i2,i3,i4)
                     call ncstats0(thi,ldfi,sigma,df,
     1                             model,lgfi,dgfi,fici)
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
                     ldfj=ldf(j1,j2,j3,j4)
                     z=nii*kldisnc1(lgfi,dgfi,fici,thj,
     1                              ldfj,sigma,df,model)
C  do not adapt on the sphere !!! 
                  ELSE
                     z=0.d0
                  END IF
                  if(z.ge.1.d0) CYCLE
                  z=w(i)*min(1.d0,2.d0-2.d0*z)
                  sw(i4,thrednr)=sw(i4,thrednr)+z
                  swy(i4,thrednr)=swy(i4,thrednr)+z*y(j1,j2,j3,j4)
               END DO
               DO i=1,n
                  if(ind(1,i).eq.0) CYCLE
                  if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
                     i4 = ind(4,i)
                     thi = th(i1,i2,i3,i4)
                     nii = ni(i1,i2,i3,i4)/lambda
                     ldfi=ldf(i1,i2,i3,i4)
                     call ncstats0(thi,ldfi,sigma,df,
     1                             model,lgfi,dgfi,fici)
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
                     ldfj=ldf(j1,j2,j3,j4)
                     z=nii*kldisnc1(lgfi,dgfi,fici,thj,
     1                              ldfj,sigma,df,model)
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
C               call rchkusr()
C            END DO
C         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni)
      RETURN
      END
      subroutine asmse3s0(y,th,ni,mask,n1,n2,n3,ns0,
     1                    lambda,ncoils,ind,w,n,starts,nstarts,
     2                    thn,ldf,sigma,swi,model)
C   model=1 takes noncentral Chi-sq values in y0
C   model=0 takes noncentral Chi values in y0
C   perform adaptive smoothing on R^3
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C
      implicit logical (a-z)
      integer n1,n2,n3,n,ind(5,n),ns0,starts(1),nstarts,ncoils,model
      logical mask(n1,n2,n3)
      real*8 y(n1,n2,n3),th(n1,n2,n3),ni(n1,n2,n3),thn(n1,n2,n3),
     1       lambda,w(n),sw0,swy0,sigma,swi(nstarts),ldf(n1,n2,n3)
      integer i,i0,i1,i2,i3,j1,j2,j3,l1,l2,l3
      real*8 z,lgfi,dgfi,fici,ns0sc,df
      real*8 kldisnc1
      external kldisnc1
      df=2.d0*ncoils
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               call lgstats(th(i1,i2,i3),sigma,df,model,ldf(i1,i2,i3))
            END DO
         END DO
      END DO
      DO i0=1,nstarts
         swi(i0)=0.d0
         DO i=starts(i0)+1,starts(i0+1)
            swi(i0)=swi(i0)+w(i)
         END DO
         swi(i0)=swi(i0)/(starts(i0+1)-starts(i0))
      END DO
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(.not.mask(i1,i2,i3)) CYCLE
               sw0=0.d0
               swy0=0.d0
               ns0sc = ns0*lambda/ni(i1,i2,i3)
               call ncstats0(th(i1,i2,i3),ldf(i1,i2,i3),sigma,df,
     1                       model,lgfi,dgfi,fici)
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
                  z=ns0*kldisnc1(lgfi,dgfi,fici,th(j1,j2,j3),
     1                            ldf(j1,j2,j3),sigma,df,model)
                  if(z.ge.ns0sc) CYCLE
                  if(z.lt.ns0sc) THEN
                     z=swi(i0)*min(1.d0,2.d0-2.d0*z/ns0sc)
                     sw0=sw0+z
                     swy0=swy0+z*y(j1,j2,j3)
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
                  z=ns0*kldisnc1(lgfi,dgfi,fici,th(j1,j2,j3),
     1                           ldf(j1,j2,j3),sigma,df,model)
                  if(z.ge.ns0sc) CYCLE
                  if(z.lt.ns0sc) THEN
                     z=swi(i0)*min(1.d0,2.d0-2.d0*z/ns0sc)
                     sw0=sw0+z
                     swy0=swy0+z*y(j1,j2,j3)
                  END IF
               END DO
               thn(i1,i2,i3) = swy0/sw0
               ni(i1,i2,i3) = sw0
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine asmse30p(y,th,ni,mask,n1,n2,n3,ns0,lambda,ncoils,
     1                ind,w,n,starts,nstarts,thn,ldf,sigma,swi,model)
C   model=1 takes noncentral Chi-sq values in y0
C   model=0 takes noncentral Chi values in y0
C   perform adaptive smoothing on R^3
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively 
C
      implicit logical (a-z)
      integer n1,n2,n3,n,ind(5,n),ns0,starts(1),nstarts,ncoils,model
      logical mask(n1,n2,n3)
      real*8 y(n1,n2,n3),th(n1,n2,n3),ni(1),thn(1),
     1       lambda,w(n),sw0,swy0,sigma,swi(nstarts),ldf(n1,n2,n3)
      integer i,i0,i1,i2,i3,j1,j2,j3,l1,l2,l3,nn
      real*8 z,ldfi,lgfi,dgfi,fici,ns0sc,df
      real*8 kldisnc1
      external kldisnc1
      df=2.d0*ncoils
      nn=n1*n2*n3
      DO i0=1,nstarts
         swi(i0)=0.d0
         DO i=starts(i0)+1,starts(i0+1)
            swi(i0)=swi(i0)+w(i)
         END DO
         swi(i0)=swi(i0)/(starts(i0+1)-starts(i0))
      END DO
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,n3,nn,ncoils,mask,y,th,ni,thn,ind,ldf,n,nstarts,
C$OMP& starts,w,swi,df)
C$OMP& FIRSTPRIVATE(ns0,model,lambda,sigma)
C$OMP& PRIVATE(sw0,swy0,i,i0,i1,i2,i3,j1,j2,j3,l1,l2,l3,z,lgfi,dgfi,
C$OMP& ldfi,fici,ns0sc)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,nn
C      DO i1=1,n1
C         DO i2=1,n2
C            DO i3=1,n3
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1
         if(.not.mask(i1,i2,i3)) CYCLE
         call lgstats(th(i1,i2,i3),sigma,df,model,ldfi)
         ldf(i1,i2,i3) = ldfi
C            END DO
C         END DO
      END DO
C$OMP END DO
C$OMP BARRIER
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,nn
C      DO i1=1,n1
C         DO i2=1,n2
C            DO i3=1,n3
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1         
               if(.not.mask(i1,i2,i3)) CYCLE
               sw0=0.d0
               swy0=0.d0
               ns0sc = ns0*lambda/ni(i)
               call ncstats0(th(i1,i2,i3),ldf(i1,i2,i3),sigma,df,
     1                       model,lgfi,dgfi,fici)
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
                  z=ns0*kldisnc1(lgfi,dgfi,fici,th(j1,j2,j3),
     1                            ldf(j1,j2,j3),sigma,df,model)
                  if(z.ge.ns0sc) CYCLE
                  if(z.lt.ns0sc) THEN
                     z=swi(i0)*min(1.d0,2.d0-2.d0*z/ns0sc)
                     sw0=sw0+z
                     swy0=swy0+z*y(j1,j2,j3)
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
                  z=ns0*kldisnc1(lgfi,dgfi,fici,th(j1,j2,j3),
     1                           ldf(j1,j2,j3),sigma,df,model)
                  if(z.ge.ns0sc) CYCLE
                  if(z.lt.ns0sc) THEN
                     z=swi(i0)*min(1.d0,2.d0-2.d0*z/ns0sc)
                     sw0=sw0+z
                     swy0=swy0+z*y(j1,j2,j3)
                  END IF
               END DO
               thn(i) = swy0/sw0
               ni(i) = sw0
C               call rchkusr()
C            END DO
C         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni)
      RETURN
      END
C
C   Kullback-leibler distance for noncentral-Chi-distributions
C   with parameters thi/sigma and thj/sigma and 2*nc degrees of freedom
C
C
C   variant with precomputed quantities
C
      subroutine lgstats(thi,sigma,df,model,lgfi)
      implicit logical (a-z)
      integer model
      real*8 thi,sigma,lgfi,df
      real*8 mu2i,z1,z2,fi
      real*8 lgammaf,digammaf
      external lgammaf, digammaf
      mu2i = thi/sigma
      if(model.eq.0) mu2i = mu2i*mu2i
      mu2i = max(0.d0,mu2i-df)
      z1 = df + mu2i
      z2 = z1 + mu2i
      fi = z1*z1/z2
      lgfi = lgammaf(fi/2.d0)
      RETURN
      END
      subroutine ncstats0(thi,lgfi0,sigma,df,model,lgfi,dgfi,fici)
      implicit logical (a-z)
      integer model
      real*8 thi,sigma,lgfi0,lgfi,dgfi,fici,df
      real*8 mu2i,z1,z2,dlci,fi,ci
      real*8 digammaf
      external digammaf
      mu2i = thi/sigma
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
      real*8 function kldisnc1(lgfi,dgfi,fici,thj,lgfj,sigma,df,
     1                          model)
C    for smoothing noncentral Chi values
C    thi,thj  current estimates
C    sigma    estimated scale parameter 
C    df= 2 * number of coils
C    model = 0   smoothing of chi values
C    model = 1   smoothing of chi^2 values
      implicit logical (a-z)
      integer model
      real*8 lgfi,dgfi,fici,thj,sigma,df,lgfj
      real*8 mu2j,fj,cj,z1,z2
C  use Chi^2 instead of Chi for KL distance
      mu2j = thj/sigma
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
     1                    thn,th2,ni2,ldf,sigma,h,vext)
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
     1       th2(1),ni2(1),lambda,sigma,h,vext(2),
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
C$OMP& FIRSTPRIVATE(lambda,sigma,cw1,cw2,cw3,df,h2)
C$OMP& PRIVATE(i1,i2,i3,j1,j2,j3,i,z,lgfi,dgfi,fici,kval,w,w0,
C$OMP& sw,sw2,swy,swy2,yj,z1,z2,z3)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
C      DO i1=1,n1
C         DO i2=1,n2
C            DO i3=1,n3
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1         
               call lgstats(th(i1,i2,i3),sigma,df,1,ldf(i1,i2,i3))
C            END DO
C         END DO
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
C      DO i1=1,n1
C         DO i2=1,n2
C           DO i3=1,n3
               if(.not.mask(i1,i2,i3)) CYCLE
               sw=0.d0
               swy=0.d0
               sw2=0.d0
               swy2=0.d0
               kval = lambda/ni(i)
               call ncstats0(th(i1,i2,i3),ldf(i1,i2,i3),sigma,df,
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
     1                             ldf(j1,j2,j3),sigma,df,1)
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
C               call rchkusr()
C            END DO
C         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni,ni2,th2)
      RETURN
      END
      
