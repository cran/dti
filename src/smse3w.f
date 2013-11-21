      subroutine adsmse3w(y,y0,th,ni,th0,ni0,fsi2,fsi02,mask,ns,n1,
     1                n2,n3,ngrad,lambda,ws0,ncores,ind,w,n,ind0,w0,
     2                n0,thn,nin,th0n,ni0n,sw,swy,thi,nii,fsi2i,
     3                ix,iy,iz,aw,aw0,nw,nw0)
C   
C  Multi-shell version (differs in dimension of th 
C  KL-distance based on all spheres and Gauss-approximation only
C  see ~polzehl/latex/1211_simplemetric/Rcode/Figuregaussapprox.r 
C  for approximative KL-distance between non-central chi distributions
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
      integer ns,n1,n2,n3,ngrad,n,n0,ind(5,n),ind0(3,n0),ncores,
     1        ix,iy,iz,nw,nw0
      logical mask(*)
      real*8 y(*),y0(*),th(ns,*),ni(ns,*),th0(ns,*),ni0(ns,*),
     1     fsi2(ns,*),fsi02(ns,*),thn(*),th0n(*),nin(*),ni0n(*),
     2     aw(nw),aw0(nw0)
C  * refers to n1*n2*n3*ngrad for y,th,ni,thn,fsi2,nin and to
C              n1*n2*n3 for y0,th0,ni0,th0n,ni0n,fsi02,mask
      real*8 w(n),w0(n0),lambda,thi(*),ws0,fsi2i(*),sw(*),swy(*),nii(*)
C  * refers to ns*ncores in thi, fsi2i, nii and to
C              ngrad*cores in sw and swy
      integer iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,k,jind,iind4,
     1        jind4,n123,n12,sthrednr,gthrednr,i4gthnr,m1,m01,ixyz
      real*8 sz,z,sw0,swy0
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      thrednr = 1
      n12 = n1*n2
      n123 = n12*n3
      ixyz = ix+(iy-1)*n1+(iz-1)*n12
      DO i=1,nw
         aw(i)=0.d0
      END DO
      DO i=1,nw0
         aw0(i)=0.d0
      END DO
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ns,n1,n2,n3,ngrad,n,n0,ind,ind0,ncoils,ncores,y,y0,
C$OMP&       th,ni,th0,ni0,w,w0,thn,th0n,nin,ni0n,thi,sw,swy,nii,
C$OMP&       lambda,mask,ws0,fsi2,fsi02,fsi2i,ix,iy,iz,nw,nw0,aw,aw0)
C$OMP& FIRSTPRIVATE(n123,n12,ixyz)
C$OMP& PRIVATE(iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,k,sz,z,
C$OMP&   sw0,swy0,jind,iind4,jind4,sthrednr,gthrednr,i4gthnr,m1,m01)
C$OMP DO SCHEDULE(GUIDED)
C  First si - images
      DO iind=1,n1*n2*n3
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         if(ixyz.eq.iind) THEN
            m1=0
            m01=0
         END IF
         sthrednr = (thrednr-1)*ns
         gthrednr = (thrednr-1)*ngrad
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1         
         if(.not.mask(iind)) CYCLE
         DO i4=1,ngrad
            sw(i4+gthrednr)=0.d0
            swy(i4+gthrednr)=0.d0
         END DO
         i4=0
         i4gthnr=gthrednr
         DO i=1,n
            if(ixyz.eq.iind) m1=m1+1
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               i4gthnr=i4+gthrednr
               iind4 = iind+(i4-1)*n123
               DO k=1,ns
                  fsi2i(k+sthrednr)=fsi2(k,iind4)
                  thi(k+sthrednr) = th(k,iind4)
                  nii(k+sthrednr) = ni(k,iind4)/lambda
               END DO
               nii(1+sthrednr)=ws0*nii(1+sthrednr)
            END IF
            j1=i1+ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            jind=j1+(j2-1)*n1+(j3-1)*n12
            if(.not.mask(jind)) CYCLE          
            j4=ind(5,i)
            jind4=jind+(j4-1)*n123
C adaptation 
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k+sthrednr)-th(k,jind4))
                  sz=sz+nii(k+sthrednr)*z*z/
     1                        (fsi2(k,jind4)+fsi2i(k+sthrednr))
               END DO
C  do not adapt on the sphere !!! 
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w(i)
            if(sz.gt.0.5d0) z=z*(2.d0-2.d0*sz)
            if(ixyz.eq.iind) aw(m1)=z
            sw(i4gthnr)=sw(i4gthnr)+z
            swy(i4gthnr)=swy(i4gthnr)+z*y(jind4)
         END DO
C  now opposite directions
         DO i=1,n
            if(ind(1,i).eq.0) CYCLE
            if(ixyz.eq.iind) m1=m1+1
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               i4gthnr=i4+gthrednr
               iind4 = iind+(i4-1)*n123
               DO k=1,ns
                  fsi2i(k+sthrednr)=fsi2(k,iind4)
                  thi(k+sthrednr) = th(k,iind4)
                  nii(k+sthrednr) = ni(k,iind4)/lambda
               END DO
               nii(1+sthrednr)=ws0*nii(1+sthrednr)
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
            jind=j1+(j2-1)*n1+(j3-1)*n12
            if(.not.mask(jind)) CYCLE          
            j4=ind(5,i)
            jind4=jind+(j4-1)*n123
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k+sthrednr)-th(k,jind4))
                  sz=sz+nii(k+sthrednr)*z*z/
     1                        (fsi2(k,jind4)+fsi2i(k+sthrednr))
               END DO
C  do not adapt on the sphere !!! 
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w(i)
            if(sz.gt.0.5d0) z=z*(2.d0-2.d0*sz)
            if(ixyz.eq.iind) aw(m1)=z
            sw(i4gthnr)=sw(i4gthnr)+z
            swy(i4gthnr)=swy(i4gthnr)+z*y(jind4)
         END DO
         DO i4=1,ngrad
            iind4 = iind+(i4-1)*n123
            thn(iind4) = swy(i4+gthrednr)/sw(i4+gthrednr)
            nin(iind4) = sw(i4+gthrednr)
         END DO
C    now the s0 image in iind
         sw0=0.d0
         swy0=0.d0
         DO k=1,ns
            thi(k+sthrednr) = th0(k,iind)
            nii(k+sthrednr) = ni0(k,iind)/lambda
            fsi2i(k+sthrednr)=fsi02(k,iind)
         END DO
         DO i=1,n0
            if(ixyz.eq.iind) m01=m01+1
            j1=i1+ind0(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind0(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind0(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            jind=j1+(j2-1)*n1+(j3-1)*n12
            if(.not.mask(jind)) CYCLE          
C adaptation 
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k+sthrednr)-th0(k,jind))
                  sz=sz+nii(k+sthrednr)*z*z/
     1                        (fsi02(k,jind)+fsi2i(k+sthrednr))
               END DO
C  do not adapt on the sphere !!! 
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w0(i)
            if(sz.gt.0.5d0) z=z*(2.d0-2.d0*sz)
            if(ixyz.eq.iind) aw0(m01)=z
            sw0=sw0+z
            swy0=swy0+z*y0(jind)
         END DO
C  now opposite directions
         DO i=1,n0
            if(ind0(1,i).eq.0) CYCLE
C
C   handle case j1-i1 < 0 which is not contained in ind 
C   using axial symmetry
C
            if(ixyz.eq.iind) m01=m01+1
            j1=i1-ind0(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2-ind0(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3-ind0(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            jind=j1+(j2-1)*n1+(j3-1)*n12
            if(.not.mask(jind)) CYCLE          
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k+sthrednr)-th0(k,jind))
                  sz=sz+nii(k+sthrednr)*z*z/
     1                        (fsi02(k,jind)+fsi2i(k+sthrednr))
               END DO
C  do not adapt on the sphere !!! 
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w0(i)
            if(sz.gt.0.5d0) z=z*(2.d0-2.d0*sz)
            if(ixyz.eq.iind) aw0(m01)=z
            sw0=sw0+z
            swy0=swy0+z*y0(jind)
         END DO
         th0n(iind) = swy0/sw0
         ni0n(iind) = sw0
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,nin,th0n,ni0n)
      RETURN
      END
