      subroutine mediansm(y,n1,n2,n3,ind,nind,work,ncores,yout)
C
C
C   3D median smoother of y with neighborhood defined by ind
C   results in yout
C   size of work needs to be 2*nind
C
      implicit logical (a-z)
      integer n1,n2,n3,nind,ind(3,nind),ncores
      real*8 y(n1,n2,n3),yout(n1,n2,n3),work(nind,ncores)
      integer i1,i2,i3,j1,j2,j3,j,k,thrednr
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i1,i2,i3,j1,j2,j3,k,thrednr)
C$OMP DO SCHEDULE(GUIDED)
      DO i1=1,n1
!$         thrednr = omp_get_thread_num()+1
         DO i2=1,n2
            DO i3=1,n3
               k=0
               DO j=1,nind
                  j1=i1+ind(1,j)
                  if(j1.le.0.or.j1.gt.n1) CYCLE
                  j2=i2+ind(2,j)
                  if(j2.le.0.or.j2.gt.n2) CYCLE
                  j3=i3+ind(3,j)
                  if(j3.le.0.or.j3.gt.n3) CYCLE
                  k=k+1
                  work(k,thrednr)=y(j1,j2,j3)
               END DO
               call qsort3(work(1,thrednr),1,k)
               IF (mod(k,2) == 0) THEN    
                  yout(i1,i2,i3) = 
     1               (work(k/2,thrednr)+work(k/2+1,thrednr))/2.d0
               ELSE
                  yout(i1,i2,i3) = work(k/2+1,thrednr)
               END IF
            END DO
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(yout)
      return
      end
      subroutine awslchi(s,th,ni,sigma,fns,L,mask,n1,n2,n3,ind,w,nw,
     1                minni,wad,sad,lambda,nthreds,iL,work,thn,sigman)
C
C  local variance estimation using (adaptive) weighted likelihood
C
C   Takes observed intensities in s and 
C     initial estimates of \sigma in sigma
C   perform adaptive smoothing on R^3
C   th containes previous estimates of E S
C   ni containes previous sum of weights
C   mask - logical mask (use if mask==TRUE)
C   n1,n2,n3 - dimensions
C   ind  - integer array dim (3,n) containing relative indices in xyz
C   w    - vector of corresponding location weights
C   nw   - number of positive weights (initial value 
C   
C   lambda   - kritical value for pairwise tests
C   thn      - new estimate sum_j w_a(j) S_j
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:3] are j1-i1,j2-i2 and j3-i3 respectively 
C   wad, sad - array for weights>0 and corresponding observed s
C
      implicit logical (a-z)
      integer n1,n2,n3,nw,ind(3,nw),nthreds,iL
      logical mask(n1,n2,n3)
      real*8 s(n1,n2,n3),th(n1,n2,n3),ni(n1*n2*n3),thn(n1*n2*n3),
     1  fns(n1,n2,n3),sigman(n1*n2*n3),lambda,w(nw),sigma(n1,n2,n3),
     2  wad(nw,nthreds),sad(nw,nthreds),L,minni,work(iL,nthreds)
      integer i1,i2,i3,j1,j2,j3,i,j,jj,n,maxit,thrednr
      real*8 z,sw,sws,sws2,sj,thi,wj,kval,fnsi,sgi,tol,low,up,ksi,
     1       xmin,fmin
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      n = n1*n2*n3
      thrednr = 1
      tol=1d-6
      maxit=100
C  precompute values of lgamma(corrected df/2) in each voxel
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& FIRSTPRIVATE(iL,L,minni,n1,n2,n3,maxit)
C$OMP& PRIVATE(i,j,i1,i2,i3,j1,j2,j3,z,sw,sws,sws2,thi,kval,
C$OMP& wj,sj,thrednr,fnsi,low,up,tol,sgi,ksi,jj,xmin,fmin)
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
         sws=0.d0
         sws2=0.d0
         sgi=sigma(i1,i2,i3)
         thi = th(i1,i2,i3)
         thn(i) = thi
         fnsi = fns(i1,i2,i3)
C   thats the estimated standard deviation of s(i1,i2,i3)
         kval = lambda/ni(i)*sgi*sgi
         jj = 0
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
            z=z*z/(fnsi+fns(j1,j2,j3))
            if(z.ge.kval) CYCLE
            wj=wj*min(1.d0,2.d0-2.d0*z/kval)
            sw=sw+wj
            sj=s(j1,j2,j3)
            sws=sws+wj*sj
            sws2=sws2+wj*sj*sj
            jj=jj+1
            wad(jj,thrednr)=wj
            sad(jj,thrednr)=sj
         END DO
         ni(i) = sw
         thn(i) = sws/sw
         if(sw.gt.minni) THEN
            ksi = sws2/sw
            low = sgi/1d1
            up = sgi*1d1
            call localmin(low,up,wad(1,thrednr),sad(1,thrednr),L,jj,
     1                    tol,maxit,work(1,thrednr),xmin,fmin)
            sigman(i)=xmin
         ELSE
            sigman(i)=sgi
         END IF
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni,sigman)
      RETURN
      END
      subroutine awslgaus(s,th,ni,sigma,mask,n1,n2,n3,ind,w,nw,
     1            minni,lambda,thn,sigman)
C
C  local variance estimation using (adaptive) weighted likelihood
C
C   Takes observed intensities in s and 
C     initial estimates of \sigma in sigma
C   perform adaptive smoothing on R^3
C   th containes previous estimates of E S
C   ni containes previous sum of weights
C   mask - logical mask (use if mask==TRUE)
C   n1,n2,n3 - dimensions
C   ind  - integer array dim (3,n) containing relative indices in xyz
C   w    - vector of corresponding location weights
C   nw   - number of positive weights (initial value 
C   
C   lambda   - kritical value for pairwise tests
C   thn      - new estimate sum_j w_a(j) S_j
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:3] are j1-i1,j2-i2 and j3-i3 respectively 
C   wad, sad - array for weights>0 and corresponding observed s
C
      implicit logical (a-z)
      integer n1,n2,n3,nw,ind(3,nw)
      logical mask(n1,n2,n3)
      real*8 s(n1,n2,n3),th(n1,n2,n3),ni(n1*n2*n3),thn(n1*n2*n3),
     1  sigman(n1*n2*n3),lambda,w(nw),sigma(n1,n2,n3),minni
      integer i1,i2,i3,j1,j2,j3,i,j,n
      real*8 z,sw,sws,sws2,sj,thi,wj,kval,sgi
      n = n1*n2*n3
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& FIRSTPRIVATE(minni,n1,n2,n3)
C$OMP& PRIVATE(i,j,i1,i2,i3,j1,j2,j3,z,sw,sws,sws2,thi,kval,wj,sj,sgi)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1         
         sgi=sigma(i1,i2,i3)
         sigman(i)=sgi
         thi = th(i1,i2,i3)
         thn(i) = thi
         if(.not.mask(i1,i2,i3)) CYCLE
         sw=0.d0
         sws=0.d0
         sws2=0.d0
C   thats the estimated standard deviation of s(i1,i2,i3)
         kval = lambda/ni(i)*sgi*sgi
         DO j=1,nw
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
            sj=s(j1,j2,j3)
            sws=sws+wj*sj
            sws2=sws2+wj*sj*sj
         END DO
         ni(i) = sw
         sws = sws / sw
         thn(i) = sws
         if(sw.gt.minni) THEN
            sws2 = (sws2/sw-sws*sws)*sw/(sw-1)
            if(sws2.gt.0.d0) THEN
               sigman(i)=sqrt(sws2)
            ELSE
               sigman(i)=sgi
            END IF
         ELSE
            sigman(i)=sgi
         END IF
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni,sigman)
      RETURN
      END
      subroutine awsvchi(y,th,ni,fns,mask,n1,n2,n3,ind,w,nw,lambda,
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
      real*8 y(n1,n2,n3),th(n1,n2,n3),ni(n1*n2*n3),thn(n1*n2*n3),
     1       sy(n1*n2*n3),lambda,w(nw),sigma,fns(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,i,j,n
      real*8 z,sw,sw2,swy,swy2,yj,thi,wj,kval,cw,fnsi
      n = n1*n2*n3
C  precompute values of lgamma(corrected df/2) in each voxel
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,j,i1,i2,i3,j1,j2,j3,z,sw,swy,sw2,swy2,thi,kval,
C$OMP& wj,yj,cw,fnsi)
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
         fnsi = fns(i1,i2,i3)
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
            z=z*z/(fnsi+fns(j1,j2,j3))
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
      subroutine awsadchi(y,th,ni,fns,mask,n1,n2,n3,ind,w,nw,lambda,
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
      real*8 y(n1,n2,n3),th(n1,n2,n3),ni(n1*n2*n3),thn(n1*n2*n3),
     1  fns(n1,n2,n3),sy(n1*n2*n3),lambda,w(nw),sigma,wad(nw,nthreds)
      integer i1,i2,i3,j1,j2,j3,i,j,n,thrednr
      real*8 z,sw,sw2,swy,swy2,yj,thi,wj,kval,cw,fnsi
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      n = n1*n2*n3
      thrednr = 1
C  precompute values of lgamma(corrected df/2) in each voxel
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,j,i1,i2,i3,j1,j2,j3,z,sw,swy,sw2,swy2,thi,kval,
C$OMP& wj,yj,cw,thrednr,fnsi)
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
         fnsi = fns(i1,i2,i3)
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
            z=z*z/(fnsi+fns(j1,j2,j3))
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
C  z  is an estimate of the standard deviation of y(i) by mean absolute deviation
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
      subroutine afmodevn(y,n1,n2,n3,mask,h,vext,sigma)
C
C   Aja-Fernandez Mode Vn (6)
C
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 y(n1,n2,n3),sigma(n1,n2,n3),h,vext(2)
      logical mask(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,ih1,ih2,ih3,ni
      real*8 m1,m2,z
      ih1=int(h)
      ih2=int(h*vext(1))
      ih3=int(h*vext(2))
      Do i1=1+ih1,n1-ih1
         Do i2=1+ih2,n2-ih2
            Do i3=1+ih3,n3-ih3
               if(mask(i1,i2,i3)) THEN
                  ni=0
                  m1=0.d0
                  m2=0.d0
                  DO j1=i1-ih1,i1+ih1
                     DO j2=i2-ih2,i2+ih2
                        DO j3=i3-ih3,i3+ih3
                           if(mask(j1,j2,j3)) THEN
                              z=y(j1,j2,j3)
                              m1=m1+z
                              m2=m2+z*z
                              ni=ni+1
                           ENDIF
                        END DO
                     END DO
                  END DO
                  m1=m1/ni
                  m2=m2/ni
                  z=m2-m1*m1
                  if(ni.gt.1) THEN
                     sigma(i1,i2,i3)=ni*z/(ni-1)
                  ELSE
                     sigma(i1,i2,i3)=0.d0
                  ENDIF
               ELSE
                  sigma(i1,i2,i3)=0.d0
               ENDIF
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine afmodem1(y,n1,n2,n3,mask,h,vext,sigma)
C
C   Aja-Fernandez Mode Vn (6)
C
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 y(n1,n2,n3),sigma(n1,n2,n3),h,vext(2)
      logical mask(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,ih1,ih2,ih3,ni
      real*8 m1
      ih1=int(h)
      ih2=int(h*vext(1))
      ih3=int(h*vext(2))
      Do i1=1+ih1,n1-ih1
         Do i2=1+ih2,n2-ih2
            Do i3=1+ih3,n3-ih3
               if(mask(i1,i2,i3)) THEN
                  ni=0
                  m1=0.d0
                  DO j1=i1-ih1,i1+ih1
                     DO j2=i2-ih2,i2+ih2
                        DO j3=i3-ih3,i3+ih3
                           if(mask(j1,j2,j3)) THEN
                              m1=m1+y(j1,j2,j3)
                              ni=ni+1
                           ENDIF
                        END DO
                     END DO
                  END DO
                  sigma(i1,i2,i3)=m1/ni
               ELSE
                  sigma(i1,i2,i3)=0.d0
               ENDIF
            END DO
         END DO
      END DO
      RETURN
      END
