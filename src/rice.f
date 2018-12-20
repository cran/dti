      subroutine ricecorr(si,w,n,nb,sbind,ns0,niter,sw,th,s2,sigma2,
     1                    fw)
      implicit none
      integer n,nb,niter(nb)
      double precision si(nb,n)
      double precision w(n),th(nb),s2(nb),sigma2,sw,fw(10000)
      logical sbind(nb)
      integer i,j,k,l,iter,ns0
      double precision z,sth,th2,ss2,thk,sii,minsi
      iter=1
      DO k=1,nb
         iter=max(iter,niter(k))
      END DO
      DO l=1,iter
         Do k=1,nb
            if(l.le.niter(k)) THEN
               sth=0.d0
               ss2=0.d0
               thk=th(k)
               z=thk/sigma2/1.d-2
               th2=thk*thk
               minsi=65535
               DO i=1,n
                  sii=si(k,i)
                  minsi=min(sii,minsi)
                  j=int(sii*z+1)
                  if(j.gt.10000) THEN
                     sth=sth+sii*w(i)
                     ss2=ss2+w(i)*((sii*sii+th2)*0.5d0-sii*thk)
                  ELSE
                     sth=sth+fw(j)*sii*w(i)
                     ss2=ss2+w(i)*((sii*sii+th2)*0.5d0-fw(j)*sii*thk)
                  END IF
               END DO
               z=minsi
               th(k)=max(z/3.d0,sth/sw)
C  just avoid extreme corrections caused by overestimated variances
               s2(k)=ss2/sw
            END IF
         END DO
         ss2=0.d0
         DO k=1,nb
            if(sbind(k)) ss2=ss2+s2(k)
         END DO
         sigma2=ss2/ns0
      END DO
      RETURN
      END
      
