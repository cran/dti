      subroutine besselq(x,n,fw)
      implicit logical (a-z) 
      integer n
      real*8 x(n),fw(n)
      external besseli
      real*8 besseli
      integer i
      DO i=1,n
         x(i)=i*1.d-2
C  avoid transform to zero for numerical stability
         fw(i)=besseli(x(i),1.d0,2.d0)/besseli(x(i),0.d0,2.d0)
      END DO
      RETURN
      END
      subroutine ricecorr(si,w,n,nb,sbind,ns0,niter,sw,th,s2,sigma2,
     1                    fw)
      implicit logical (a-z) 
      integer n,nb,niter(nb)
      real*8 si(nb,n)
      real*8 w(n),th(nb),s2(nb),sigma2,sw,fw(10000)
      logical sbind(nb)
      integer i,j,k,l,iter,ns0
      real*8 z,sth,th2,ss2,thk,sii,minsi
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
                  j=sii*z+1
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
      real*8 function Lhalf(x)
      real*8 x,mxh,bi(10)
      real*8 bessliex
      external bessliex
      mxh = -0.5d0*x
C      Lhalf = (1.d0-x)*besseli(mxh,0.d0,2.d0)-x*besseli(mxh,1.d0,2.d0)
      Lhalf = (1.d0-x)*bessliex(mxh,0.d0,2.d0,bi)-
     1               x*bessliex(mxh,1.d0,2.d0,bi)
      Return
      End
      real*8 function erice(x)
      real*8 x
      real*8 Lhalf
      external Lhalf
      if(x.gt.5d1) THEN
         erice = x
      ELSE
         erice = 1.253314137d0*Lhalf(-0.5d0*x*x)
      END IF
      RETURN
      END
      subroutine rcstep(x,n)
      implicit logical (a-z)
      integer n
      real*8 x(n)
      integer i,counter
      real*8 xt,xta,yt,erice
      external erice
C ####################################################
C  Attempt to parallelize causes errors
C  " ****snapping into wrong generation "
C  probably coming from http://svn.r-project.org/R/trunk/src/main/memory.c
C ####################################################
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(x,n)
C$OMP& PRIVATE(i,counter,xt,xta,yt)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         counter=0
         xta = 0
         xt  = x(i)
         yt  = xt
         DO while (abs(xt-xta).gt.4d-4.and.(counter.lt.50))
            counter = counter + 1
            xta = xt
            if(xt.le.5d1) THEN
               xt = max(0.d0,xt+yt-erice(xt))
            ELSE
               xt = yt
            END IF
C            if(counter.eq.100) THEN
C               call intpr("maxit in i",10,i,1)
C               call dblepr("yt",2,yt,1)
C               call dblepr("xt",2,xt,1)
C               call dblepr("xta",3,xta,1)
C            END IF
         END DO
         x(i) = xt
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(x)
      RETURN
      END
      
      
